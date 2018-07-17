proc sv_job {args} {
	global appdir
	upvar job_logdir job_logdir
	set method cg
	set distrreg chr
	set pre ""
	set opts {}
	set split 1
	set deps {}
	set threads 2
	set cleanup 1
	set regmincoverage 3
	set opts {}
	set cmdline [list cg sv $args]
	cg_options sv args {
		-method {
			set method $value
		}
		-regionfile {
			set regionfile $value
		}
		-regmincoverage {
			set regmincoverage $value
		}
		-pre {
			set pre $value
		}
		-split {
			set split $value
		}
		-threads - -t {
			set threads $value
		}
		-cleanup {
			set cleanup $value
		}
		-distrreg {
			if {$value ni {1 chr chromosome 0}} {
				error "unknown value $value for -distrreg, must be one of: chr or 1, 0"
			}
			set distrreg $value
		}
		default {
			lappend opts $key $value
		}
	} {bamfile refseq}
	set bamfile [file_absolute $bamfile]
	set refseq [file_absolute $refseq]
	set destdir [file dir $bamfile]
	if {[info exists regionfile]} {
		set regionfile [file_absolute $regionfile]
	} else {
		set regionfile [bam2reg_job -mincoverage $regmincoverage $bamfile]
	}
	# logfile
	job_logfile $destdir/sv_{$method}_[file tail $bamfile] $destdir $cmdline \
		{*}[versions bwa bowtie2 samtools gatk picard java gnusort8 lz4 os]
	# check if regionfile is supported
	catch {sv_${method}_job} temp
	if {[regexp {with options:(.*)} $temp temp temp]} {
		if {![inlist [split $temp ,] -regionfile]} {
			putslog "sv_$method does not support distrreg"
			return [sv_${method}_job {*}$opts -pre $pre \
				-split $split -threads $threads -cleanup $cleanup $bamfile $refseq]
		}
	}
	# run
	if {$distrreg in {0 {}}} {
		sv_${method}_job {*}$opts -regionfile $regionfile -pre $pre \
			-split $split -threads $threads -cleanup $cleanup $bamfile $refseq
	} else {
		# check what the resultfiles are for the method
		set resultfiles [sv_${method}_job -resultfiles 1 {*}$opts \
			-regionfile $regionfile -pre $pre \
			-split $split $bamfile $refseq]
		set skips [list -skip $resultfiles]
		# if {[jobtargetexists $resultfiles [list $refseq $bamfile $regionfile]]} return
		foreach {svfile} $resultfiles break
		set file [file tail $bamfile]
		set root [file_rootname $varfile]
		# start
		## Create sequencing region files
		set keeppwd [pwd]
		cd $destdir
		set indexdir [gzroot $varallfile].index
		file mkdir $indexdir
		set chromosomes [distrreg_regs $distrreg $refseq]
		set basename [gzroot [file tail $varallfile]]
		set regfiles {}
		foreach chromosome $chromosomes {
			lappend regfiles $indexdir/$basename-$chromosome.bed
		}
		job [gzroot $varallfile]-distrreg-beds {*}$skips -deps {
			$regionfile
		} -targets $regfiles -vars {
			regionfile chromosomes appdir basename indexdir
		} -code {
			cg select -f {chromosome begin end} $regionfile | $appdir/bin/distr2chr $indexdir/$basename-
			foreach chromosome $chromosomes {
				if {[file exists $indexdir/$basename-$chromosome]} {
					file rename -force $indexdir/$basename-$chromosome $indexdir/$basename-$chromosome.bed
				} else {
					file_write $indexdir/$basename-$chromosome.bed ""
				}
			}
			file delete $indexdir/$basename-chromosome
		}
		set todo {}
		# Produce variant calls
		set ibam $indexdir/[file tail $bamfile]
		mklink $bamfile $ibam
		mklink $bamfile.bai $ibam.bai
		defcompressionlevel 1
		foreach chromosome $chromosomes regfile $regfiles {
			lappend todo [var_${method}_job {*}$opts {*}$skips -rootname $root-$chromosome -regionfile $regfile \
				-split $split -threads $threads -cleanup $cleanup $ibam $refseq]
		}
		defcompressionlevel 9
		# concatenate results
		set pos 0
		foreach resultfile $resultfiles {
			set list [list_subindex $todo $pos]
			set deps $list
			set ainfolist {}
			foreach el $list {
				set analysisinfo [lindex [jobglob [gzroot $el].analysisinfo]]
				if {$analysisinfo ne ""} {
					lappend ainfolist $analysisinfo
				}
				
			}
			set analysisinfo [lindex $ainfolist 0]
			lappend deps {*}$ainfolist
			job $resultfile  {*}$skips -deps $list -rmtargets $list -targets {
				$resultfile
			} -vars {
				analysisinfo list
			} -code {
				if {[llength $analysisinfo]} {
					file copy -force $analysisinfo [gzroot $target].analysisinfo
				}
				if {[file extension [gzroot $target]] in ".vcf .gvcf"} {
					cg vcfcat -i 1 -o $target {*}[jobglob {*}$list]
				} else {
					cg cat -c f {*}[jobglob {*}$list] | cg lz4 -c 9 > $target.temp
					file rename $target.temp $target
					cg_lz4index $target
				}
				foreach file $list {
					file delete $file
					if {[file extension $file] eq ".lz4"} {file delete $file.lz4i}
					file delete -force [gzroot $file].index
					file delete [gzroot $file].analysisinfo
				}
			}
			incr pos
		}
		cleanup_job cleanup-var_distrreg_[file tail $bamfile] $indexdir $resultfiles
		cd $keeppwd
		return $resultfiles
	}
}

proc cg_sv {args} {
	set args [job_init {*}$args]
	sv_job {*}$args
	job_wait
}
