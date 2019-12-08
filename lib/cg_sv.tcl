proc sv_job {args} {
	global appdir
	upvar job_logdir job_logdir
	set method manta
	set refseq {}
	set distrreg chr
	set opts {}
	set split 1
	set deps {}
	set threads 2
	set cleanup 1
	set regmincoverage 3
	set opts {}
	set cmdline [list cg sv $args]
	set resultfile {}
	cg_options sv args {
		-method {
			set method $value
		}
		-refseq {
			set refseq $value
		}
		-regionfile {
			set regionfile $value
		}
		-regmincoverage {
			set regmincoverage $value
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
			if {$value in {1 chr chromosome 0}} {
				set distrreg $value
			} elseif {[isint $value]} {
				set distrreg $value
			} elseif {[file exists $value]} {
				set distrreg [file_absolute $value]
			} else {
				error "unknown value $value for -distrreg, must be a (region) file or one of: chr or 1, 0"
			}
		}
		default {
			lappend opts $key $value
		}
	} {bamfile resultfile} 1 2
	set bamfile [file_absolute $bamfile]
	set refseq [refseq $refseq]
	set destdir [file dir $bamfile]
	# logfile
	job_logfile $destdir/sv_{$method}_[file tail $bamfile] $destdir $cmdline \
		{*}[versions bwa bowtie2 samtools gatk picard java gnusort8 zst os]
	# check if regionfile is supported
	catch {sv_${method}_job} temp
	if {[regexp {with options:(.*)} $temp temp temp]} {
		if {![inlist [split $temp ,] -regionfile]} {
			putslog "sv_$method does not support -regionfile, so cannot be run distributed, -distrreg and -regionfile ignored"
			return [sv_${method}_job {*}$opts \
				-split $split -threads $threads -cleanup $cleanup \
				-refseq $refseq $bamfile $resultfile]
		}
	}
	# get or create regionfile
	if {[info exists regionfile]} {
		set regionfile [file_absolute $regionfile]
	} else {
		set regionfile [bam2reg_job -mincoverage $regmincoverage $bamfile]
	}
	# run
	if {$distrreg in {0 {}}} {
		sv_${method}_job {*}$opts -regionfile $regionfile \
			-split $split -threads $threads -cleanup $cleanup	-refseq $refseq $bamfile $resultfile
	} else {
		# check what the resultfiles are for the method
		set resultfiles [sv_${method}_job -resultfiles 1 {*}$opts \
			-regionfile $regionfile \
			-split $split	-refseq $refseq $bamfile $resultfile]
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
		set regions [distrreg_regs $distrreg $refseq]
		set basename [gzroot [file tail $varallfile]]
		set regfiles {}
		foreach region $regions {
			lappend regfiles $indexdir/$basename-$region.bed
		}
		job [gzroot $varallfile]-distrreg-beds {*}$skips -deps {
			$regionfile
		} -targets $regfiles -vars {
			regionfile regions appdir basename indexdir
		} -code {
			set header [cg select -h $regionfile]
			set poss [tsv_basicfields $header 3]
			set header [list_sub $header $poss]
			# puts "cg select -f \'$header\' $regionfile | $appdir/bin/distrreg $indexdir/$basename- \'$regions\' 0 1 2"
			cg select -f $header $regionfile | $appdir/bin/distrreg $indexdir/$basename- .bed 0 $regions 0 1 2
		}
		set todo {}
		# Produce variant calls
		set ibam $indexdir/[file tail $bamfile]
		set indexext [indexext $bamfile]
		mklink $bamfile $ibam
		mklink $bamfile.$indexext $ibam.$indexext
		defcompressionlevel 1
		foreach region $regions regfile $regfiles {
			lappend todo [var_${method}_job {*}$opts {*}$skips -rootname $root-$region -regionfile $regfile \
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
					cg cat -c f {*}[lsort -dict [jobglob {*}$list]] {*}[compresspipe $target] > $target.temp
					file rename $target.temp $target
					cg_zindex $target
				}
				foreach file $list {
					file delete $file
					if {[file extension $file] eq ".lz4"} {file delete $file.lz4i}
					if {[file extension $file] eq ".zst"} {file delete $file.zsti}
					file delete -force [gzroot $file].index
					file delete [gzroot $file].analysisinfo
				}
			}
			incr pos
		}
		cleanup_job cleanup-sv_${method}_[file tail $bamfile] $indexdir $resultfiles
		cd $keeppwd
		return $resultfiles
	}
}

proc cg_sv {args} {
	set args [job_init {*}$args]
	sv_job {*}$args
	job_wait
}
