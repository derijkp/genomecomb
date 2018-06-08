proc var_distrreg_regs {regfile refseq} {
	if {$regfile in "chr chromosome 1"} {
		return [exec cut -f 1 $refseq.fai]
	}
	set f [gzopen $regfile]
	set header [tsv_open $f]
	set poss [tsv_basicfields $header 3]
	set result {}
	while {[gets $f line] != -1} {
		lappend result [join [list_sub [split $line \t] $poss] _]
	}
	return $result
}

proc var_distrreg_job {args} {
	global appdir
	upvar job_logdir job_logdir
	set method gatk
	set distrreg chr
	set pre ""
	set opts {}
	set split 1
	set deps {}
	set regionfile {}
	set threads 2
	set cleanup 1
	set regmincoverage 3
	set opts {}
	set cmdline [list cg var_distrreg $args]
	cg_options var_distrreg args {
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
	if {$regionfile ne ""} {
		set regionfile [file_absolute $regionfile]
	} else {
		set regionfile [bam2reg_job -mincoverage $regmincoverage $bamfile]
	}
	# logfile
	job_logfile $destdir/var_distrreg_[file tail $bamfile] $destdir $cmdline \
		{*}[versions bwa bowtie2 samtools gatk picard java gnusort8 lz4 os]
	# run
	if {$distrreg in {0 {}}} {
		var_${method}_job {*}$opts -regionfile $regionfile -pre $pre \
			-split $split -threads $threads -cleanup $cleanup $bamfile $refseq
	} else {
		# check what the resultfiles are for the method
		set resultfiles [var_${method}_job -resultfiles 1 {*}$opts \
			-regionfile $regionfile -pre $pre \
			-split $split $bamfile $refseq]
		if {[llength [jobglob $resultfiles]] == [llength $resultfiles]} {
			putslog "var_distrreg: resultfiles ($resultfiles) already made, skipping"
			return $resultfiles
		}
		foreach {varfile sregfile varallfile} $resultfiles break
		set file [file tail $bamfile]
		set root [file_rootname $varfile]
		# start
		## Create sequencing region files
		set keeppwd [pwd]
		cd $destdir
		set indexdir [gzroot $varallfile].index
		file mkdir $indexdir
		set chromosomes [var_distrreg_regs $distrreg $refseq]
		set basename [gzroot [file tail $varallfile]]
		set regfiles {}
		foreach chromosome $chromosomes {
			lappend regfiles $indexdir/$basename-$chromosome.bed
		}
		job [gzroot $varallfile]-distrreg-beds -deps {
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
			lappend todo [var_${method}_job {*}$opts -rootname $root-$chromosome -regionfile $regfile \
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
			job $resultfile -deps $list -rmtargets $list -targets {
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
				foreach file $deps {
					file delete $file
					if {[file extension $file] eq ".lz4"} {file delete $file.lz4i}
					file delete [gzroot $file].index
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

proc cg_var_distrreg {args} {
	set args [job_init {*}$args]
	var_distrreg_job {*}$args
	job_wait
}
