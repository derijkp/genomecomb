proc var_job {args} {
	upvar job_logdir job_logdir
	set cmdline "[list cd [pwd]] \; [list cg var {*}$args]"
	global appdir
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
	cg_options var args {
		-method {
			set method $value
		}
		-regionfile {
			set regionfile $value
		}
		-regmincoverage {
			set regmincoverage $value
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
		default {
			lappend opts $key $value
		}
	} {bamfile refseq} 2 2
	set bamfile [file_absolute $bamfile]
	set refseq [file_absolute $refseq]
	set destdir [file dir $bamfile]
	if {$regionfile ne ""} {
		set regionfile [file_absolute $regionfile]
	} else {
		set regionfile [bam2reg_job -mincoverage $regmincoverage $bamfile]
	}
	# logfile
	job_logfile $destdir/var_${method}_[file tail $bamfile] $destdir $cmdline \
		{*}[versions bwa bowtie2 samtools gatk picard java gnusort8 zst os]
	# run
	if {$distrreg in {0 {}}} {
		var_${method}_job {*}$opts -regionfile $regionfile -pre $pre \
			-split $split -threads $threads -cleanup $cleanup $bamfile $refseq
	} else {
		# check what the resultfiles are for the method
		set resultfiles [var_${method}_job -resultfiles 1 {*}$opts \
			-regionfile $regionfile -pre $pre \
			-split $split $bamfile $refseq]
		set skips [list -skip $resultfiles]
		# if {[jobtargetexists $resultfiles [list $refseq $bamfile $regionfile]]} return
		foreach {varfile sregfile varallfile} $resultfiles break
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
		mklink $bamfile $ibam
		mklink $bamfile.bai $ibam.bai
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
			job $resultfile  {*}$skips -deps $deps -rmtargets $list -targets {
				$resultfile
			} -vars {
				list
			} -code {
				set analysisinfo [gzroot $dep].analysisinfo
				if {[file exists $analysisinfo]} {
					file copy -force $analysisinfo [gzroot $target].analysisinfo
					exec touch [gzroot $target].analysisinfo
				}
				if {[file extension [gzroot $target]] in ".vcf .gvcf"} {
					cg vcfcat -i 1 -o $target {*}[jobglob {*}$list]
				} else {
					cg cat -c f {*}[jobglob {*}$list] {*}[compresspipe $target] > $target.temp
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
		cleanup_job cleanup-var_${method}_[file tail $bamfile] $indexdir $resultfiles
		cd $keeppwd
		return $resultfiles
	}
}

proc cg_var {args} {
	set args [job_init {*}$args]
	var_job {*}$args
	job_wait
}
