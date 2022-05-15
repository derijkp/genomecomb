proc var_job {args} {
	upvar job_logdir job_logdir
	global appdir
	set cmdline "[list cd [pwd]] \; [list cg var {*}$args]"
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
	set datatype {}
	set cmdopts {}
	set var_opts {}
	set resultfile {}
	set hap_bam 0
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
			if {$value eq "regionfile"} {
				set distrreg regionfile
			} else {
				set distrreg [distrreg_checkvalue $value]
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
		-datatype {
			set datatype $value
		}
		-opts {
			set opts $value
		}
		-hap_bam {
			set hap_bam [true $value]
		}
		default {
			lappend var_opts $key $value
		}
	} {bamfile refseq resultfile} 2 3
	set bamfile [file_absolute $bamfile]
	set refseq [file_absolute $refseq]
	if {$resultfile eq ""} {
		set resultfile [file dir $bamfile]/${pre}var-${method}-[file_rootname $bamfile].tsv.zst
	}
	set destdir [file dir $resultfile]
	# logfile
	set tools {gatk picard java gnusort8 zst os}
	catch {lappend tools {*}[var_${method}_tools]}
	set tools [list_remdup $tools]
	job_logfile $destdir/var_${method}_[file tail $bamfile] $destdir $cmdline \
		{*}[versions {*}$tools]
	# check if regionfile is supported
	set cmdopts {}
	catch {var_${method}_job} temp
	set supportsregionfile 0
	set supportsregion 0
	if {[regexp {with options: *([^\n]*)} $temp temp methodoptions]} {
		set methodoptions [split $methodoptions ,]
		if {[inlist $methodoptions -regionfile]} {
			lappend cmdopts -regionfile $regionfile
			set supportsregionfile 1
		}
		if {[inlist $methodoptions -region]} {
			set supportsregion 1
		}
		if {![inlist $methodoptions -hap_bam]} {
			set hap_bam 0
		}
	} else {
		set hap_bam 0
	}
	if {$hap_bam} {
		# limit regdistr to chr as smaller distribution than chr can cause errors when making a haplotyped bam:
		# alignments can stick out of the regions, or even cover multiple regions
		# -> they could get in multiple times, or get incorrectly sorted
		if {$distrreg ni "0 chr"} {set distrreg chr}
		# give option only when supported and needed via var_opts
		lappend var_opts -hap_bam 1
	}
	catch {
		# see if method wants to change distrreg from requested
		set distrreg [var_${method}_distrreg $distrreg]
	}
	if {!$supportsregionfile && !$supportsregion} {
		if {$distrreg ni {0 {}}} {
			putslog "var_$method does not support -regionfile or -region, so cannot be run distributed, -distrreg and -regionfile ignored"
		}
		set distrreg 0
	}
	if {$supportsregionfile} {
		if {$regionfile ne ""} {
			set regionfile [file_absolute $regionfile]
		} else {
			set regionfile [bam2reg_job -mincoverage $regmincoverage -distrreg $distrreg -refseq $refseq $bamfile]
		}
	}
	# run
	if {$distrreg in {0 {}}} {
		var_${method}_job {*}$var_opts -opts $opts {*}$cmdopts -datatype $datatype -pre $pre \
			-split $split -threads $threads -cleanup $cleanup $bamfile $refseq $resultfile
	} else {
		# check what the resultfiles are for the method
		set resultfiles [var_${method}_job -resultfiles 1 {*}$var_opts -opts $opts \
			-pre $pre \
			-datatype $datatype \
			-split $split $bamfile $refseq $resultfile]
		set skips [list -skip [list_remove  $resultfiles {}]]
		# if {[jobtargetexists $resultfiles [list $refseq $bamfile $regionfile]]} return
		foreach {varfile sregfile varallfile vcffile} $resultfiles break
		set file [file tail $bamfile]
		set root [file_rootname $varfile]
		# start
		## Create sequencing region files
		set workdir [workdir $varfile]
		file mkdir $workdir
		if {$distrreg eq "regionfile"} {
			if {$regionfile eq ""} {
				error "used -distrreg regionfile and -regionfile is not given"
			}
			if {![file exists $regionfile]} {
				error "-regionfile $regionfile does not exist"
			}
			set regions {}
			set f [gzopen $regionfile]
			set header [tsv_open $f]
			set poss [tsv_basicfields $header 3]
			while {[gets $f line] != -1} {
				foreach {chr b e} [list_sub [split $line \t] $poss] break
				lappend regions ${chr}-${b}-${e}
			}
			close $f
			# make sure it is going to work using -region only
			set supportsregionfile 0
		} else {
			# at some point I thought of supporting -regionfile here for methods that do not support it using -region
			# however, forcing working on the typically small regions in -regionfile is a bad idea 
			# as several methods do not support -regionfile because working on such small regions would give 
			# bad results (longshot and medaka use haplotype info for improving varcall)
			set regions [list_remove [distrreg_regs $distrreg $refseq] unaligned]
			set basename [gzroot [file tail $varallfile]]
			if {$supportsregionfile} {
				set regfiles {}
				foreach region $regions {
					lappend regfiles $workdir/$basename-$region.bed
				}
				job [gzroot $varallfile]-distrreg-beds {*}$skips -deps {
					$regionfile
				} -targets $regfiles -vars {
					regionfile regions appdir basename workdir
				} -code {
					set header [cg select -h $regionfile]
					set poss [tsv_basicfields $header 3]
					set header [list_sub $header $poss]
					# puts "cg select -f \'$header\' $regionfile | $appdir/bin/distrreg $workdir/$basename- \'$regions\' 0 1 2 1 \#"
					cg select -f $header $regionfile | $appdir/bin/distrreg $workdir/$basename- .bed 0 $regions 0 1 2 1 \#
				}
			}
		}
		set todo {}
		# Produce variant calls
		defcompressionlevel 1
		if {$supportsregionfile && !$hap_bam} {
			foreach region $regions regfile $regfiles {
				lappend todo [var_${method}_job {*}$var_opts -opts $opts {*}$skips \
					-datatype $datatype \
					-regionfile $regfile \
					-split $split -threads $threads -cleanup $cleanup $bamfile $refseq $workdir/var-$root-$region.tsv.zst]
			}
		} else {
			# if making a new bam file (currently only longshot), we need to process the unmapped reads as well
			if {$hap_bam} {
				lappend regions unmapped
			}
			# run per region
			foreach region $regions {
				lappend todo [var_${method}_job {*}$var_opts -opts $opts {*}$skips \
					-datatype $datatype \
					-region $region \
					-split $split -threads $threads -cleanup $cleanup \
					$bamfile $refseq $workdir/var-$root-$region.tsv.zst]
			}
		}
		defcompressionlevel 9
		# concatenate results
		set pos -1
		set cleanupfiles {}
		foreach resultfile $resultfiles {
			incr pos
			if {$resultfile eq ""} continue
			set list [bsort [list_subindex $todo $pos]]
			set deps $list
			job var_combineresults-[file tail $resultfile] {*}$skips -deps $deps -rmtargets $list -targets {
				$resultfile
			} -vars {
				list regionfile method bamfile
			} -code {
				set analysisinfo [analysisinfo_file $dep]
				if {[file exists $analysisinfo]} {
					set root $method-[file_rootname [file tail $bamfile]]
					analysisinfo_copy $analysisinfo [analysisinfo_file $target] \
						[list varcaller_region [file tail $regionfile] \
						sample $root]
					exec touch [analysisinfo_file $target]
				}
				# using bsort because
				# files are xxx-100 -> would sort reverse of what we want because of the -
				if {[file extension [gzroot $target]] in ".vcf .gvcf"} {
					cg vcfcat -i 1 -o $target {*}[bsort [jobglob {*}$list]]
				} elseif {[file extension [gzroot $target]] in ".bam"} {
					cg sam_catmerge \
						-sort nosort \
						-index 1 \
						-deletesams 1 \
						$target {*}[bsort [jobglob {*}$list]]
				} elseif {[lindex [split [file tail $target] -_] 0] eq "sreg"} {
					cg cat -m 1 -c f {*}[bsort [jobglob {*}$list]] | cg regjoin {*}[compresspipe $target] > $target.temp
					file rename -force $target.temp $target
					cg_zindex $target
				} else {
					cg cat -m 1 -c f {*}[bsort [jobglob {*}$list]] {*}[compresspipe $target] > $target.temp
					file rename -force $target.temp $target
					cg_zindex $target
				}
			}
			lappend cleanupfiles {*}$list
		}
		foreach bam $bams {
			bam_index_job $bam
		}
		cleanup_job -forcedirs 1 -delassociated 1 cleanup-var_${method}_[file tail $varfile] $cleanupfiles $resultfiles
		return $resultfiles
	}
}

proc cg_var {args} {
	set args [job_init {*}$args]
	var_job {*}$args
	job_wait
}
