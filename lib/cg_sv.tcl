proc sv_job {args} {
	global appdir
	upvar job_logdir job_logdir
	set method manta
	set refseq {}
	set distrreg 0
	set opts {}
	set split 1
	set deps {}
	set regionfile {}
	set threads 2
	set cleanup 1
	set regmincoverage 3
	set opts {}
	set preset {}
	set cmdline [list cg sv $args]
	set resultfile {}
	cg_options sv args {
		-method {
			set method $value
		}
		-preset {
			set preset $value
		}
		-refseq {
			set refseq $value
		}
		-regionfile {
			set regionfile [file_absolute $value]
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
			set distrreg [distrreg_checkvalue $value]
		}
		default {
			lappend opts $key $value
		}
	} {bamfile resultfile} 1 2
	set bamfile [file_absolute $bamfile]
	set refseq [refseq $refseq]
	set destdir [file dir $bamfile]
	# logfile
	job_logfile $destdir/sv_${method}_[file tail $bamfile] $destdir $cmdline \
		{*}[versions bwa bowtie2 samtools gatk picard java gnusort8 zst os]
	# check if regionfile or region is supported
	set cmdopts {}
	catch {sv_${method}_job} temp
	set supportsregionfile 0
	set supportsregion 0
	if {[regexp {with options:(.*)} $temp temp methodoptions]} {
		set methodoptions [split $methodoptions ,]
		if {[inlist $methodoptions -regionfile]} {
			lappend cmdopts -regionfile $regionfile
			set supportsregionfile 1
		}
		if {[inlist $methodoptions -region]} {
			set supportsregion 1
		}
	}
	catch {
		# see if method wants to change distrreg from requested
		set distrreg [sv_${method}_distrreg $distrreg]
	}
	if {!$supportsregionfile && !$supportsregion} {
		if {$distrreg ni {0 {}}} {
			putslog "sv_$method does not support -regionfile, so cannot be run distributed, -distrreg and -regionfile ignored"
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
		sv_${method}_job {*}$opts	-preset $preset {*}$cmdopts \
			-split $split -threads $threads -cleanup $cleanup	\
			-refseq $refseq $bamfile $resultfile
	} else {
		# check what the resultfiles are for the method
		set resultfiles [sv_${method}_job -resultfiles 1 {*}$opts	-preset $preset \
			-split $split	-refseq $refseq $bamfile $resultfile]
		set skips [list -skip $resultfiles]
		foreach {resultfile resultanalysisfile resultvcffile} $resultfiles break
		set file [file tail $bamfile]
		set root [file_rootname $resultfile]
		# start
		## Create sequencing region files
		set keeppwd [pwd]
		cd $destdir
		set workdir [workdir $resultfile]
		set regions [list_remove [distrreg_regs $distrreg $refseq] unaligned]
		set basename [file_root [file tail $resultfile]]
		set regfiles {}
		if {$supportsregionfile} {
			foreach region $regions {
				lappend regfiles $workdir/$basename-$region.bed
			}
			if {[info exists regionfile]} {
				job [gzroot $resultfile]-distrreg-beds {*}$skips -deps {
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
			} else {
				# if regionfile is supported, but not defined, 
				# make a regionfile for each region inregions and use for distrreg
				foreach region $regions {
					distrreg_reg2bed $workdir/$basename-$region.bed $region $refseq
				}
			}
		} else {
			# we'll just run per region
		}
		# Produce SV calls
		set ibam $workdir/[file tail $bamfile]
		set indexext [indexext $bamfile]
		mklink $bamfile $ibam
		mklink $bamfile.$indexext $ibam.$indexext
		defcompressionlevel 1
		set todo {}
		set sample [file_rootname $resultfile]
		if {$supportsregionfile} {
			foreach region $regions regfile $regfiles {
				set target $workdir/sv-$root-$region.tsv.zst
				lappend todo [sv_${method}_job {*}$opts {*}$skips	\
					-split $split -threads $threads -cleanup $cleanup \
					-preset $preset \
					-regionfile $regfile \
					-refseq $refseq \
					-sample $sample \
					$ibam $target]
			}
		} else {
			# run per region
			foreach region $regions {
				set target $workdir/sv-$root-$region.tsv.zst
				lappend todo [sv_${method}_job {*}$opts {*}$skips \
					-split $split -threads $threads -cleanup $cleanup \
					-preset $preset \
					-region $region \
					-refseq $refseq \
					-sample $sample \
					$bamfile $target]

			}
		}
		defcompressionlevel 9
		# concatenate results
		set pos 0
		set cleanupfiles $regfiles
		foreach resultfile $resultfiles {
			set list [list_subindex $todo $pos]
			set deps $list
			job sv_combineresults-[file tail $resultfile] {*}$skips -deps $list -rmtargets $list -targets {
				$resultfile
			} -vars {
				analysisinfo list method regfile distrreg sample
			} -code {
				if {![auto_load sv_${method}_sortdistrreg]} {
					set sort 0
					set sortpipe {}
				} else {
					set sort [sv_${method}_sortdistrreg]
					set sortpipe {| cg select -s - }
				}
				set ext [file extension [gzroot $target]]
				set analysisinfo [analysisinfo_file $dep]
				if {[file exists $analysisinfo] && $ext ni ".analysisinfo"} {
					analysisinfo_copy $analysisinfo [analysisinfo_file $target] [list \
						svcaller_region [file tail $regfile] \
						svcaller_distrreg [file tail $distrreg] \
						sample $sample]
					exec touch [analysisinfo_file $target]
				}
				if {$ext in ".vcf .gvcf"} {
					set sample [file_sample $target]
					cg vcfcat -i 0 -s $sort -o $target -sample $sample {*}[bsort [jobglob {*}$list]]
				} elseif {$ext in ".analysisinfo"} {
					analysisinfo_copy $dep $target [list \
						svcaller_region [file tail $regfile] \
						svcaller_distrreg [file tail $distrreg] \
						sample $sample]
				} else {
					cg cat -c m -m 1 -sample $sample {*}[bsort [jobglob {*}$list]] {*}$sortpipe {*}[compresspipe $target] > $target.temp
					file rename -force $target.temp $target
					cg_zindex $target
				}
			}
			foreach file $list {
				lappend cleanupfiles $file
				if {[file extension $file] eq ".lz4"} {lappend cleanupfiles $file.lz4i}
				if {[file extension $file] eq ".zst"} {lappend cleanupfiles $file.zsti}
				lappend cleanupfiles [gzroot $file].temp
				lappend cleanupfiles [analysisinfo_file $file]
			}
			incr pos
		}
		if {[llength $regfiles]} {
			cleanup_job cleanup-sv_${method}_[file tail $bamfile] [list {*}$cleanupfiles] $resultfiles
		}
		cd $keeppwd
		return $resultfiles
	}
}

proc cg_sv {args} {
	set args [job_init {*}$args]
	sv_job {*}$args
	job_wait
}
