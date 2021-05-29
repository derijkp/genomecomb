proc sv_job {args} {
	global appdir
	upvar job_logdir job_logdir
	set method manta
	set refseq {}
	set distrreg 0
	set opts {}
	set split 1
	set deps {}
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
	# check if regionfile is supported
	catch {sv_${method}_job} temp
	if {[regexp {with options:(.*)} $temp temp temp] && ![inlist [split $temp ,] -regionfile]} {
		set distrreg 0
	}
	if {
		($distrreg in {0 {}} && ![info exists regionfile])
		|| ([catch {sv_${method}_job} temp] && [regexp {with options:(.*)} $temp temp temp] && ![inlist [split $temp ,] -regionfile])
	} {
		if {[info exists regionfile] || $distrreg ni {0 {}}} {
			putslog "sv_$method does not support -regionfile, so cannot be run distributed, -distrreg and -regionfile ignored"
		}
		return [sv_${method}_job {*}$opts	-preset $preset \
			-split $split -threads $threads -cleanup $cleanup \
			-refseq $refseq $bamfile $resultfile]
	}			
	# run
	if {$distrreg in {0 {}}} {
		sv_${method}_job {*}$opts	-preset $preset -regionfile $regionfile \
			-split $split -threads $threads -cleanup $cleanup	-refseq $refseq $bamfile $resultfile
	} else {
		# check what the resultfiles are for the method
		set resultfiles [sv_${method}_job -resultfiles 1 {*}$opts	-preset $preset \
			-split $split	-refseq $refseq $bamfile $resultfile]
		set skips [list -skip $resultfiles]
		foreach {resultfile resultanalysisfile resultvcffile} $resultfiles break
		set file [file tail $bamfile]
		set root [file_rootname $file]
		# start
		## Create sequencing region files
		set keeppwd [pwd]
		cd $destdir
		set workdir [workdir $resultfile]
		set regions [distrreg_regs $distrreg $refseq]
		set basename [file_root [file tail $resultfile]]
		set regfiles {}
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
			foreach region $regions {
				distrreg_reg2bed $workdir/$basename-$region.bed $region $refseq
			}
		}
		# Produce variant calls
		set ibam $workdir/[file tail $bamfile]
		set indexext [indexext $bamfile]
		mklink $bamfile $ibam
		mklink $bamfile.$indexext $ibam.$indexext
		defcompressionlevel 1
		set todo {}
		foreach region $regions regfile $regfiles {
			set target [file root $ibam]-$region.tsv.zst
			lappend todo [sv_${method}_job {*}$opts {*}$skips	-preset $preset -regionfile $regfile \
				-split $split -threads $threads -cleanup $cleanup -refseq $refseq $ibam $target]
		}
		defcompressionlevel 9
		# concatenate results
		set pos 0
		foreach resultfile $resultfiles {
			set list [list_subindex $todo $pos]
			set deps $list
			set ainfolist {}
			foreach el $list {
				set analysisinfo [lindex [jobglob [analysisinfo_file $el]]]
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
					file_copy $analysisinfo [analysisinfo_file $target]
				}
				if {[file extension [gzroot $target]] in ".vcf .gvcf"} {
					cg vcfcat -i 1 -o $target {*}[jobglob {*}$list]
				} else {
					cg cat -c f {*}[bsort [jobglob {*}$list]] {*}[compresspipe $target] > $target.temp
					file rename $target.temp $target
					cg_zindex $target
				}
				foreach file $list {
					file delete $file
					if {[file extension $file] eq ".lz4"} {file delete $file.lz4i}
					if {[file extension $file] eq ".zst"} {file delete $file.zsti}
					file delete -force [gzroot $file].temp
					file delete [analysisinfo_file $file]
				}
			}
			incr pos
		}
		cleanup_job cleanup-sv_${method}_[file tail $bamfile] [list {*}$regfiles] $resultfiles
		cd $keeppwd
		return $resultfiles
	}
}

proc cg_sv {args} {
	set args [job_init {*}$args]
	sv_job {*}$args
	job_wait
}
