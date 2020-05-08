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
			set distrreg [distrreg_checkvalue $value]
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
		default {
			lappend var_opts $key $value
		}
	} {bamfile refseq} 2 2
	set bamfile [file_absolute $bamfile]
	set refseq [file_absolute $refseq]
	set destdir [file dir $bamfile]
	# logfile
	set tools {gatk picard java gnusort8 zst os}
	catch {lappend tools {*}[var_${method}_tools]}
	set tools [list_remdup $tools]
	job_logfile $destdir/var_${method}_[file tail $bamfile] $destdir $cmdline \
		{*}[versions {*}$tools]
	set cmdopts {}
	# check if regionfile is supported
	catch {var_${method}_job} temp
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
	if {!$supportsregionfile && !$supportsregion} {set distrreg 0}
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
			-split $split -threads $threads -cleanup $cleanup $bamfile $refseq
	} else {
		# check what the resultfiles are for the method
		set resultfiles [var_${method}_job -resultfiles 1 {*}$var_opts -opts $opts \
			-pre $pre \
			-datatype $datatype \
			-split $split $bamfile $refseq]
		set skips [list -skip [list_remove  $resultfiles {}]]
		# if {[jobtargetexists $resultfiles [list $refseq $bamfile $regionfile]]} return
		foreach {varfile sregfile varallfile vcffile} $resultfiles break
		set file [file tail $bamfile]
		set root [file_rootname $varfile]
		# start
		## Create sequencing region files
		set keeppwd [pwd]
		cd $destdir
		set indexdir [gzroot $varfile].index
		file mkdir $indexdir
		set regions [list_remove [distrreg_regs $distrreg $refseq] unaligned]
		set basename [gzroot [file tail $varallfile]]
		if {$supportsregionfile} {
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
				# puts "cg select -f \'$header\' $regionfile | $appdir/bin/distrreg $indexdir/$basename- \'$regions\' 0 1 2 1 \#"
				cg select -f $header $regionfile | $appdir/bin/distrreg $indexdir/$basename- .bed 0 $regions 0 1 2 1 \#
			}
		}
		set todo {}
		# Produce variant calls
		set ibam $indexdir/[file tail $bamfile]
		set indexext [indexext $bamfile]
		mklink $bamfile $ibam
		mklink $bamfile.$indexext $ibam.$indexext
		defcompressionlevel 1
		if {$supportsregionfile} {
			foreach region $regions regfile $regfiles {
				lappend todo [var_${method}_job {*}$var_opts -opts $opts {*}$skips \
					-datatype $datatype \
					-rootname $root-$region -regionfile $regfile \
					-split $split -threads $threads -cleanup $cleanup $ibam $refseq]
			}
		} else {
			foreach region $regions {
				lappend todo [var_${method}_job {*}$var_opts -opts $opts {*}$skips \
					-datatype $datatype \
					-rootname $root-$region -region $region \
					-split $split -threads $threads -cleanup $cleanup $ibam $refseq]
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
			job var_combineresults-$resultfile {*}$skips -deps $deps -rmtargets $list -targets {
				$resultfile
			} -vars {
				list regionfile method bamfile
			} -code {
				set analysisinfo [gzroot $dep].analysisinfo
				if {[file exists $analysisinfo]} {
					set root $method-[file_rootname [file tail $bamfile]]
					analysisinfo_copy $analysisinfo [gzroot $target].analysisinfo \
						[list varcaller_region [file tail $regionfile] \
						sample $root]
					exec touch [gzroot $target].analysisinfo
				}
				# using bsort because
				# files are xxx-100 -> would sort reverse of what we want because of the -
				if {[file extension [gzroot $target]] in ".vcf .gvcf"} {
					cg vcfcat -i 1 -o $target {*}[bsort [jobglob {*}$list]]
				} elseif {[lindex [split [file tail $target] -_] 0] eq "sreg"} {
					cg cat -c f {*}[bsort [jobglob {*}$list]] | cg regjoin {*}[compresspipe $target] > $target.temp
					file rename $target.temp $target
					cg_zindex $target
				} else {
					cg cat -c f {*}[bsort [jobglob {*}$list]] {*}[compresspipe $target] > $target.temp
					file rename $target.temp $target
					cg_zindex $target
				}
			}
			lappend cleanupfiles {*}$list
		}
		lappend cleanupfiles $indexdir
		cleanup_job -forcedirs 1 -delassociated 1 cleanup-var_${method}_[file tail $varfile] $cleanupfiles $resultfiles
		cd $keeppwd
		return $resultfiles
	}
}

proc cg_var {args} {
	set args [job_init {*}$args]
	var_job {*}$args
	job_wait
}
