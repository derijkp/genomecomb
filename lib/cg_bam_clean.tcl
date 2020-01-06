proc bam_clean_job {args} {
	upvar job_logdir job_logdir
	set skips {}
	set sort 0
	set removeduplicates 0
	set realign 0
	set clipamplicons {}
	set regionfile {}
	set threads 2
	set refseq {}
	set resultfile ""
	set inputformat ""
	set outputformat ""
	set distrreg 0
	set keep 0
	cg_options bam_clean args {
		-sort {set sort $value}
		-refseq {set refseq $value}
		-removeduplicates {set removeduplicates $value}
		-realign {set realign $value}
		-clipamplicons {set clipamplicons $value}
		-threads {set threads $value}
		-regionfile {set regionfile $value}
		-outputformat {set outputformat $value}
		-distrreg {
			set distrreg $value
		}
		-keep {set keep $value}
		-skip {
			lappend skips -skip $value
		}
	} {sourcefile resultfile} 1 2
	set inputformat [gzroot [ext2format $sourcefile bam {bam cram sam}]]
	set inputcompressed [gzext $sourcefile]
	if {$resultfile ne ""} {
		if {$outputformat ne ""} {error "cannot use -outputformat option when resultfile is specified"}
		set outputformat [ext2format $resultfile bam {bam cram sam}]
	} elseif {$outputformat eq ""} {
		set outputformat bam
	}
	set dir [file dir $sourcefile]
	set file [file tail $sourcefile]
	set root [file_rootname $file pre]
	# precalc name and steps
	set steps 0
	set temproot $root
	if {$sort} {
		set temproot s$temproot
		if {$sort != 2} {incr steps}
	}
	if {$removeduplicates ne "0"} {
		set temproot d$temproot
		incr steps
	}
	if {$realign ne "0"} {
		set temproot r$temproot
		incr steps
	}
	if {$clipamplicons ne ""} {
		set temproot c$temproot
		incr steps
	}
	if {$resultfile eq ""} {
		set resultfile $dir/$pre$temproot.$outputformat
	}
	if {$steps == 0} {
		if {$outputformat eq $inputformat} {
			if {$resultfile eq $sourcefile} {
				return $resultfile
			} else {
				analysisinfo_write $dep $target bamclean genomecomb bamclean_version [version genomecomb]
				job bamclean-$root -deps {$sourcefile} -targets {$resultfile} -code {
					hardcopy $dep $target
				}
			}
		} else {
			job bamclean-$root -deps {$sourcefile} -targets {$resultfile} -vars {
				inputformat outputformat refseq
			} -code {
				analysisinfo_write $dep $target bamclean genomecomb bamclean_version [version genomecomb]
				if {[file size $dep] > 0} {
					exec {*}[convert_pipe -.$inputformat -.$outputformat -refseq $refseq] < $dep > $target
				} else {
					file_write $target ""
				}
			}
		}
		return $resultfile
	}
	# make pipe
	set stack [get ::stacktraceonerror 0]
	set deps {}
	lappend deps $sourcefile
	set pipe {}
	set addanalysisinfo [list bamclean genomecomb bamclean_version [version genomecomb]]
	set optsio {}
	if {$inputcompressed eq ""} {
		lappend optsio < $sourcefile
	} else {
		lappend pipe {*}[gzcat $sourcefile] $sourcefile
	}
	# start jobs
	# sort using default
	set curstep 0
	if {$sort ne "0"} {
		if {$sort != 2} {
			incr curstep
			if {$curstep != $steps} {
				set curoutputformat bam
				set compressionlevel 1
			} else {
				set curoutputformat $outputformat
				set compressionlevel [defcompressionlevel 5]
			}
			if {[llength $pipe]} {lappend pipe |}
			set sortm [methods_bam_sort $sort]
			lappend addanalysisinfo bamsort $sortm bamsort_version [version $sortm]
			lappend pipe cg bam_sort -stack $stack -method $sort \
				-inputformat $inputformat -outputformat $curoutputformat \
				-compressionlevel $compressionlevel -threads $threads -refseq $refseq
		}
	}
	if {$removeduplicates ne "0"} {
		incr curstep
		if {$curstep != $steps} {
			set curoutputformat bam
			set compressionlevel 1
		} else {
			set curoutputformat $outputformat
			set compressionlevel [defcompressionlevel 5]
		}
		if {[llength $pipe]} {lappend pipe |}
		set removeduplicatesm [methods_bam_markduplicates $removeduplicates]
		lappend addanalysisinfo removeduplicates $removeduplicatesm removeduplicates_version [version $removeduplicatesm]
		lappend pipe cg bam_markduplicates -stack $stack -method $removeduplicates \
			-inputformat $inputformat -outputformat $curoutputformat \
			-compressionlevel $compressionlevel -threads $threads -refseq $refseq
	}
	if {$realign ne "0"} {
		if {$regionfile eq ""} {set regionfile 3}
		if {[isint $regionfile]} {
			# extract regions with coverage >= $regionfile (for cleaning)
			set regionfile [bam2reg_job -mincoverage $regionfile {*}$skips \
				-refseq $refseq \
				$sourcefile]
		}
		lappend deps $regionfile
		incr curstep
		if {$curstep != $steps} {
			set curoutputformat bam
			set compressionlevel 1
		} else {
			set curoutputformat $outputformat
			set compressionlevel [defcompressionlevel 5]
		}
		if {[llength $pipe]} {lappend pipe |}
		set realign [methods_realign $realign]
		if {$realign eq "gatk"} {set realignm gatk3} else {set realignm $realign}
		lappend addanalysisinfo realign $realign realign_version [version $realignm]
		lappend pipe cg realign -stack $stack -method $realign -regionfile $regionfile -refseq $refseq \
			-inputformat $inputformat -outputformat $curoutputformat \
			-compressionlevel $compressionlevel -threads $threads
	}
	if {$clipamplicons ne ""} {
		incr curstep
		if {$curstep != $steps} {
			set curoutputformat bam
			set compressionlevel 1
		} else {
			set curoutputformat $outputformat
			set compressionlevel [defcompressionlevel 5]
		}
		if {[llength $pipe]} {lappend pipe |}
		lappend addanalysisinfo clipamplicons genomecomb clipamplicons_version [version genomecomb] clipampliconsfile [file tail $clipamplicons]
		lappend pipe cg sam_clipamplicons -stack $stack \
			-inputformat $inputformat -outputformat $curoutputformat -refseq $refseq \
			-compressionlevel $compressionlevel \
			$clipamplicons
	}
	lappend pipe {*}$optsio > $resultfile.temp
	set resultfileindex $resultfile.[indexext $resultfile]
	if {!$keep} {
		set rmtargets [list $sourcefile $sourcefile.[indexext $sourcefile]]
	} else {
		set rmtargets {}
	}
	job bamclean-$root -deps $deps -targets {
		$resultfile
	} -rmtargets $rmtargets -vars {
		pipe sourcefile resultfile keep addanalysisinfo
	} -code {
		# analysisinfo_pipe_file $dep $target
eputsvars addanalysisinfo
		analysisinfo_write $dep $target {*}$addanalysisinfo
		if {![sam_empty $dep]} {
			catch_exec {*}$pipe
			file rename -force -- $resultfile.temp $resultfile
			# analysisinfo_pipe_stop
		} else {
			hardcopy $dep $resultfile
		}
		if {!$keep} {file delete $sourcefile $sourcefile.[indexext $sourcefile]}
	}
	job bamclean_index-$root -optional 1 {*}$skips -deps {$resultfile} -targets {$resultfileindex} -code {
		putslog "making $target"
		exec samtools index $dep >@ stdout 2>@ stderr
	}
	return $resultfile
}

proc cg_bam_clean args {
	set args [job_init {*}$args]
	bam_clean_job {*}$args
	job_wait
}
