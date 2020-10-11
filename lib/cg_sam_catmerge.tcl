proc sam_catmerge_job {args} {
	# job_logdir
	upvar job_logdir job_logdir
	set cmdline "[list cd [pwd]] \; [list cg var {*}$args]"
	set threads 1
	set force 0
	set deletesams 0
	set optional 0
	set index 0
	set skips {}
	set sort coordinate
	set distrreg 0
	set refseq ""
	set mergesort 0
	set maxopenfiles {}
	cg_options sam_catmerge args {
		-name {
			set name $value
		}
		-sort {
			if {$value in "1 c"} {set value coordinate}
			if {$value ni "coordinate name nosort"} {error "-sort must be coordinate, name or nosort"}
			set sort $value
		}
		-mergesort {
			set mergesort $value
		}
		-maxopenfiles {
			set maxopenfiles $value
		}
		-index {
			set index $value
		}
		-deletesams {
			set deletesams $value
		}
		-threads {
			set threads $value
		}
		-force {
			set force $value
		}
		-optional {
			set optional $value
		}
		-refseq {
			set refseq [refseq $value]
		}
		-distrreg {
			set distrreg [distrreg_checkvalue $value]
		}
		-skips {
			set skips $value
		}
		-skip {
			lappend skips -skip $value
		}
	} {resultfile samfile} 1 ... {
		merge sam files by concatenating (no problem with max open files) and then sorting them
	}
	set samfiles [list $samfile {*}$args]
	if {![info exists name]} {
		set name [job_relfile2name sam_catmerge- $resultfile]
	}
	set resultfile [file_absolute $resultfile]
	set outputformat [ext2format $resultfile bam {bam cram sam}]
	set outputformat [gzroot $outputformat]
	if {$outputformat eq "sam"} {set index 0}
	job_logfile [file dir $resultfile]/sam_catmerge [file dir $resultfile] $cmdline \
		{*}[versions samtools]
	# run

	set workdir [workdir $resultfile]
	set tempresultfile $workdir/[file tail $resultfile]

	if {$deletesams} {
		set rmfiles {}
		foreach file $samfiles {
			lappend rmfiles $file [index_file $file] [analysisinfo_file $file]
		}
		set rmfiles [list_remove $rmfiles {}]
		job_cleanup_add $tempresultfile $workdir
	} else {
		set rmfiles {}
	}
	set analysisinfofile [analysisinfo_file $resultfile]
	set regions [distrreg_regs $distrreg $refseq]
	if {![llength $regions]} {
		set targets [list $resultfile $analysisinfofile]
		set regresults [list $resultfile]
	} else {
		set basename [file_root $resultfile]
		set regresults {}
		set tempregresults {}
		foreach region $regions {
			lappend regresults [file_root $resultfile]-$region[file_ext $resultfile]
			lappend tempregresults [file_root $tempresultfile]-$region[file_ext $resultfile]
		}
		set targets $regresults
		lappend targets [gzroot [lindex $regresults 0]].analysisinfo
	}
	if {$sort eq "nosort" && $mergesort} {error "cannot combine -sort nosort with -mergesort 1"}
	if {![llength $regions] && $sort eq "name"} {error "cannot combine distrreg with name sorting"}
	set deps $samfiles
	job $name -optional $optional -force $force -cores $threads {*}$skips \
	-deps $samfiles -rmtargets $rmfiles -targets $targets -vars {
		samfiles rmfiles regresults tempregresults resultfile regions mergesort sort refseq threads maxopenfiles outputformat workdir tempresultfile
	} -code {
		set testsam [lindex $samfiles 0]
		if {![llength $regions]} {
			analysisinfo_write $testsam $resultfile sammerge genomecomb sammerge_version [version genomecomb] sammerge_sort $sort sammerge_mergesort $mergesort
		} else {
			foreach regresult $regresults {
				analysisinfo_write $testsam $regresult sammerge genomecomb sammerge_version [version genomecomb] sammerge_sort $sort sammerge_mergesort $mergesort
			}
		}
		if {[gziscompressed $testsam]} {
			set header [exec cg zcat $testsam | samtools view --no-PG -H]
		} else {
			set header [exec samtools view --no-PG -H $testsam]
		}
		if {[file_ext $resultfile] eq ".cram"} {
			set refseq [refseq $refseq]
			set header [sam_header_addm5 $header $refseq]
		}
		set outcmd [convert_pipe -.sam $resultfile -refseq $refseq -threads $threads]
		if {$sort eq "nosort"} {
			if {$mergesort} {error "cannot combine -sort nosort with -mergesort 1"}
			set pipe {}
			set opencmd {}
			if {$outputformat eq "cram"} {
				set incmd [list samcat -header $header {*}$deps]
			} else {
				set incmd [list samcat {*}$deps]
			}
			if {![llength $regions]} {
				if {$outcmd ne ""} {set outcmd [list | {*}$outcmd]}
				lappend outcmd >
				exec {*}$incmd {*}$outcmd $tempresultfile 2>@ stderr
			} else {
				exec {*}$incmd | distrreg [file_root $tempresultfile] [file_ext $resultfile] 1 $regions 2 3 3 0 @ $outcmd 2>@ stderr
			}
		} else {
			if {[regexp @HD $header]} {
				regsub {@HD[^\n]+} $header "@HD	VN:1.6	SO:$sort" header
			} else {
				set header "@HD	VN:1.6	SO:$sort\n$header"
			}
			set headerlines [llength [split $header \n]]
			# fputsvars ~/tmp/temp header deps regions resultfile outcmd headerlines threads
			if {!$mergesort} {
				if {$sort eq "coordinate"} {
					set sortopt {-k3,3N -k4,4N -k1,1N -k2,2N}
				} else {
					set sortopt {-k1,1N}
				}
				if {![llength $regions]} {
					if {$outcmd ne ""} {set outcmd [list | {*}$outcmd]}
					lappend outcmd >
					exec samcat -header $header {*}$deps \
						| gnusort8 --header-lines $headerlines --parallel $threads -T [scratchdir] -t \t -s {*}$sortopt \
						{*}$outcmd $tempresultfile
				} else {
					if {[llength $outcmd]} {lappend outcmd >}
					exec samcat -header $header {*}$deps \
						| gnusort8 --header-lines $headerlines --parallel $threads -T [scratchdir] -t \t -s {*}$sortopt \
						| distrreg [file_root $tempresultfile]- [file_ext $resultfile] 1 $regions 2 3 3 0 @ $outcmd 2>@ stderr
				}
			} else {
				if {$sort eq "coordinate"} {
					set sortopt {2 3}
				} else {
					set sortopt {0}
				}
				if {![llength $regions]} {
					if {$outcmd ne ""} {set outcmd [list | {*}$outcmd]}
					lappend outcmd >
					set finaloutcmd [list {*}$outcmd $tempresultfile]
				} else {
					if {[llength $outcmd]} {lappend outcmd >}
					set finaloutcmd [list | distrreg [file_root $tempresultfile]- [file_ext $resultfile] 1 $regions 2 3 3 0 @ $outcmd]
				}
				set maxopenfiles [maxopenfiles $maxopenfiles]
				set len [llength $deps]
				if {$len <= $maxopenfiles} {
					exec cg mergesorted -headerline 0 -commentchar @ -sortpos $sortopt \
						{*}$deps {*}$finaloutcmd
				} else {
					set workdir [scratchdir]/merge
					file delete -force $workdir
					file mkdir $workdir
					set todo $deps
					set len [llength $todo]
					set delete 0
					set num 1
					while {$len > $maxopenfiles} {
						set pos 0
						set newtodo {}
						while {$pos < $len} {
							incr num
							set part [lrange $todo $pos [expr {$pos+$maxopenfiles-1}]]
							incr pos $maxopenfiles
							# puts [list ../bin/tsv_paste {*}$deps]
							if {[llength $part] > 1} {
								set parttarget $workdir/paste.temp$num.zst
								exec cg mergesorted -commentchar @ -headerline 0 -header $header -sortpos $sortopt \
									{*}$part | cg zst -compressionlevel 1 > $parttarget.temp.zst
								# exec samtools merge {*}$sortopt -t $threads $parttarget.temp {*}$part
								if {$delete} {file delete {*}$part}
								file rename -force -- $parttarget.temp.zst $parttarget
							} elseif {!$delete} {
								set part [lindex $part 0]
								set parttarget $workdir/paste.temp$num[gzext $part]
								mklink $part $parttarget
							} else {
								set parttarget $workdir/paste.temp$num[gzext $part]
								set part [lindex $part 0]
								file rename $part $parttarget
							}
							lappend newtodo $parttarget
						}
						set delete 1
						set todo $newtodo
						set len [llength $todo]
					}
					exec cg mergesorted -commentchar @ -headerline 0 -header $header -sortpos $sortopt \
						{*}$todo {*}$finaloutcmd
					if {$delete} {file delete {*}$todo}
				}
			}
		}
		if {![llength $regions]} {
			result_rename $tempresultfile $resultfile
		} else {
			foreach tempregresult $tempregresults regresult $regresults {
				result_rename $tempregresult $regresult
			}
		}
		foreach rmfile $rmfiles {file delete -force $rmfile}
	}
	if {$index} {
		bam_index_job {*}$skips -optional $optional $resultfile
	}
	return $regresults
}

proc cg_sam_catmerge {args} {
	set args [job_init {*}$args]
	unset job_logdir
	set result [sam_catmerge_job {*}$args]
	job_wait
	return $result
}
