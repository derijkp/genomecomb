proc sam_header_addm5 {header {refseq {}}} {
	if {![regexp {M5:} $header]} {
		dbdir $refseq
		if {![file exists [get ::env(REF_PATH) /complgen/refseq]/mapping.tsv]} {
			error "Could not find reference md5 mapping file (for cram): specify dbdir, reference or use REF_PATH"
		}
		set temp [file_read $::env(REF_PATH)/mapping.tsv]
		foreach line [split $temp \n] {
			set line [split $line \t]
			if {[llength $line] < 2} continue
			set a([lindex $line 1]) [lindex $line 2]
		}
		set newheader {}
		foreach line [split $header \n] {
			if {[regexp ^@SQ $line]} {
				if {![regexp {SN:([^ \t]+)} $line temp name]} {
					error "@SQ field without SN: for line $line"
				}
				set name [chr_clip $name]
				if {![info exists a($name)]} {
					set name chr$name
					if {![info exists a($name)]} {
						error "no md5 mapping found for sequence $name in $::env(REF_PATH)/mapping.tsv"
					}
				}
				append line "\tM5:$a($name)"
				lappend newheader $line
			} else {
				lappend newheader $line
			}
		}
		set header [join $newheader \n]
	}
	return $header
}

proc sam_catmerge_job {args} {
	set threads 1
	set force 0
	set deletesams 0
	set optional 0
	set index 0
	set skips {}
	set sort coordinates
	set distrreg 0
	set refseq ""
	cg_options sam_catmerge args {
		-name {
			set name $value
		}
		-sort {
			if {$value ni "coordinates c nosort merge"} {error "-sort must be coordinates, nosort or merge"}
			set sort $value
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
			set distrreg $value
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
		set name sam_catmerge-[file tail $resultfile]
	}
	set resultfile [file_absolute $resultfile]
	set outputformat [ext2format $resultfile bam {bam cram sam}]
	set outputformat [gzroot $outputformat]
	if {$outputformat eq "sam"} {set index 0}
	# job_logdir
	upvar job_logdir job_logdir
	if {![info exists job_logdir]} {
		job_logdir $resultfile.log_jobs
	}
	# run
	if {$deletesams} {
		set rmfiles $samfiles
		foreach file $samfiles {lappend rmfiles [gzroot $file].analysisinfo}
	} else {
		set rmfiles {}
	}
	set analysisinfofile [gzroot $resultfile].analysisinfo
	if {$distrreg in {0 {}}} {
		set regions {}
		set targets [list $resultfile $analysisinfofile]
		set regresults [list $resultfile]
	} else {
		set regions [distrreg_regs $distrreg $refseq]
		set basename [file_root $resultfile]
		set regresults {}
		foreach region $regions {
			lappend regresults [file_root $resultfile]-$region[file_ext $resultfile]
		}
		set targets $regresults
		lappend targets [gzroot [lindex $regresults 0]].analysisinfo
	}
	job $name -optional $optional -force $force -cores $threads {*}$skips \
	-deps $samfiles \
	-rmtargets $rmfiles \
	-targets $targets -vars {
		threads sort rmfiles outputformat refseq regresults distrreg resultfile regions
	} -code {
		puts "making $target ..."
		analysisinfo_write $dep $target cat_merge [version genomecomb]
		if {$sort eq "nosort"} {
			set pipe {}
			set opencmd {}
			if {$outputformat eq "bam"} {
				set incmd [list samcat {*}$deps]
				set outcmd [list samtools view --threads $threads -b - >]
			} elseif {$outputformat eq "cram"} {
				set refseq [refseq $refseq]
				set header [exec samtools view -H [lindex $deps 0]]
				set header [sam_header_addm5 $header $refseq]
				set incmd [list samcat -header $header {*}$deps]
				set outcmd [list samtools view -h --threads $threads -C -T $refseq - < ]
			} elseif {[gziscompressed $target]} {
				set incmd [list samcat {*}$deps]
				set outcmd [list {*}[lrange [compresspipe $target] 1 end] >]
			} else {
				set incmd [list samcat {*}$deps]
				set outcmd {}
			}
			if {![llength $regions]} {
				if {[llength $outcmd]} {set outcmd [list | {*}$outcmd]} else {set outcmd >}
				exec {*}$incmd {*}$outcmd $target.temp 2>@ stderr
			} else {
				exec {*}$incmd | distrreg [file_root $resultfile] [file_ext $resultfile].temp 1 $regions 2 3 3 0 @ $outcmd 2>@ stderr
			}
		} else {
			set header [exec samtools view -H [lindex $deps 0]]
			set headerlines [llength [split $header \n]]
			if {[regexp @HD $header]} {
				regsub {@HD[^\n]+} $header "@HD	VN:1.6	SO:coordinate" header
			} else {
				set header "@HD	VN:1.6	SO:coordinate\n$header"
			}
			if {$outputformat eq "cram"} {
				set refseq [refseq $refseq]
				set header [sam_header_addm5 $header $refseq]
				set outcmd [list samtools view -h --threads $threads -C -T $refseq - >]
			} elseif {$outputformat eq "bam"} {
				set outcmd [list samtools view -h --threads $threads -b - >]
			} elseif {[gziscompressed $target]} {
				set outcmd [list {*}[lrange [compresspipe $target] 1 end] >]
			} else {
				set outcmd {}
			}
			set sortopts {}
			if {$sort eq "merge" && [llength $deps] < [expr {[maxopenfiles]-2}]} {
				if {![llength $regions]} {
					if {[llength $outcmd]} {set outcmd [list | {*}$outcmd]} else {set outcmd >}
					exec mergesorted @ 0 $header {2 3} {*}$deps \
						{*}$outcmd $target.temp
				} else {
					exec mergesorted @ 0 $header {2 3} {*}$deps \
						| distrreg [file_root $resultfile]- [file_ext $resultfile].temp 1 $regions 2 3 3 0 @ $outcmd 2>@ stderr
						exec cat temp | distrreg [file_root $resultfile]- [file_ext $resultfile].temp 1 $regions 2 3 3 0 @ $outcmd 2>@ stderr
				}
			} else {
				if {![llength $regions]} {
					if {[llength $outcmd]} {set outcmd [list | {*}$outcmd]} else {set outcmd >}
					exec samcat -header $header {*}$deps \
						| gnusort8 --header-lines $headerlines --parallel $threads -T [scratchdir] -t \t -s -k3,3B -k4,4B \
						{*}$outcmd $target.temp
				} else {
					exec samcat -header $header {*}$deps \
						| gnusort8 --header-lines $headerlines --parallel $threads -T [scratchdir] -t \t -s -k3,3B -k4,4B \
						| distrreg [file_root $resultfile]- [file_ext $resultfile].temp 1 $regions 2 3 3 0 @ $outcmd 2>@ stderr
				}
			}
		}
		if {$regresults eq ""} {
			file rename -force -- $target.temp $target
		} else {
			foreach regresult $regresults {
				file rename -force -- $regresult.temp $regresult
			}
		}
		foreach dep $rmfiles {
			file delete $dep
		}
	}
	if {$index} {
		set target $resultfile.[indexext $resultfile]
		job $name-index {*}$skips -optional $optional -force $force -deps {
			$resultfile
		} -targets {
			$target
		} -code {
			exec samtools index $dep
		}
		
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
