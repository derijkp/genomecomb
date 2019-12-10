proc sam_header_addm5 header {
	if {![regexp {M5:} $header]} {
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
	set index 1
	set skips {}
	set sort coordinates
	cg_options sam_catmerge args {
		-name {
			set name $value
		}
		-sort {
			if {$value ni "coordinates names nosort c n"} {error "-sort must be coordinates, names or nosort"}
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
		-outputformat - -aliformat {
			switch $value {
				BAM - bam {set outputformat bam}
				SAM - sam {set outputformat sam}
				CRAM - cram {set outputformat cram}
				default {error "unknown outputformat $value"}
			}
		}
		-skips {
			set skips $value
		}
	} {resultfile samfile} 1 ... {
		merge sam files by concatenating (no problem with max open files) and then sorting them
	}
	set samfiles [list $samfile {*}$args]
	if {![info exists name]} {
		set name sam_catmerge-[file tail $resultfile]
	}
	set resultfile [file_absolute $resultfile]
	if {![info exists outputformat]} {
		if {[file extension $resultfile] eq ".sam"} {
			set outputformat sam
			set index 0
		} else {
			set outputformat bam
		}
	}
	# job_logdir
	upvar job_logdir job_logdir
	if {![info exists job_logdir]} {
		job_logdir $resultfile.log_jobs
	}
	# run
	if {[string index $sort 0] eq "n"} {set sortopt "-n"} else {set sortopt ""}
	if {$deletesams} {
		set rmfiles $samfiles
		foreach file $samfiles {lappend rmfiles [gzroot $file].analysisinfo}
	} else {
		set rmfiles {}
	}
	job $name -optional $optional -force $force -cores $threads \
	-deps $samfiles \
	-rmtargets $rmfiles \
	-targets {
		$resultfile $resultfile.analysisinfo
	} {*}$skips -vars {
		threads sort sortopt rmfiles outputformat
	} -code {
		puts "making $target"
		analysisinfo_write $dep $target cat_merge [version genomecomb]
		if {[catch {
			if {$sort eq "nosort"} {
				if {$outputformat eq "bam"} {
					exec samcat {*}$deps | samtools view --threads $threads -b -o $target.temp - 2>@ stderr
				} elseif {$outputformat eq "cram"} {
					set header [exec samtools view -H [lindex $deps 0]]
					set header [sam_header_addm5 $header]
					exec samcat -header $header {*}$deps | samtools view -h --threads $threads -C -T $::env(REFSEQ) -o $target.temp - 2>@ stderr
				} else {
					exec samcat {*}$deps > $target.temp 2>@ stderr
				}
			} elseif {$outputformat eq "cram"} {
				set header [exec samtools view -H [lindex $deps 0]]
				set header [sam_header_addm5 $header]
				set o [open "| samtools view -h --threads $threads -C -T $::env(REFSEQ) > $target.temp -" w]
				puts $o [string trim $header]
				exec samcat -header {} {*}$deps | gnusort8 -T [scratchdir] -t \t -s -k3,3B -k4,4B --parallel $threads >@ $o
				close $o
			} else {
				exec samcat {*}$deps | samtools sort {*}$sortopt --threads $threads -T [scratchfile] -O $outputformat -o $target.temp 2>@ stderr
			}
		} msg]} {
			error $msg
		}
		file rename -force -- $target.temp $target
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
	return $resultfile
}

proc cg_sam_catmerge {args} {
	set args [job_init {*}$args]
	unset job_logdir
	set result [sam_catmerge_job {*}$args]
	job_wait
	return $result
}
