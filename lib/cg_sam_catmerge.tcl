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
		-outputformat {
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
			# exec samcat {*}$deps | bamsort SO=coordinate tmpfile=[scratchfile] index=1 indexfilename=$target.bai inputformat=sam > $target.temp 2>@ stderr
			if {$sort eq "nosort"} {
				if {$outputformat eq "bam"} {
					exec samcat {*}$deps | samtools view --threads $threads -b -o $target.temp - 2>@ stderr
				} elseif {$outputformat eq "cram"} {
					exec samcat {*}$deps | samtools view --threads $threads -c -o $target.temp - 2>@ stderr
				} else {
					exec samcat {*}$deps > $target.temp 2>@ stderr
				}
			} else {
				exec samcat {*}$deps | samtools sort {*}$sortopt --threads $threads -T [scratchfile] -O $outputformat -o $target.temp 2>@ stderr
			}
		} msg]} {
			error $msg
		}
		file rename -force $target.temp $target
		foreach dep $rmfiles {
			file delete $dep
		}
	}
	if {$index} {	
		job $name-index {*}$skips -optional $optional -force $force -deps {$resultfile} -targets {$resultfile.bai} -code {
			exec samtools index $dep
		}
		
	}
}

proc cg_sam_catmerge {args} {
	set args [job_init {*}$args]
	unset job_logdir
	set result [sam_catmerge_job {*}$args]
	job_wait
	return $result
}
