proc sam_catmerge_job {args} {
	upvar job_logdir job_logdir
	set threads 1
	set force 0
	set deletesams 0
	set optional 0
	set index 1
	set skips {}
	set sort coordinates
	set outputformat BAM
	cg_options sam_catmerge args {
		-name {
			set name $value
		}
		-sort {
			if {$value ni "coordinates names c n"} {error "-sort must be coordinates or names"}
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
				BAM - bam {set outputformat BAM}
				SAM - sam {set outputformat SAM}
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
	if {[string index $sort 0] eq "n"} {set sortopt "-n"} else {set sortopt ""}
	if {$deletesams} {set rmsamfiles $samfiles} else {set rmsamfiles {}}
	job $name -optional $optional -force $force -cores $threads \
	-deps $samfiles \
	-rmtargets $rmsamfiles \
	-targets {
		$resultfile $resultfile.analysisinfo
	} {*}$skips -vars {
		threads sortopt rmsamfiles outputformat
	} -code {
		puts "making $target"
		analysisinfo_write $dep $target
		if {[catch {
			# exec samcat {*}$deps | bamsort SO=coordinate tmpfile=[scratchfile] index=1 indexfilename=$target.bai inputformat=sam > $target.temp 2>@ stderr
			exec samcat {*}$deps | samtools sort {*}$sortopt --threads $threads -T [scratchfile] -O $outputformat -o $target.temp 2>@ stderr
		} msg]} {
			error $msg
		}
		file rename -force $target.temp $target
		foreach dep $rmsamfiles {
			file delete $dep [gzroot $dep].analysisinfo
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
