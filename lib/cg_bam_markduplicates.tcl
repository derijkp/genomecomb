proc bam_markduplicates_job {args} {
	upvar job_logdir job_logdir
	set method samtools
	set threads 1
	set skips {}
	cg_options bam_sort args {
		-method {
			if {$value ni {1 picard biobambam samtools sam}} {error "bam_markduplicates: unsupported -method $value"}
			set method $value
		}
		-skip {
			lappend skips -skip $value
		}
		-threads {
			set threads $value
		}
	} {src dest}
	if {$method eq "1"} {set method samtools}
	set tail [file tail $src]
	putslog "removing duplicates $tail"
	set oformat [string range [file extension $dest] 1 end]
	if {$method eq "picard"} {
		if {$oformat eq "cram"} {error "cram output not supported by bam_markduplicates using picard method"}
		job bamremdup-$tail -mem [job_mempercore 10G 2] -cores 2 -deps {
			$src
		} -targets {
			$dest $dest.analysisinfo
		} -vars {} {*}$skips -code {
			analysisinfo_write $dep $target removeduplicates picard removeduplicates_version [version picard]
			puts "removing duplicates"
			file mkdir [scratchdir]/picard
			picard MarkDuplicates	I=$dep	O=$target.temp METRICS_FILE=$target.dupmetrics TMP_DIR=[scratchdir]/picard 2>@ stderr >@ stdout
			file rename -force $target.temp $target
		}
	} elseif {$method eq "biobambam"} {
		if {$oformat eq "cram"} {error "cram output not supported by bam_markduplicates using biobambam method"}
		job bamremdup-$tail -deps {
			$src
		} -targets {
			$dest $dest.analysisinfo
		} -vars {} {*}$skips -code {
			analysisinfo_write $dep $target removeduplicates biobambam removeduplicates_version [version biobambam]
			biobambam bammarkduplicates2 I=$dep	O=$target.temp M=$target.dupmetrics rmdup=0 markthreads=1 tmpfile=[scratchfile] 2>@ stderr >@ stdout
			file rename -force $target.temp $target
		}
	} else {
		job bamremdup-$tail -deps {
			$src
		} -targets {
			$dest $dest.analysisinfo
		} -vars {oformat} {*}$skips -code {
			analysisinfo_write $dep $target removeduplicates samtools removeduplicates_version [version samtools]
			exec samtools markdup --output-fmt $oformat -l 500 $dep	$target.temp 2>@ stderr >@ stdout
			file rename -force $target.temp $target
		}
	}
}

proc cg_bam_markduplicates {args} {
	set args [job_init {*}$args]
	set result [bam_markduplicates_job {*}$args]
	job_wait
	return $result
}
