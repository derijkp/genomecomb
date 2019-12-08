proc bam_markduplicates_job {args} {
	upvar job_logdir job_logdir
	set method samtools
	set threads 1
	set skips {}
	cg_options bam_markduplicates args {
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
	job bamremdup-$tail {*}$skips -deps {
		$src
	} -targets {
		$dest $dest.analysisinfo
	} -vars {
		oformat method threads
	} -code {
		cg_bam_markduplicates -method $method -outputformat $oformat -threads $threads $dep $target
	}
}

proc cg_bam_markduplicates {args} {
	upvar job_logdir job_logdir
	set method samtools
	set inputformat -
	set outputformat -
	set sourcefile -
	set resultfile -
	set threads 1
	set skips {}
	cg_options bam_markduplicates args {
		-method {
			if {$value ni {1 picard biobambam samtools sam}} {error "bam_markduplicates: unsupported -method $value"}
			set method $value
		}
		-inputformat - -if {
			set inputformat $value
		}
		-outputformat - -of {
			set outputformat $value
		}
		-threads {
			set threads $value
		}
	} {sourcefile resultfile} 0 2
	if {$inputformat eq "-"} {set inputformat [ext2format $sourcefile bam {bam cram sam}]}
	if {$outputformat eq "-"} {set outputformat [ext2format $resultfile bam {bam cram sam}]}
	if {$method eq "1"} {set method samtools}
	set tail [file tail $sourcefile]
	putslog "removing duplicates $tail"
	if {$method eq "picard"} {
		if {$outputformat eq "cram"} {error "cram output not supported by bam_markduplicates using picard method"}
		analysisinfo_write $sourcefile $resultfile removeduplicates picard removeduplicates_version [version picard]
		puts "removing duplicates"
		set scratchdir [scratchdir]/picard
		file mkdir $scratchdir
		set opts {}
		set optsio {}
		lappend opts COMPRESSION_LEVEL=[defcompressionlevel 5]
		if {$sourcefile eq "-"} {
			set o [open $scratchdir/source.$format w]
			fconfigure $o -translation binary
			fcopy stdin $o
			set sourcefile $scratchdir/source.$format
		}
		if {$resultfile eq "-"} {
			set tempresult [tempfile]
		} else {
			set tempresult [filetemp $resultfile 0 1]
		}
		picard MarkDuplicates	I=$sourcefile	O=$tempresult \
			METRICS_FILE=$resultfile.dupmetrics TMP_DIR=[scratchdir]/picard \
			{*}$opts {*}$optsio 2>@ stderr
		if {[file exists $tempresult]} {
			file rename -force $tempresult $resultfile
		} else {
			file copy $tempresult >@ stdout
		}
	} elseif {$method eq "biobambam"} {
		if {$outputformat eq "cram"} {error "cram output not supported by bam_markduplicates using biobambam method"}
		analysisinfo_write $sourcefile $resultfile removeduplicates biobambam removeduplicates_version [version biobambam]
		set opts {}
		set optsio {}
		lappend opts level=[defcompressionlevel -1]
		if {$resultfile eq "-"} {
			lappend optsio >@ stdout
		} else {
			set tempresult [filetemp $resultfile 0 1]
			lappend opts O=$tempresult
		}
		if {$sourcefile eq "-"} {
			lappend optsio <@ stdin
		} else {
			lappend opts I=$sourcefile
		}
		biobambam bammarkduplicates2 M=$resultfile.dupmetrics rmdup=0 markthreads=1 tmpfile=[scratchfile] \
			{*}$opts {*}$optsio 2>@ stderr
		if {$resultfile ne "-"} {
			file rename -force $tempresult $resultfile
		}
	} else {
		analysisinfo_write $sourcefile $resultfile removeduplicates samtools removeduplicates_version [version samtools]
		set opts {} ; set optsio {}
		set compressionlevel [defcompressionlevel -1]
		if {$compressionlevel != -1} {lappend opts --output-fmt-option level=$compressionlevel}
		if {$resultfile ne "-"} {
			set tempresult [filetemp $resultfile 0 1]
		} else {
			set tempresult -
			lappend optsio >@ stdout
		}
		if {$sourcefile eq "-"} {
			lappend optsio <@ stdin
		}
		puts stderr [list samtools markdup --output-fmt $outputformat -l 500 \
			{*}$opts $sourcefile $tempresult {*}$optsio]
		exec samtools markdup --output-fmt $outputformat -l 500 \
			{*}$opts $sourcefile $tempresult {*}$optsio 2>@ stderr
		if {$resultfile ne "-"} {
			file rename -force $tempresult $resultfile
		}
	}
}

#proc cg_bam_markduplicates {args} {
#	set args [job_init {*}$args]
#	set result [bam_markduplicates_job {*}$args]
#	job_wait
#	return $result
#}
