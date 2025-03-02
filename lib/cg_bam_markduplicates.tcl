proc methods_bam_markduplicates {args} {
	set cmd bam_markduplicates
	set supportedmethods {samtools sam picard biobambam}
	if {[llength $args]} {
		set value [lindex $args 0]
		if {$value eq "1"} {set value [lindex $supportedmethods 0]}
		if {$value ni $supportedmethods} {error "$cmd: unsupported -method $value"}
		return $value
	} else {
		return $supportedmethods
	}
}

proc bam_markduplicates_job {args} {
	upvar job_logdir job_logdir
	set method samtools
	set threads 1
	set skips {}
	cg_options bam_markduplicates args {
		-method {
			set method [methods_bam_markduplicates $value]
		}
		-skip {
			lappend skips -skip $value
		}
		-threads {
			set threads $value
		}
	} {src dest}
	set destanalysisinfo [analysisinfo_file $dest]
	if {$method eq "1"} {set method samtools}
	set tail [file tail $src]
	putslog "removing duplicates $tail"
	set oformat [string range [file extension $dest] 1 end]
	job bamremdup-$tail {*}$skips -deps {
		$src
	} -targets {
		$dest $destanalysisinfo
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
	set refseq {}
	cg_options bam_markduplicates args {
		-method {
			set method [methods_bam_markduplicates $value]
		}
		-inputformat - -if {
			set inputformat $value
		}
		-outputformat - -of {
			set outputformat $value
		}
		-refseq {
			set refseq $value
		}
		-threads {
			set threads $value
		}
	} {sourcefile resultfile} 0 2
	if {$inputformat eq "-"} {set inputformat [ext2format $sourcefile bam {bam cram sam}]}
	if {$outputformat eq "-"} {set outputformat [ext2format $resultfile bam {bam cram sam}]}
	set tail [file tail $sourcefile]
	if {$method eq "picard"} {
		if {$outputformat eq "cram"} {error "cram output not supported by bam_markduplicates using picard method"}
		analysisinfo_write $sourcefile $resultfile removeduplicates picard removeduplicates_version [version picard]
		set scratchdir [scratchdir]/picard
		file mkdir $scratchdir
		set opts {}
		set optsio {}
		lappend opts COMPRESSION_LEVEL=[defcompressionlevel 5]
		if {$sourcefile eq "-"} {
			set newsourcefile $scratchdir/source.$inputformat
			set o [open $newsourcefile w]
			copybinary stdin $o
			close $o
			set sourcefile $newsourcefile
		}
		if {$resultfile eq "-"} {
			set tempresult [tempfile].$outputformat
			set metricsfile [tempfile]
		} else {
			set tempresult [filetemp $resultfile 0 1]
			set metricsfile $resultfile.dupmetrics
		}
		picard MarkDuplicates	I=$sourcefile	O=$tempresult \
			METRICS_FILE=$metricsfile TMP_DIR=[scratchdir]/picard \
			{*}$opts {*}$optsio
		if {$resultfile ne "-"} {
			file rename -force -- $tempresult $resultfile
		} else {
			file2stdout $tempresult
		}
	} elseif {$method eq "biobambam"} {
		if {$outputformat eq "cram"} {error "cram output not supported by bam_markduplicates using biobambam method"}
		analysisinfo_write $sourcefile $resultfile removeduplicates biobambam removeduplicates_version [version biobambam]
		set opts {}
		set optsio {}
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
		lappend opts level=[defcompressionlevel -1]
		biobambam bammarkduplicates2 M=$resultfile.dupmetrics rmdup=0 markthreads=1 tmpfile=[scratchfile] \
			{*}$opts {*}$optsio 2>@ stderr
		if {$resultfile ne "-"} {
			file rename -force -- $tempresult $resultfile
		}
	} else {
		analysisinfo_write $sourcefile $resultfile removeduplicates samtools removeduplicates_version [version samtools]
		set opts {} ; set optsio {}
		set compressionlevel [defcompressionlevel -1]
		set outputfmtoption {}
		if {$outputformat eq "cram"} {
			lappend outputfmtoption reference=[refseq $refseq]
		}
		if {$outputfmtoption ne ""} {
			lappend opts --output-fmt-option [join $outputfmtoption ,]
		}
		if {$resultfile ne "-"} {
			set tempresult [filetemp $resultfile 0 1]
			if {[file extension [gzroot -.$outputformat]] eq ".sam"} {
				lappend optsio {*}[convert_pipe -.sam $tempresult -endpipe 1]
				set outputformat sam
				set output -
			} else {
				set output $tempresult
			}
		} else {
			set output -
			if {[file extension [gzroot -.$outputformat]] eq ".sam"} {
				lappend optsio {*}[convert_pipe -.sam -.$outputformat -endpipe 1]
				set outputformat sam
			} else {
				lappend optsio >@ stdout
			}
		}
		if {$compressionlevel != -1 && $outputformat ne "sam"} {
			lappend outputfmtoption level=$compressionlevel
		}
		if {$sourcefile eq "-"} {
			lappend optsio <@ stdin
		}
		# puts stderr [list samtools markdup --output-fmt $outputformat -l 500 --threads $threads -T [scratchfile] \
			{*}$opts $sourcefile $output {*}$optsio]
		catch_exec samtools markdup --output-fmt $outputformat -l 500 --threads $threads -T [scratchfile] \
			{*}$opts $sourcefile $output {*}$optsio 2>@ stderr
		if {$resultfile ne "-"} {
			file rename -force -- $tempresult $resultfile
		}
	}
}
