proc bam_sort_job {args} {
	upvar job_logdir job_logdir
	set method samtools
	set sort coordinate
	set inputformat bam
	set threads 1
	set skips {}
	set refseq {}
	cg_options bam_sort args {
		-method {
			if {$value ni {biobambam samtools alreadysorted}} {error "bamsort: unsupported -method $value"}
			set method $value
		}
		-sort {
			if {$value ni {coordinate name hash}} {error "bamsort: unsupported -sort $value"}
			set sort $value
		}
		-inputformat {
			set inputformat $value
		}
		-threads {
			set threads $value
		}
		-refseq {
			set refseq $value
		}
		-skip {
			lappend skips -skip $value
		}
	} {sourcefile resultfile}
	set analysisinfo [gzroot $resultfile].analysisinfo
	job bamsort-[file tail $resultfile] -deps {
		$sourcefile
	} -targets {
		$resultfile $analysisinfo
	} -vars {
		method sort inputformat threads
	} {*}$skips -cores $threads -code {
		file delete $target.temp
		if {$method eq "alreadysorted"} {
			hardlink $dep $target
		} else {
			bam_sort -method $method -sort $sort -threads $threads -refseq $refseq -inputformat $inputformat $dep $target.temp
			file rename -force -- $target.temp $target
		}
	}
}

proc cg_bam_sort {args} {
	set method samtools
	set sort coordinate
	set inputformat -
	set outputformat -
	set threads 1
	set sourcefile -
	set resultfile -
	set refseq {}
	cg_options bam_sort args {
		-method {
			if {$value eq "1"} {set value samtools}
			if {$value ni {biobambam samtools}} {error "bamsort: unsupported -method $value"}
			set method $value
		}
		-sort {
			if {$value ni {coordinate name hash}} {error "bamsort: unsupported -sort $value"}
			set sort $value
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
		-refseq {
			set refseq $value
		}
	} {sourcefile resultfile} 0 2 {
		sort a bamfile
	}
	if {$inputformat eq "-"} {set inputformat [ext2format $sourcefile bam {bam cram sam}]}
	if {$outputformat eq "-"} {set outputformat [ext2format $resultfile bam {bam cram sam}]}
	analysisinfo_write $sourcefile $resultfile bamsort $method bamsort_version [version $method]
	if {$method eq "biobambam"} {
		set opts {}
		set optsio {}
		if {$sort eq "name"} {
			set sort queryname
		}
		if {$sourcefile ne "-"} {
			lappend opts I=$sourcefile
		} else {
			lappend optsio <@ stdin
		}
		if {$outputformat eq "sam"} {
			lappend optsio level=0 \| samtools view -h
			if {$resultfile ne "-"} {
				set tempresult [filetemp $resultfile]
				lappend optsio > $tempresult
			} else {
				lappend optsio >@ stdout
			}
		} else {
			lappend opts level=[defcompressionlevel -1]
			if {$resultfile ne "-"} {
				set tempresult [filetemp $resultfile]
				lappend opts O=$tempresult
			} else {
				lappend optsio >@ stdout
			}
		}
		if {$sort eq "coordinate" && $resultfile ne "-"} {
			lappend opts index=1 indexfilename=$tempresult.bai
		}
		catch_exec bamsort inputformat=$inputformat SO=$sort tmpfile=[scratchfile] {*}$opts {*}$optsio
		if {[info exists tempresult]} {
			catch {file rename -force -- $tempresult.bai $resultfile.bai}
			file rename -force -- $tempresult $resultfile
		}
	} else {	
		set opts {} ; set optsio {}
		if {$sort eq "name"} {
			lappend opts -n
		}
		set compressionlevel [defcompressionlevel -1]
		if {$compressionlevel != -1} {lappend opts -l $compressionlevel}
		if {$resultfile ne "-"} {
			set tempresult [filetemp $resultfile 0 1]
			lappend opts -o $tempresult
		} else {
			lappend optsio >@ stdout
		}
		if {$sourcefile ne "-"} {
			lappend opts $sourcefile
		} else {
			lappend optsio <@ stdin
		}
		if {$outputformat eq "cram"} {
			lappend opts --output-fmt-option reference=[refseq $refseq]
		}
		if {[catch {exec samtools sort --threads $threads -T [scratchfile] \
			-O [string toupper $outputformat] \
			 {*}$opts {*}$optsio 2>@ stderr
		} msg]} {
			error $msg
		}
		if {[info exists tempresult]} {
			file rename -force -- $tempresult $resultfile
		}
	}
}