proc bam_sort_job {args} {
	upvar job_logdir job_logdir
	set method samtools
	set sort coordinate
	set inputformat bam
	set threads 1
	set skips {}
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
		analysisinfo_write $dep $target bamsort $method bamsort_version [version $method]
		file delete $target.temp
		if {$method eq "alreadysorted"} {
			hardlink $dep $target
		} else {
			bam_sort -method $method -sort $sort -threads $threads -inputformat $inputformat $dep $target.temp
			file rename -force $target.temp $target
		}
	}
}

proc bam_sort {args} {
	set method samtools
	set sort coordinate
	set inputformat bam
	set threads 1
	cg_options bam_sort args {
		-method {
			if {$value ni {biobambam samtools}} {error "bamsort: unsupported -method $value"}
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
	} {sourcefile resultfile}
	if {$method eq "biobambam"} {
		set tempresult [filetemp $resultfile]
		catch_exec bamsort I=$sourcefile inputformat=$inputformat SO=$sort tmpfile=[scratchfile] index=1 indexfilename=$tempresult.bai O=$tempresult
		file rename -force $tempresult.bai $resultfile.bai
		file rename -force $tempresult $resultfile
	} else {	
		set opts {}
		if {$sort eq "name"} {
			lappend opts -n
		}
		if {[catch {version samtools 1}]} {
			# version < 1
			if {[catch {exec samtools sort {*}$opts $sourcefile $resultfile.temp 2>@ stdout} msg]} {
				error $msg
			}
			if {[file exists $resultfile.temp.bam]} {
				file rename -force $resultfile.temp.bam $resultfile
			} else {
				file rename -force $resultfile.temp $resultfile
			}
		} else {
			set oformat [string toupper [string range [file extension $resultfile] 1 end]]
			if {$oformat ni "BAM CRAM SAM"} {
				set oformat [string toupper [string range [file extension [file root $resultfile]] 1 end]]
			}
			if {$oformat ni "BAM CRAM SAM"} {
				set oformat BAM
			}
			if {[catch {exec samtools sort {*}$opts --threads $threads -T [scratchfile] \
				-O $oformat -o $resultfile.temp[file extension $resultfile] $sourcefile 2>@ stderr} msg]} {
				error $msg
			}
			file rename -force $resultfile.temp[file extension $resultfile] $resultfile
		}
	}
}
