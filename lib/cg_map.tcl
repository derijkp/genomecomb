proc map_readgroupdata {readgroupdata sample} {
	array set a [list PL illumina LB solexa-123 PU $sample SM $sample]
	if {$readgroupdata ne ""} {
		array set a $readgroupdata
	}
	set readgroupdata [array get a]
}

proc methods_map {args} {
	set cmd map
	set supportedmethods [list bwa]
	set start [expr {[string length cg_${cmd}] + 1}]
	foreach name [array names ::auto_index cg_${cmd}_*] {
		list_addnew supportedmethods [string range $name $start end]
	}
	if {[llength $args]} {
		set value [lindex $args 0]
		if {$value eq "1"} {set value [lindex $supportedmethods 0]}
		if {$value ni $supportedmethods} {error "$cmd: unsupported -method $value"}
		return $value
	} else {
		return $supportedmethods
	}
}

proc map_mem {method {mem {}}} {
	if {[auto_load map_mem_$method]} {
		return [map_mem_$method $mem $threads]
	}
	if {$mem eq ""} {set mem 5G}
	return $mem
}

proc map_job {args} {
	upvar job_logdir job_logdir
	set cmdline "[list cd [pwd]] \; [list cg map {*}$args]"
	set paired 1
	set preset {}
	set readgroupdata {}
	set skips {}
	set threads 2
	set keepsams 0
	set sort coordinates
	set mergesort 0
	set fixmate 1
	set method bwa
	set mem {}
	cg_options map args {
		-method {
			set method [methods_map $value]
		}
		-paired - -p {
			set paired $value
		}
		-x - -preset - -p {
			set preset $value
		}
		-readgroupdata {
			set readgroupdata $value
		}
		-fixmate {
			set fixmate $value
		}
		-sort {
			if {$value eq "1"} {set value coordinates}
			if {$value ni "coordinates c nosort name"} {error "-sort must be coordinates, nosort or name"}
			set sort $value
		}
		-mergesort {
			set mergesort $value
		}
		-skips {
			set skips $value
		}
		-keepsams {
			set keepsams $value
		}
		-threads - -t {
			set threads $value
		}
		-mem {
			set mem $value
		}
	} {result refseq sample fastqfile1} 4 ... {
		align reads in fastq files to a reference genome
	}
	set files [list $fastqfile1 {*}$args]
	set result [file_absolute $result]
	set refseq [refseq $refseq]
	set resultdir [file dir $result]
	if {![info exists job_logdir]} {
		job_logdir [file dir $result]/log_jobs
	}
	job_logfile $resultdir/map_${method}_[file tail $result] $resultdir $cmdline \
		{*}[versions minimap2]
	# start
	array set a [list PL illumina LB solexa-123 PU $sample SM $sample]
	if {$readgroupdata ne ""} {
		array set a $readgroupdata
	}
	set readgroupdata [array get a]

	dbdir [file dir $refseq]
	set resultbase [file root $result]
	set samfiles {}
	set asamfiles {}
	set num 1
	if {!$paired} {
		foreach file $files {
			set name [file root [file tail $file]]
			set target $resultbase-$name.sam.zst
			lappend samfiles $target
			set analysisinfo [gzroot $target].analysisinfo
			lappend asamfiles $analysisinfo
			job map_${method}-$sample-$name -mem [job_mempercore [map_mem $method $mem] $threads] -cores $threads \
			-skip [list $resultbase.bam] {*}$skips \
			-deps {
				$refseq $file
			} -targets {
				$target $analysisinfo
			} -vars {
				method sort mergesort preset sample readgroupdata fixmate paired threads refseq file
			} -code {
				set tempfile [filetemp $target 1 1]
				if {!$mergesort || $sort eq "nosort"} {
					cg map_${method} -paired $paired	-preset $preset \
						-readgroupdata $readgroupdata -fixmate $fixmate \
						-threads $threads \
						$tempfile $refseq $sample $file
				} else {
					exec cg map_${method} -paired $paired	-preset $preset \
						-readgroupdata $readgroupdata -fixmate $fixmate \
						-threads $threads \
						-.sam $refseq $sample $file \
						| cg _sam_sort_gnusort {*}[compresspipe -.*.sam.zst 1] > $tempfile
				}
				result_rename $tempfile $target
			}
		}
	} else {
		if {[expr {[llength $files]%2}]} {
			error "mapping needs even number of files for paired analysis"
		}
		foreach {file1 file2} $files {
			set name [file root [file tail $file1]]
			set target $resultbase-$name.sam
			lappend samfiles $target
			set analysisinfo [gzroot $target].analysisinfo
			lappend asamfiles $analysisinfo
			job map_${method}-$sample-$name -mem [job_mempercore 5G $threads] -cores $threads \
			-skip [list $resultbase.bam] {*}$skips \
			-deps {
				$refseq $file1 $file2
			} -targets {
				$target $analysisinfo
			} -vars {
				method preset sample readgroupdata fixmate paired threads refseq file1 file2
			} -code {
				set tempfile [filetemp $target 1 1]
				if {!$mergesort || $sort eq "nosort"} {
					cg map_${method} -paired $paired	-preset $preset	\
						-readgroupdata $readgroupdata -fixmate $fixmate \
						-threads $threads \
						$tempfile $refseq $sample $file1 $file2
				} else {
					cg map_${method} -paired $paired	-preset $preset	\
						-readgroupdata $readgroupdata -fixmate $fixmate \
						-threads $threads \
						-.sam $refseq $sample $file1 $file2 \
						| cg _sam_sort_gnusort {*}[compresspipe -.*.sam.zst 1] > $tempfile
				}
				result_rename $tempfile $target
			}
		}
	}
	if {!$mergesort || $sort eq "nosort"} {
		sam_catmerge_job {*}$skips -sort $sort -name map_${method}2bam-$sample -refseq $refseq \
			-deletesams [string is false $keepsams] -threads $threads -- $result {*}$samfiles
	} else {
		sam_merge_job {*}$skips -sort $sort -name map_${method}2bam-$sample -refseq $refseq \
			-deletesams [string is false $keepsams] -threads $threads -- $result {*}$samfiles
	}
}

proc cg_map {args} {
	set args [job_init {*}$args]
	unset job_logdir
	set return [map_job {*}$args]
	job_wait
	return $return
}
