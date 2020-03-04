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

proc map_mem {method {mem {}} {threads 1}} {
	if {[auto_load map_mem_$method]} {
		return [map_mem_$method $mem $threads]
	}
	if {$mem eq ""} {set mem 5G}
	job_mempercore $mem $threads
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
	set sort coordinate
	set mergesort 0
	set maxopenfiles {}
	set fixmate 1
	set method bwa
	set mem {}
	set compressionlevel {}
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
			if {$value eq "1"} {set value coordinate}
			if {$value ni "coordinate c nosort name"} {error "-sort must be coordinate, nosort or name"}
			set sort $value
		}
		-mergesort {
			set mergesort $value
		}
		-maxopenfiles {
			set maxopenfiles $value
		}
		-skips {
			set skips $value
		}
		-skip {
			lappend skips -skip $value
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
		-compressionlevel {
			set compressionlevel $value
		}
	} {result refseq sample fastqfile1} 4 ... {
		align reads in fastq files to a reference genome
	}
	set files [list $fastqfile1 {*}$args]
	set result [file_absolute $result]
	set refseq [refseq $refseq]
	set resultdir [file dir $result]
	file mkdir $resultdir
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
	#
	dbdir [file dir $refseq]
	set resultbase [file root $result]
	set samfiles {}
	set asamfiles {}
	set num 1
	if {($paired && [llength $files] == 2) || (!$paired && [llength $files] == 1)} {
		set analysisinfo [analysisinfo_file $result]
		set file [lindex $files 0]
		job map_${method}-[file tail $result] {*}$skips \
			-mem [map_mem $method $mem $threads] -cores $threads \
		-deps [list {*}$files $refseq] -targets {
			$result $analysisinfo
		} -vars {
			result method sort preset sample readgroupdata fixmate paired threads refseq files compressionlevel
		} -code {
			set tempfile [filetemp_ext $result]
			if {$sort eq "nosort"} {
				catch_exec cg map_${method} -paired $paired	-preset $preset \
					-readgroupdata $readgroupdata -fixmate $fixmate \
					-threads $threads \
					$tempfile $refseq $sample {*}$files
			} else {
				analysisinfo_write [lindex $files 0] $result aligner $method aligner_version [version $method] reference [file2refname $refseq] aligner_paired $paired aligner_sort gnusort aligner_sort_version [version gnusort8]
				if {[file_ext $result] eq ".cram"} {set addm5 1} else {set addm5 0}
				catch_exec cg map_${method} -paired $paired	-preset $preset \
					-readgroupdata $readgroupdata -fixmate $fixmate \
					-threads $threads \
					-.sam $refseq $sample {*}$files \
					| cg _sam_sort_gnusort $sort $threads $refseq $addm5 \
					{*}[convert_pipe -.sam $tempfile -compressionlevel $compressionlevel -refseq $refseq -threads $threads -endpipe 1]
			}
			result_rename $tempfile $result
		}
		return $result
	}
	set indexdir [gzroot $result].index
	file mkdir $indexdir
	if {!$paired} {
		foreach file $files {
			set name [file root [file tail $file]]
			set target $resultbase-$name.sam.zst
			lappend samfiles $target
			set analysisinfo [analysisinfo_file $target]
			lappend asamfiles $analysisinfo
			job map_${method}-$sample-$name -mem [map_mem $method $mem $threads] -cores $threads \
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
					if {[file_ext $target] eq ".cram"} {set addm5 1} else {set addm5 0}
					exec cg map_${method} -paired $paired	-preset $preset \
						-readgroupdata $readgroupdata -fixmate $fixmate \
						-threads $threads \
						-.sam $refseq $sample $file \
						| cg _sam_sort_gnusort $sort $threads $refseq $addm5 {*}[compresspipe -.*.sam.zst 1] > $tempfile
				}
				result_rename $tempfile $target
			}
		}
	} else {
		if {[expr {[llength $files]%2}]} {
			error "mapping needs even number of files for paired analysis"
		}
		foreach {file1 file2} $files {
			set name [file tail [file_root $file1]]
			set target $indexdir/[file_root [file tail $result]]-$name.sam.zst
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
				method mergesort preset sample readgroupdata fixmate paired threads refseq file1 file2
			} -code {
				set tempfile [filetemp_ext $target]
				if {!$mergesort || $sort eq "nosort"} {
					cg map_${method} -paired $paired	-preset $preset	\
						-readgroupdata $readgroupdata -fixmate $fixmate \
						-threads $threads \
						$tempfile $refseq $sample $file1 $file2
				} else {
					if {[file_ext $target] eq ".cram"} {set addm5 1} else {set addm5 0}
					cg map_${method} -paired $paired	-preset $preset	\
						-readgroupdata $readgroupdata -fixmate $fixmate \
						-threads $threads \
						-.sam $refseq $sample $file1 $file2 \
						| cg _sam_sort_gnusort $sort $threads $refseq $addm5 {*}[compresspipe -.*.sam.zst 1] > $tempfile
				}
				result_rename $tempfile $target
			}
		}
	}
	sam_catmerge_job {*}$skips -sort $sort -mergesort $mergesort -maxopenfiles $maxopenfiles -name map_${method}2bam-$sample -refseq $refseq \
		-deletesams [string is false $keepsams] -threads $threads -- $result {*}$samfiles
	return $result
}

proc cg_map {args} {
	set args [job_init {*}$args]
	unset job_logdir
	set return [map_job {*}$args]
	job_wait
	return $return
}
