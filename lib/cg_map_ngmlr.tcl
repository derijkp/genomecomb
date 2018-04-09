proc ngmlr_refseq_job {refseq preset} {
	upvar job_logdir job_logdir
	if {![jobtargetexists [list $refseq-ht-13-2.2.ngm $refseq-enc.2.ngm] [list $refseq]]} {
		set tail [file tail $refseq]
		job ngmlr_2refseq-[file tail $refseq] -deps {$refseq} -targets {
			$refseq-ht-13-2.2.ngm $refseq-enc.2.ngm
		} -vars {preset} -code {
			set temp [tempfile]
			file_write $temp ""
			if {[catch {
				exec -ignorestderr ngmlr -x $preset -r $dep -q $temp 2>@ stderr
			} e]} {
				error $e
			}
		}
	}
	return [list $refseq $refseq-ht-13-2.2.ngm $refseq-enc.2.ngm]
}

proc map_ngmlr_job {args} {
	upvar job_logdir job_logdir
	set paired 1
	set keepargs $args
	set preset {}
	set readgroupdata {}
	set skips {}
	set keepsams 0
	set threads 2
	cg_options map_ngmlr args {
		-x - -preset - -p {
			set preset $value
		}
		-readgroupdata {
			set readgroupdata $value
		}
		-keepsams {
			set keepsams [true $value]
		}
		-skips {
			set skips $value
		}
		-threads - -t {
			set threads $value
		}
		-m - -maxopenfiles {
			set ::maxopenfiles $value
		}
	} {result refseq sample fastqfile1} 4 ... {
		align reads in fastq files to a reference genome using ngmlr
	}
	if {$preset eq ""} {set preset ont}
	set files [list $fastqfile1 {*}$args]
	set result [file_absolute $result]
	set refseq [file_absolute $refseq]
	set resultdir [file dir $result]
	file mkdir $result.index
	if {![info exists job_logdir]} {
		job_logdir $result.index/log_jobs
	}
	job_logfile $resultdir/map_ngmlr_[file tail $result] $resultdir \
		[list cg map_ngmlr_ {*}$keepargs] \
		{*}[versions ngmlr]
	#
	array set a [list PL illumina LB solexa-123 PU $sample SM $sample]
	if {$readgroupdata ne ""} {
		array set a $readgroupdata
	}
	set readgroupdata [array get a]
	set refdeps [ngmlr_refseq_job $refseq $preset]
	set ngmlr_refseq [list_pop refdeps 0]
	set resultbase [file root $result]
	set samfiles {}
	set num 1
	foreach file $files {
		set name [file root [file tail $file]]
		set target $result.index/$name.sam
		lappend samfiles $target
		set analysisinfo [gzroot $target].analysisinfo
		job ngmlr-$sample-$name -mem 5G -cores $threads -deps [list $ngmlr_refseq $file {*}$refdeps] \
		-targets {
			$target $analysisinfo
		} -vars {
			threads preset readgroupdata sample
		} -skip [list $resultbase.bam] {*}$skips -code {
			puts "making $target"
			foreach {ngmlr_refseq fastq} $deps break
			analysisinfo_write $fastq $target aligner ngmlr aligner_version [version ngmlr] reference [file2refname $ngmlr_refseq] aligner_paired 0
			set rg {}
			foreach {key value} $readgroupdata {
				lappend rg "$key:$value"
			}
			exec ngmlr -x $preset -t $threads -r $ngmlr_refseq -q $fastq > $target.temp 2>@ stderr
			file rename -force $target.temp $target
		}
	}
	# if {$keepsams} {set deletesams 0} else {set deletesams 1}
	if {$keepsams} {set rmsamfiles {}} else {set rmsamfiles $samfiles}
	job ngmlr_sort2bam-$sample -cores $threads -deps $samfiles -rmtargets $rmsamfiles -targets {
		$result $result.analysisinfo
	} -vars {
		resultbase rmsamfiles threads
	} {*}$skips -code {
		puts "making $target"
		analysisinfo_write $dep $target
		if {[catch {
			exec samcat {*}$deps | samtools sort -T [scratchfile] -@ $threads -o $target.temp 2>@ stderr
		}]} {
			error $msg
		}
		file rename -force $target.temp $target
		foreach dep $rmsamfiles {
			file delete $dep [gzroot $dep].analysisinfo
		}
	}
	job ngmlr_bamindex-$sample -deps {$result} -targets {$result.bai} -code {
		exec samtools index $dep
	}
}

proc cg_map_ngmlr {args} {
	set args [job_init {*}$args]
	unset job_logdir
	set result [map_ngmlr_job {*}$args]
	job_wait
	return $result
}
