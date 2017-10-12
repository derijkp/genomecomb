proc bwarefseq_job {refseq} {
	upvar job_logdir job_logdir
	set bwarefseq $refseq.bwa/[file tail $refseq]
	if {[file exists $bwarefseq]} {return $bwarefseq}
	set tail [file tail $refseq]
	set targets [list $refseq.bwa $refseq.bwa/$tail]
	foreach ext {amb ann bwt pac sa} {lappend targets $refseq.bwa/$tail.$ext}
	if {[jobtargetexists $targets $refseq]} return
	job bwa2refseq-[file tail $refseq] -deps {$refseq} -targets $targets -code {
		file delete -force $target.temp
		file mkdir $target.temp
		mklink $dep $target.temp/[file tail $dep]
		exec bwa index $target.temp/[file tail $dep] 2>@ stderr
		file delete -force $target
		file rename $target.temp $target
	}
	return $bwarefseq
}

proc map_bwa_job {args} {
	upvar job_logdir job_logdir
	set paired 1
	set readgroupdata {}
	set pre {}
	set skips {}
	set threads 2
	cg_options map_bwa args {
		-paired - -p {
			set paired $value
		}
		-readgroupdata {
			set readgroupdata $value
		}
		-pre {
			set pre $value
		}
		-skips {
			set skips $value
		}
		-threads - -t {
			set threads $value
		}
	} {result refseq sample fastqfile1} 4
	set files [list $fastqfile1 {*}$args]
	if {![info exists job_logdir]} {
		job_logdir [file dir $result]/log_jobs
	}
	array set a [list PL illumina LB solexa-123 PU $sample SM $sample]
	if {$readgroupdata ne ""} {
		array set a $readgroupdata
	}
	set readgroupdata [array get a]
	set bwarefseq [bwarefseq_job $refseq]
	set resultbase [file root $result]
	set samfiles {}
	set num 1
	if {!$paired} {
		foreach file $files {
			set name [file root [file tail $file]]
			set target $resultbase-$name.sam
			lappend samfiles $target
			job bwa-$sample-$name -mem 5G -cores $threads \
			-deps [list $bwarefseq $file] -targets {$target} -vars {readgroupdata sample paired threads} \
			-skip [list $resultbase.bam] {*}$skips -code {
				puts "making $target"
				foreach {bwarefseq fastq} $deps break
				set rg {}
				foreach {key value} $readgroupdata {
					lappend rg "$key:$value"
				}
				exec bwa mem -t $threads -R @RG\\tID:$sample\\t[join $rg \\t] $bwarefseq $fastq > $target.temp 2>@ stderr
				file rename -force $target.temp $target
			}
		}
	} else {
		foreach {file1 file2} $files {
			set name [file root [file tail $file1]]
			set target $resultbase-$name.sam
			lappend samfiles $target
			job bwa-$sample-$name -mem 5G -cores $threads -deps [list $bwarefseq $file1 $file2] -targets {$target} -vars {readgroupdata sample paired threads} \
			-skip [list $resultbase.bam] {*}$skips -code {
				puts "making $target"
				foreach {bwarefseq fastq1 fastq2} $deps break
				set rg {}
				foreach {key value} $readgroupdata {
					lappend rg "$key:$value"
				}
				exec bwa mem -t $threads -M -R @RG\\tID:$sample\\t[join $rg \\t] $bwarefseq $fastq1 $fastq2 > $target.temp 2>@ stderr
				file rename -force $target.temp $target
				file delete bwa1.fastq bwa2.fastq
			}
		}
	}
	job bwa2bam-$sample -deps $samfiles -rmtargets $samfiles -targets {$result $result.bai} {*}$skips -vars {resultbase} -code {
		puts "making $target"
			if {[catch {
				exec samcat {*}$deps | bamsort SO=coordinate tmpfile=[scratchfile] index=1 indexfilename=$target.bai inputformat=sam > $target.temp 2>@ stderr
			}]} {
				error $msg
			}
			file rename -force $target.temp $target
		file delete {*}$deps
	}
}

proc cg_map_bwa {args} {
	set args [job_init {*}$args]
	unset job_logdir
	map_bwa_job {*}$args
	job_wait
}
