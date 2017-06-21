proc bwarefseq_job {refseq} {
	upvar job_logdir job_logdir
	set bwarefseq $refseq.bwa/[file tail $refseq]
	if {[file exists $bwarefseq]} {return $bwarefseq}
	set tail [file tail $refseq]
	set targets [list $refseq.bwa $refseq.bwa/$tail]
	foreach ext {amb ann bwt pac sa} {lappend targets $refseq.bwa/$tail.$ext}
	if {[jobtargetexists $targets $refseq]} return
	job bwa2refseq-[file tail $refseq] -deps $refseq -targets $targets -code {
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
	oargs map_bwa_job {result refseq files sample
		{paired 1}
		{readgroupdata {}}
		{pre {}}
		{skips {}}
	} $args
	upvar job_logdir job_logdir
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
			job bwa-$sample-$name -mem 5G -deps [list $bwarefseq $file] -targets {$target} -vars {readgroupdata sample paired} \
			-skip $resultbase.bam {*}$skips -code {
				puts "making $target"
				foreach {bwarefseq fastq} $deps break
				set rg {}
				foreach {key value} $readgroupdata {
					lappend rg "$key:$value"
				}
				exec bwa mem -t 2 -R @RG\\tID:$sample\\t[join $rg \\t] $bwarefseq $fastq > $target.temp 2>@ stderr
				file rename -force $target.temp $target
			}
		}
	} else {
		foreach {file1 file2} $files {
			set name [file root [file tail $file1]]
			set target $resultbase-$name.sam
			lappend samfiles $target
			job bwa-$sample-$name -mem 5G -deps [list $bwarefseq $file1 $file2] -targets {$target} -vars {readgroupdata sample paired} \
			-skip $resultbase.bam {*}$skips -code {
				puts "making $target"
				foreach {bwarefseq fastq1 fastq2} $deps break
				set rg {}
				foreach {key value} $readgroupdata {
					lappend rg "$key:$value"
				}
				exec bwa mem -t 2 -M -R @RG\\tID:$sample\\t[join $rg \\t] $bwarefseq $fastq1 $fastq2 > $target.temp 2>@ stderr
				file rename -force $target.temp $target
				file delete bwa1.fastq bwa2.fastq
			}
		}
	}
	job bwa2bam-$sample -deps $samfiles -rmtargets $samfiles -targets $result {*}$skips -vars {resultbase} -code {
		puts "making $target"
		if {[catch {version samtools 1}]} {
			if {[catch {
				exec samcat {*}$deps | samtools view -S -h -b - | samtools sort - $target.temp 2>@ stderr
			}]} {
				error $msg
			}
			file rename -force $target.temp.bam $target
		} else {
			if {[catch {
				exec samcat {*}$deps | samtools view -S -h -b - | samtools sort - > $target.temp 2>@ stderr
			}]} {
				error $msg
			}
			file rename -force $target.temp $target
		}
		file delete {*}$deps
	}
	job bwa_index-$sample -deps $result -targets $result.bai {*}$skips -code {
		exec samtools index $dep >@ stdout 2>@ stderr
		puts "making $target"
	}
}

