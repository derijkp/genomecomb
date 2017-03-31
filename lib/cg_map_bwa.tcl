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
	job bwa-$sample -mem 5G -deps [list $bwarefseq {*}$files] -targets {$resultbase.sam} -vars {readgroupdata sample paired} \
	-skip $resultbase.bam {*}$skips -code {
		puts "making $target"
		set bwarefseq [list_shift deps]
		set rg {}
		foreach {key value} $readgroupdata {
			lappend rg "$key:$value"
		}
		set tempdir [scratchdir]
		if {!$paired} {
			if {[llength $deps] > 1} {
				exec cat {*}$deps > $tempdir/bwa1.fastq
			} else {
				mklink $deps $tempdir/bwa1.fastq
			}
			exec bwa mem -t 2 -R @RG\\tID:$sample\\t[join $rg \\t] $bwarefseq $tempdir/bwa1.fastq > $target.temp 2>@ stderr
		} else {
			set files1 {}
			set files2 {}
			foreach {file1 file2} $deps {
				lappend files1 $file1
				lappend files2 $file2
			}
			if {[llength $files1] > 1} {
				exec cat {*}$files1 > $tempdir/bwa1.fastq
				exec cat {*}$files2 > $tempdir/bwa2.fastq
			} else {
				mklink $file1 $tempdir/bwa1.fastq
				mklink $file2 $tempdir/bwa2.fastq
			}
			exec bwa mem -t 2 -M -R @RG\\tID:$sample\\t[join $rg \\t] $bwarefseq $tempdir/bwa1.fastq $tempdir/bwa2.fastq > $target.temp 2>@ stderr
		}
		file rename -force $target.temp $target
		file delete bwa1.fastq bwa2.fastq
	}
	job bwa2bam-$sample -deps $resultbase.sam -targets $result {*}$skips -vars {resultbase} -code {
		puts "making $target"
		exec samtools view -S -h -b -o $resultbase.ubam $resultbase.sam >@ stdout 2>@ stderr
		samtools_sort $resultbase.ubam $target.temp
		file rename -force $target.temp $target
		file delete $resultbase.ubam
		file delete $resultbase.sam
	}
	job bwa_index-$sample -deps $result -targets $result.bai {*}$skips -code {
		exec samtools index $dep >@ stdout 2>@ stderr
		puts "making $target"
	}
}

