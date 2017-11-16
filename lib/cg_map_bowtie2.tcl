proc bowtie2refseq_job {refseq} {
	upvar job_logdir job_logdir
	set bowtie2refseq $refseq.bowtie2/[file tail $refseq]
	job bowtie2refseq-[file tail $refseq] -deps {$refseq} -targets {$refseq.bowtie2 $bowtie2refseq.1.bt2} \
	-vars {refseq} -code {
		file mkdir $refseq.bowtie2.temp
		mklink $refseq $refseq.bowtie2.temp/[file tail $refseq]
		exec bowtie2-build $refseq $refseq.bowtie2.temp/[file tail $refseq]
		file rename -force $refseq.bowtie2.temp $refseq.bowtie2
	}
	return $bowtie2refseq
}

proc map_bowtie2_job {args} {
	upvar job_logdir job_logdir
	set paired 1
	set readgroupdata {}
	set pre {}
	set skips {}
	cg_options map_bowtie2 args {
		-paired {
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
	} {result refseq sample fastqfile1} 4
	set files [list $fastqfile1 {*}$args]
	if {![info exists job_logdir]} {
		job_logdir [file dir $result]/log_jobs
	}
	array set a [list PL illumina LB solexa-123 PU $sample SM $sample]
	if {$readgroupdata ne ""} {
		array set a $readgroupdata
	}
	set resultbase [file root $result]
	set readgroupdata [array get a]
	upvar job_logdir job_logdir
	set bowtie2refseq [bowtie2refseq_job $refseq]
	job bowtie2-$sample -deps [list $bowtie2refseq {*}$files] \
	-targets {$resultbase.sam $resultbase.sam.analysisinfo} \
	-vars {paired bowtie2refseq readgroupdata sample} \
	-skip [list $resultbase.bam] {*}$skips -code {
		puts "making $target"
		analysisinfo_write $dep2 $target aligner bowtie2 aligner_version [version bowtie2] reference [file2refname $bowtie2refseq] aligner_paired $paired
		list_shift deps
		set rg {}
		foreach {key value} $readgroupdata {
			lappend rg --rg "$key:$value"
		}
		set temptarget [filetemp $target]
		if {$paired} {
			set files1 {}
			set files2 {}
			foreach {file1 file2} $deps {
				lappend files1 $file1
				lappend files2 $file2
			}
			exec bowtie2 -p 2 --sensitive -x $bowtie2refseq -1 [join $files1 ,] -2 [join $files2 ,] \
			--rg-id "$sample" {*}$rg \
			-S $temptarget >@ stdout 2>@ stderr
		} else {
			exec bowtie2 -p 2 --sensitive -x $bowtie2refseq -U [join $$deps ,] \
			--rg-id "$sample" {*}$rg \
			-S $temptarget >@ stdout 2>@ stderr
		}
		file rename -force $temptarget $target
	}
	set analysisinfo [gzroot $result].analysisinfo
	job bowtie2_bam-$sample -deps {$resultbase.sam} \
	-targets {$result $analysisinfo} -vars {resultbase} {*}$skips \
	-code {
		puts "making $target"
		analysisinfo_write $dep $target
		catch {exec samtools view -S -h -b -o $resultbase.ubam $resultbase.sam >@ stdout 2>@ stderr}
		catch {bam_sort $resultbase.ubam $target.temp}
		file rename -force $target.temp $target
		file delete $resultbase.ubam $resultbase.ubam.analysisinfo
		file delete $resultbase.sam $resultbase.sam.analysisinfo
	}
	job bowtie2_index-$sample -deps {$result} -targets {$result.bai} {*}$skips -code {
		exec samtools index $dep >@ stdout 2>@ stderr
		puts "making $target"
	}
}

proc cg_map_bowtie2 {args} {
	set args [job_init {*}$args]
	unset job_logdir
	map_bowtie2_job {*}$args
	job_wait
}
