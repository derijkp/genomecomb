proc bowtie2refseq_job {refseq} {
	upvar job_logdir job_logdir
	set bowtie2refseq $refseq.bowtie2/[file tail $refseq]
	job bowtie2refseq-[file tail $refseq] -deps {$refseq} -targets {$refseq.bowtie2 $bowtie2refseq.1.bt2} \
	-vars {refseq} -code {
		file mkdir $refseq.bowtie2.temp
		mklink $refseq $refseq.bowtie2.temp/[file tail $refseq]
		exec bowtie2-build $refseq $refseq.bowtie2.temp/[file tail $refseq]
		file rename -force -- $refseq.bowtie2.temp $refseq.bowtie2
	}
	return $bowtie2refseq
}

proc map_bowtie2_job {args} {
	upvar job_logdir job_logdir
	set paired 1
	set readgroupdata {}
	set pre {}
	set skips {}
	set threads 2
	set fixmate 1
	cg_options map_bowtie2 args {
		-paired {
			set paired $value
		}
		-readgroupdata {
			set readgroupdata $value
		}
		-fixmate {
			set fixmate $value
		}
		-threads {
			set threads $value
			# not used (yet)
		}
		-pre {
			set pre $value
		}
		-skips {
			set skips $value
		}
		-aliformat {
			# currently ignored
			if {$value ne "bam"} {puts stderr "map_bowtie2 ignores aliformat, bam will always be made"}
			set aliformat $value
		}
	} {result refseq sample fastqfile1} 4 ... {
		align reads in fastq files to a reference genome using bowtie2
	}
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
			if {[expr {[llength $deps]%2}]} {
				error "bowtie2 needs even number of files for paired analysis"
			}
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
		file rename -force -- $temptarget $target
	}
	set analysisinfo [gzroot $result].analysisinfo
	job bowtie2_bam-$sample -deps {$resultbase.sam} \
	-targets {$result $analysisinfo} -vars {resultbase} {*}$skips \
	-code {
		puts "making $target"
		analysisinfo_write $dep $target fixmate samtools fixmate_version [version samtools]
		set tempfixmate [filetemp $target 1 1]
		set tempresult [filetemp $target 1 1]
		catch_exec samtools fixmate -O bam $resultbase.sam $tempfixmate
		cg_bam_sort $tempfixmate $tempresult
		file rename -force -- $tempresult $target
		file delete $tempfixmate $tempfixmate.analysisinfo
		file delete $resultbase.sam $resultbase.sam.analysisinfo
	}
	set indexext [indexext $result]
	job bowtie2_index-$sample -deps {$result} -targets {$result.$indexext} {*}$skips -code {
		exec samtools index $dep >@ stdout 2>@ stderr
		puts "making $target"
	}
}

proc cg_map_bowtie2 {args} {
	set args [job_init {*}$args]
	unset job_logdir
	set result [map_bowtie2_job {*}$args]
	job_wait
	return $result
}
