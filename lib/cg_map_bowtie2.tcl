proc refseq_bowtie2_job {refseq} {
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

proc cg_refseq_bowtie2 args {
	set args [job_init {*}$args]
	set return [refseq_bowtie2_job {*}$args]
	job_wait
	return $return
}

proc refseq_bowtie2 {refseq} {
	upvar job_logdir job_logdir
	set refseq [file_absolute $refseq]
	set bowtie2refseq $refseq.bowtie2/[file tail $refseq]
	if {![file exists $bowtie2refseq]} {
		error "The bowtie2 version of the refseq does not exist (should be at $bowtie2refseq)
You can create it using:
cg refseq_bowtie2 \'$refseq\'"
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
		-x - -preset - -p {
			# not used
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
	set readgroupdata [map_readgroupdata $readgroupdata $sample]
	set resultbase [file root $result]
	set refseq [refseq $refseq]
	dbdir [file dir $refseq]
	set bowtie2refseq [refseq_bowtie2 $refseq]
	set outpipe [convert_pipe -.sam $result -endpipe 1 -refseq $refseq]
	analysisinfo_write $fastqfile1 $result aligner bowtie2 aligner_version [version bowtie2] reference [file2refname $bowtie2refseq] aligner_paired $paired
	set rg {}
	foreach {key value} $readgroupdata {
		lappend rg --rg "$key:$value"
	}
	if {$paired} {
		if {[expr {[llength $files]%2}]} {
			error "bowtie2 needs even number of files for paired analysis"
		}
		if {$fixmate} {
			set fixmate "| samtools fixmate -m -O sam - -"
		}
		set files1 {}
		set files2 {}
		foreach {file1 file2} $files {
			lappend files1 $file1
			lappend files2 $file2
		}
		exec bowtie2 -p 2 --sensitive -x $bowtie2refseq -1 [join $files1 ,] -2 [join $files2 ,] \
			--rg-id "$sample" {*}$rg \
			{*}$fixmate {*}$outpipe 2>@ stderr
	} else {
		exec bowtie2 -p 2 --sensitive -x $bowtie2refseq -U [join $files ,] \
			--rg-id "$sample" {*}$rg \
			{*}$outpipe 2>@ stderr
	}
}

proc cg_map_bowtie2 {args} {
	set args [job_init {*}$args]
	unset job_logdir
	set result [map_bowtie2_job {*}$args]
	job_wait
	return $result
}
