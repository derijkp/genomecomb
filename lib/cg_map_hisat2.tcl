proc refseq_hisat2_job {refseq} {
	upvar job_logdir job_logdir
	set hisat2refseq $refseq.hisat2
	if {[file exists $hisat2refseq]} {return $hisat2refseq}
	if {[jobtargetexists [list $hisat2refseq] $refseq]} return
	job [job_relfile2name hisat2_2refseq- $refseq] -deps {
		$refseq
	} -targets {
		$hisat2refseq
	} -vars {
		refseq hisat2refseq
	} -code {
		file delete $hisat2refseq.temp
		file mkdir $hisat2refseq.temp
		set tail [file tail $refseq]
		mklink $refseq $hisat2refseq.temp/$tail
		set temp [catch_exec hisat2-build $refseq $hisat2refseq.temp/$tail]
		if {[regexp {loaded/built the index for 0 target sequence\(s\)} $temp]} {
			error "could not properly index $dep: contains no sequences"
		}
		file rename -- $target.temp $target
	}
	return $hisat2refseq
}

proc cg_refseq_hisat2 args {
	set args [job_init {*}$args]
	set return [refseq_hisat2_job {*}$args]
	job_wait
	return $return
}

proc refseq_hisat2 {refseq} {
	upvar job_logdir job_logdir
	set refseq [file_absolute $refseq]
	set hisat2refseq $refseq.hisat2/[file tail $refseq]
	if {![file exists $hisat2refseq]} {
		error "The hisat2 index of the refseq does not exist (should be at $hisat2refseq)
You can create it using:
cg refseq_hisat2 \'$refseq\'"
	}
	return $hisat2refseq
}

proc map_mem_hisat2 {mem threads preset deps} {
	return 5G
}

proc cg_map_hisat2 {args} {
	if {[info exists ::cgextraopts(hisat2)]} {set extraopts $::cgextraopts(hisat2)} else {set extraopts {}}
	set paired 1
	set keepargs $args
	set readgroupdata {}
	set threads 2
	set mem 5G
	set fixmate 1
	set aliformat bam
	cg_options map_hisat2 args {
		-paired - -p {
			set paired $value
		}
		-readgroupdata {
			set readgroupdata $value
		}
		-fixmate {
			set fixmate $value
		}
		-threads - -t {
			set threads $value
		}
		-mem {
			set mem $value
		}
		-extraopts {
			lappend extraopts {*}$value
		}
	} {result refseq sample fastqfile1} 4 5 {
		align reads in fastq files to a reference genome using hisat2
	}
	set files [list $fastqfile1 {*}$args]
	set result [file_absolute $result]
	set refseq [refseq $refseq]
	#
	set readgroupdata [map_readgroupdata $readgroupdata $sample]
	set hisat2refseq [refseq_hisat2 $refseq]
	set outpipe [convert_pipe -.sam $result -endpipe 1 -refseq $refseq]
	analysisinfo_write $fastqfile1 $result aligner hisat2 aligner_version [version hisat2] reference [file2refname $hisat2refseq] aligner_paired $paired
	if {!$paired} {
		putslog "making $result"
		set rg {}
		foreach {key value} [sam_readgroupdata_fix $readgroupdata] {
			lappend rg "$key:$value"
		}
		if {[catch {
			exec hisat2 -t $threads \
				--rg-id @RG\\tID:$sample\\t[join $rg \\t] \
				-x $hisat2refseq -U [join $files ,] \
				{*}$extraopts \
				{*}$outpipe
		} msg]} {
			if {[regexp ERROR: $msg] || $::errorCode ne "NONE"} {
				puts stderr $msg
				error $msg
			}
		}
		puts stderr $msg
	} else {
		if {$fixmate} {
			set fixmate "| samtools fixmate -m -O sam - -"
		}
		if {[expr {[llength $files]%2}]} {
			error "hisat2 needs even number of files for paired analysis"
		}
		set files1 {}
		set files2 {}
		foreach {file1 file2} $files {
			lappend files1 $file1
			lappend files2 $file2
		}
		putslog "making $result"
		set rg {}
		foreach {key value} [sam_readgroupdata_fix $readgroupdata] {
			lappend rg "$key:$value"
		}
		if {[catch {
			exec hisat2 --threads $threads \
				--rg-id @RG\\tID:$sample\\t[join $rg \\t] \
				-x $hisat2refseq \
				-1 [join $files1 ,] -2 [join $files2 ,] \
				{*}$extraopts \
				{*}$fixmate {*}$outpipe
		} msg]} {
			if {[regexp ERROR: $msg] || $::errorCode ne "NONE"} {
				puts stderr $msg
				error $msg
			}
		}
		puts stderr $msg
	}
}
