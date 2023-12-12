proc version_STAR {} {
	set version [exec STAR --version]
	regsub -- {_20..-..-..$} $version {} version
	return $version
}

proc version_star {} {
	version_STAR
}

proc refseq_star_job {refseq {gtffile {}} {threads 4}} {
	upvar job_logdir job_logdir
	set refseq [refseq $refseq]
	set gtffile [file_absolute $gtffile]
	set starrefseq $refseq.star
	if {[file exists $starrefseq]} {return $starrefseq}
	if {[jobtargetexists [list $starrefseq] $refseq]} return
	job [job_relfile2name star_2refseq- $refseq] -deps {
		$refseq
	} -targets {
		$starrefseq
	} -vars {
		refseq starrefseq gtffile threads
	} -code {
		if {[gziscompressed $gtffile]} {
			set tempfile [tempfile].gtf
			exec cg zcat $gtffile > $tempfile
			set gtffile $tempfile
		}
		file delete -force $starrefseq.temp
		file mkdir $starrefseq.temp
		set tail [file tail $refseq]
		mklink $refseq $starrefseq.temp/$tail
		set extraopts {}
		if {$gtffile ne ""} {
			lappend extraopts --sjdbGTFfile $gtffile
		}
		set temp [catch_exec STAR \
			--runMode genomeGenerate \
			--genomeFastaFiles $refseq \
			--genomeDir $starrefseq.temp \
			--runThreadN $threads \
			{*}$extraopts \
		]
		if {[regexp {loaded/built the index for 0 target sequence\(s\)} $temp]} {
			error "could not properly index $dep: contains no sequences"
		}
		file rename -- $target.temp $target
	}
	return $starrefseq
}

proc cg_refseq_star args {
	set args [job_init {*}$args]
	set return [refseq_star_job {*}$args]
	job_wait
	return $return
}

proc refseq_star {refseq} {
	upvar job_logdir job_logdir
	set refseq [file_absolute $refseq]
	set starrefseq $refseq.star
	if {![file exists $starrefseq]} {
		error "The star index of the refseq does not exist (should be at $starrefseq)
You can create it using:
cg refseq_star \'$refseq\'"
	}
	return $starrefseq
}

proc map_mem_star {mem threads preset deps} {
	return 30G
}

proc cg_map_star {args} {
	set extraopts {}
	set paired 1
	set keepargs $args
	set readgroupdata {}
	set threads 2
	set mem 30G
	set fixmate 1
	set aliformat bam
	set ali_keepcomments {}
	set preset {}
	cg_options map_star args {
		-paired - -p {
			set paired $value
		}
		-preset {
			if {$value eq "2p"} {
				lappend extraopts --twopassMode Basic
			}
			set preset $value
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
		-ali_keepcomments {
			# not used yet
			set ali_keepcomments [true $value]
		}
		-extraopts {
			set extraopts $value
		}
		-mem {
			set mem $value
		}
	} {result refseq sample fastqfile1} 4 5 {
		align reads in fastq files to a reference genome using star
	}
	set opts {}
	foreach {key value} [specialopts -star] {
		switch $key {
			default {
				if {[string index $key 1] != "-"} {
					set key -$key
				}
				lappend extraopts $key $value
			}
		}
	}
	set files [list $fastqfile1 {*}$args]
	if {$paired && [llength $files] != 2} {
		error "cg map_star error: for paired read alignment 2 fastq files must be given"
	}
	if {!$paired && [llength $files] != 1} {
		error "cg map_star error: for single read alignment 1 fastq files must be given"
	}
	set result [file_absolute $result]
	set refseq [refseq $refseq]
	#
	set readgroupdata [map_readgroupdata $readgroupdata $sample]
	set starrefseq [refseq_star $refseq]
	set outpipe [convert_pipe -.sam $result -endpipe 1 -refseq $refseq]
	analysisinfo_write $fastqfile1 $result sample [file tail $sample] aligner star aligner_version [version star] reference [file2refname $starrefseq] aligner_paired $paired
	putslog "making $result"
	set rg {}
	foreach {key value} [sam_readgroupdata_fix $readgroupdata] {
		lappend rg "$key:$value"
	}
	if {!$paired} {
		set files1 {}
		foreach {file} $files {
			lappend rgids ID:$sample\ [join $rg " "]
			lappend files1 [file_absolute $file]
		}
		set keepdir [pwd]
		set scratch [scratchdir]/STAR
		catch {file delete -force $scratch}
		mkdir $scratch
		cd $scratch
		if {[catch {
			exec STAR \
				--outTmpDir $scratch/STARTMP \
				--runThreadN $threads \
				--genomeDir $starrefseq \
				--readFilesIn [join $files1,] \
				--readFilesCommand {*}[gzcat [lindex $files1 0]] \
				--outSAMattrRGline {*}[join $rgids " , "] \
				--outSAMunmapped Within \
				--outStd SAM \
				{*}$extraopts \
				{*}$outpipe
		} msg]} {
			cd $keepdir
			if {[regexp ERROR: $msg] || $::errorCode ne "NONE"} {
				puts stderr $msg
				error $msg
			}
		}
		cd $keepdir
		file delete -force $scratch
		puts stderr $msg
	} else {
		if {$fixmate} {
			set fixmate "| samtools fixmate -m -O sam - -"
		}
		if {[expr {[llength $files]%2}]} {
			error "star needs even number of files for paired analysis"
		}
		set files1 {}
		set files2 {}
		set rgids {}
		foreach {file1 file2} $files {
			lappend files1 [file_absolute $file1]
			lappend files2 [file_absolute $file2]
			lappend rgids ID:$sample\ [join $rg " "]
		}
		putslog "making $result"
		set keepdir [pwd]
		set scratch [scratchdir]/STAR
		catch {file delete -force $scratch}
		mkdir $scratch
		cd $scratch
		if {[catch {
			exec STAR \
				--outTmpDir $scratch/STARTMP \
				--runThreadN $threads \
				--genomeDir $starrefseq \
				--readFilesIn [join $files1 ,] [join $files2 ,] \
				--readFilesCommand {*}[gzcat [lindex $files1 0]] \
				--outSAMattrRGline {*}[join $rgids " , "] \
				--outSAMunmapped Within \
				--outStd SAM \
				{*}$extraopts \
				{*}$outpipe
		} msg]} {
			cd $keepdir
			if {[regexp ERROR: $msg] || $::errorCode ne "NONE"} {
				puts stderr $msg
				error $msg
			}
		}
		cd $keepdir
		file delete -force $scratch
		puts stderr $msg
	}
}
