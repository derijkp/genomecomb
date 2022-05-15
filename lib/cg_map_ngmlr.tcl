proc refseq_ngmlr_job {refseq preset} {
	upvar job_logdir job_logdir
	set tail [file tail $refseq]
	set ngmlrrefseqdir $refseq.ngmlr.$preset
	set ngmlrrefseq $ngmlrrefseqdir/$tail
	if {[jobtargetexists $ngmlrrefseq]} {
		return $ngmlrrefseq
	}
	if {![jobtargetexists [list $ngmlrrefseq-ht-13-2.2.ngm $refseq-enc.2.ngm] [list $ngmlrrefseq]]} {
		job [job_relfile2name ngmlr_2refseq- $refseq] -deps {
			$refseq
		} -targets {
			$ngmlrrefseq
		} -vars {
			refseq ngmlrrefseqdir ngmlrrefseq preset
		} -code {
			file mkdir $ngmlrrefseqdir.temp
			set tail [file tail $refseq]
			mklink $refseq $ngmlrrefseqdir.temp/$tail
			set temp [tempfile]
			file_write $temp ""
			if {[catch {
				exec -ignorestderr ngmlr -x $preset -r $ngmlrrefseqdir.temp/$tail -q $temp 2>@ stderr
			} e]} {
				error $e
			}
			file rename $ngmlrrefseqdir.temp $ngmlrrefseqdir
		}
	}
	return $ngmlrrefseq
}

proc cg_refseq_ngmlr args {
	set args [job_init {*}$args]
	set return [refseq_ngmlr_job {*}$args]
	job_wait
	return $return
}

proc refseq_ngmlr {refseq preset} {
	upvar job_logdir job_logdir
	set refseq [file_absolute $refseq]
	set ngmlrrefseq $refseq.ngmlr.$preset/[file tail $refseq]
	if {![file exists $ngmlrrefseq]} {
		error "The ngmlr version of the refseq (preset $preset) does not exist (should be $ngmlrrefseq)
You can create it using:
cg refseq_ngmlr \'$refseq\' $preset"
	}
	return $ngmlrrefseq
}

proc cg_map_ngmlr {args} {
	upvar job_logdir job_logdir
	set paired 1
	set keepargs $args
	set preset {}
	set readgroupdata {}
	set threads 2
	set aliformat bam
	cg_options map_ngmlr args {
		-paired {
			if {$value} {error "ngmlr does not support paired read alignment"}
		}
		-x - -preset - -p {
			set preset $value
		}
		-readgroupdata {
			set readgroupdata $value
		}
		-fixmate {
			# not used
		}
		-threads - -t {
			set threads $value
		}
	} {result refseq sample fastqfile} 4 4 {
		align reads in fastq files to a reference genome using ngmlr
	}
	if {$preset eq ""} {set preset ont}
	set result [file_absolute $result]
	set refseq [refseq $refseq]
	set readgroupdata [map_readgroupdata $readgroupdata $sample]
	set ngmlr_refseq [refseq_ngmlr $refseq $preset]
	set outpipe [convert_pipe -.sam $result -endpipe 1 -refseq $refseq]
	putslog "making $result"
	analysisinfo_write $fastqfile $result aligner ngmlr aligner_version [version ngmlr] reference [file2refname $ngmlr_refseq] aligner_paired 0
	set rg {}
	foreach {key value} [sam_readgroupdata_fix $readgroupdata] {
		lappend rg "$key:$value"
	}
	exec ngmlr -x $preset -t $threads -r $ngmlr_refseq -q $fastqfile {*}$outpipe 2>@ stderr
}
