proc refseq_ngmlr_job {refseq preset} {
	upvar job_logdir job_logdir
	if {![jobtargetexists [list $refseq-ht-13-2.2.ngm $refseq-enc.2.ngm] [list $refseq]]} {
		set tail [file tail $refseq]
		job ngmlr_2refseq-[file_part $refseq end] -deps {$refseq} -targets {
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

proc cg_refseq_ngmlr args {
	set args [job_init {*}$args]
	set return [refseq_ngmlr_job {*}$args]
	job_wait
	return $return
}

proc refseq_ngmlr {refseq preset} {
	upvar job_logdir job_logdir
	set refseq [file_absolute $refseq]
	if {![file exists $refseq] || ![file exists $refseq-ht-13-2.2.ngm] || ![file exists $refseq-enc.2.ngm]} {
		error "The ngmlr version of the refseq does not exist (should be at $refseq-ht-13-2.2.ngm and $refseq-enc.2.ngm)
You can create it using:
cg refseq_ngmlr \'$refseq\' $preset"
	}
	return $refseq
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
	foreach {key value} $readgroupdata {
		lappend rg "$key:$value"
	}
	exec ngmlr -x $preset -t $threads -r $ngmlr_refseq -q $fastqfile {*}$outpipe 2>@ stderr
}
