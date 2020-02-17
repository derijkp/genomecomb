proc distrreg_job {args} {
	upvar job_logdir job_logdir
	set maxopenfiles {}
	set sortfields {}
	set addcomment 1
	set refseq {}
	set threads 1
	cg_options distrreg args {
		-refseq {
			set refseq $value
		}
		-threads {
			set threads $value
		}
	} {file resultprefix resultsuffix regions} 4 4 {
		distribute (sorted) tsv files into multiple files per region
	}
	set targets {}
	foreach region $regions {
		lappend targets $resultprefix$region$resultsuffix
	}
	job distrreg-[file tail $file] -deps {$file} -targets $targets -vars {
		refseq threads file resultprefix resultsuffix regions
	} -code {
		cg_distrreg -refseq $refseq -threads $threads $file $resultprefix $resultsuffix $regions
	}
	return $targets
}

proc cg_distrreg {args} {
	# job_logdir
	set maxopenfiles {}
	set sortfields {}
	set addcomment 1
	set refseq {}
	set threads 1
	cg_options distrreg args {
		-refseq {
			set refseq $value
		}
		-threads {
			set threads $value
		}
	} {file resultprefix resultsuffix regions} 4 4 {
		distribute (sorted) tsv files into multiple files per region
	}
	set f [gzopen $file]
	set header [tsv_open $f]
	gzclose $f
	set poss [tsv_basicfields $header 3]
	set incmd [convert_pipe $file -.tsv -refseq $refseq -threads $threads]
	set outcmd [convert_pipe $file $resultprefix$resultsuffix -refseq $refseq -threads $threads]
	if {$incmd eq ""} {
		exec distrreg $resultprefix $resultsuffix 1 $regions {*}$poss 1 \# $outcmd < $file
	} else {
		exec $incmd | distrreg $resultprefix $resultsuffix 1 $regions {*}$poss 1 \# $outcmd
	}
}
