proc distrreg_job {args} {
	upvar job_logdir job_logdir
	set maxopenfiles {}
	set sortfields {}
	set addcomment 1
	set refseq {}
	set threads 1
	set other {}
	set skips {}
	cg_options distrreg args {
		-refseq {
			set refseq $value
		}
		-threads {
			set threads $value
		}
		-other {
			set other $value
		}
		-skip {
			lappend skips -skip $value
		}
	} {file resultprefix resultsuffix regions} 4 4 {
		distribute (sorted) tsv files into multiple files per region
	}
	set targets {}
	foreach region $regions {
		lappend targets $resultprefix$region$resultsuffix
	}
	job [job_relfile2name distrreg- $file] {*}$skips -deps {$file} -targets $targets -vars {
		refseq threads file resultprefix resultsuffix regions other
	} -code {
		cg_distrreg -refseq $refseq -threads $threads -other $other $file $resultprefix $resultsuffix $regions
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
	set other {}
	cg_options distrreg args {
		-refseq {
			set refseq $value
		}
		-threads {
			set threads $value
		}
		-other {
			set other $value
		}
	} {file resultprefix resultsuffix regions} 4 4 {
		distribute (sorted) tsv files into multiple files per region
	}
	foreach region $regions {
		set target $resultprefix$region$resultsuffix
		analysisinfo_write $file $target
	}
	set f [gzopen $file]
	set header [tsv_open $f]
	gzclose $f
	set poss [tsv_basicfields $header 3]
	set incmd [convert_pipe $file -.tsv -refseq $refseq -threads $threads]
	set outcmd [convert_pipe -.tsv $resultprefix$resultsuffix -refseq $refseq -threads $threads]
	if {$incmd eq ""} {
		# puts [list distrreg $resultprefix $resultsuffix 1 $regions {*}$poss 1 \# $outcmd < $file]
		catch_exec distrreg $resultprefix $resultsuffix 1 $regions {*}$poss 1 \# $outcmd $other < $file 2>@ stderr
	} else {
		# puts [list {*}$incmd | distrreg $resultprefix $resultsuffix 1 $regions {*}$poss 1 \# $outcmd $other]
		catch_exec {*}$incmd | distrreg $resultprefix $resultsuffix 1 $regions {*}$poss 1 \# $outcmd $other 2>@ stderr
	}
}
