proc bam2reg_job {args} {
	upvar job_logdir job_logdir
	set compress 1
	set skips {}
	cg_options bam2reg args {
		-mincoverage {
			set mincov $value
		}
		-compress {
			set compress $value
		}
		-skip {
			lappend skips -skip $value
		}
	} {bamfile mincoverage target} 1 3
	set bamfile [file_absolute $bamfile]
	set pre [lindex [split $bamfile -] 0]
	set dir [file dir $bamfile]
	set file [file tail $bamfile]
	set root [file_rootname $file]
	if {![info exists target] && [info exists mincoverage] && [info exists mincov]} {
		set target $mincoverage
		set mincoverage $mincov
	}
	if {![info exists mincoverage]} {
		if {[info exists mincov]} {set mincoverage $mincov} else {set mincoverage 5}
	}
	if {![info exists target]} {
		set target $dir/sreg-cov$mincoverage-$root.tsv
		if {$compress} {append target .lz4}
	}
	if {![info exists job_logdir]} {
		job_logdir $target.log_jobs
	}
	job cov$mincoverage-$root -optional 1 {*}$skips -deps {
		$bamfile
	} -targets {
		$target
	} -vars {
		mincoverage
	} -code {
		set compress [compresspipe $target]
		set temptarget [filetemp $target]
		exec cg regextract -min $mincoverage $dep {*}$compress > $temptarget
		file rename -force $temptarget $target
		if {[file extension $target] eq ".lz4"} {cg_lz4index $target}
	}
	return $target
}

proc cg_bam2reg {args} {
	set args [job_init {*}$args]
	unset job_logdir
	set result [bam2reg_job {*}$args]
	job_wait
	return $result
}
