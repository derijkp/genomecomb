proc bam2covstats_job {bamfile regionfile {suffix {}}} {
	upvar job_logdir job_logdir
	set bamfile [file_absolute $bamfile]
	set dir [file dir $bamfile]
	set file [file tail $bamfile]
	set root [join [lrange [split [file root $file] -] 1 end] -]
	if {$root eq ""} {set root $file}
#	job bam2coverage-$root -deps $bamfile -targets {$dir/coverage-$root $dir/coverage-$root/coverage-$root.FINISHED} -vars {root} -code {
#		cg bam2coverage $dep $target/coverage-$root
#	}
	job make_histo-$root -deps {$bamfile $bamfile.bai $regionfile} -targets $dir/$root.histo${suffix} -vars {regionfile} -code {
		set tempfile [filetemp $target]
		cg bam_histo $regionfile $dep {1 5 10 20 50 100 200 500 1000} > $tempfile
		file rename -force $tempfile $target
	}
}

proc cg_coverage_report {args} {
	set args [job_init {*}$args]
	set suffix {}
	cg_options coverage_report args {
		-s - --suffix {
			set suffix $value
		}
	} {regionfile bamfile} 1
	swet bams [list $bamfile {*}$args]
	foreach b $bams {
		bam2covstats_job $b $regionfile $suffix
	}

	#job coverage_report-$experiment -deps [list $regfile {*}$histofiles] -targets [list coverage_${experiment}_avg.tsv coverage_${experiment}_frac_above_50.tsv ] -code {
	#	exec python2.6 /complgen2/mastr-procedure/coverage_mastrs.py
	#}

	job_wait
}
