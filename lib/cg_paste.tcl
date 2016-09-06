proc tsv_paste_job {target files} {
	set target [file_absolute $target]
	job_logdir [file dir $target]/log_jobs
	job paste-[file tail $target] -deps $files -targets {$target} -code {
		exec tsv_paste {*}$deps > $target.temp 2>@ stderr
		file rename -force $target.temp $target
	}
}

proc cg_paste {args} {
	set args [job_init -silent {*}$args]
	cg_options paste args {
		-o - --outputfile {
			set outputfile $value
		}
	} 1
	if {[info exists outputfile]} {
		tsv_paste_job $outputfile $args
	} else {
		exec tsv_paste {*}$args >@ stdout 2>@ stderr
	}
}

