proc cg_qjobs {args} {
	set basedir [file_absolute [pwd]]
	set options {}
	cg_options qjobs args {
		-u {lappend options -u $value}
	} {} 0
	set result {}
	set xml [exec qstat -xml -pri {*}$options]
	set list [regexp -all -inline {<job_list .*?</job_list>} $xml]
	foreach el $list {
		set data [regexp -all -inline {<([^>]+)>([^>]*?)</[^>]+>} $el]
		set a(tasks) ""
		set a(run) ?
		set a(runversion) ?
		foreach {temp key value} $data {
			set a($key) $value
		}
		set resultline $a(JB_job_number),$a(tasks)
		if {[regexp {^j([^#]+)\.([0-9_-]+)\#} $a(JB_name) temp a(run) a(runversion)]} {
		} elseif {[regexp {^j([^.]+)\.([0-9_-]+)\.} $a(JB_name) temp a(run) a(runversion)]} {
		} else {
			regexp {^j(.+)\.([0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]_[0-9][0-9]-[0-9][0-9]-[0-9][0-9]-[0-9][0-9][0-9])\.} $a(JB_name) temp a(run) a(runversion)
		}
		foreach field {JB_job_number tasks state JB_submission_time JAT_start_time JB_priority JAT_prio JB_owner queue_name slots run runversion JB_name} {
			lappend resultline [get a($field) .]
		}
		lappend result $resultline
	}
	set result [lsort -dictionary -index 0 $result]
	puts [join {id tasks state submissiontime starttime priority JAT_prio owner queue slots run runversion name} \t]
	foreach line $result {
		puts [join [lrange $line 1 end] \t]
	}
}
