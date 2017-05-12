proc cg_qjobs {args} {
	set basedir [file_absolute [pwd]]
	set options {}
	cg_options qjobs args {
		-u {lappend options -u $value}
	} {} 0
	set xml [exec qstat -xml -pri {*}$options]
	set temp [string range $xml 21 end]
	regsub -all {[ \t]} $temp {} temp
	regsub -all {([ \n]*</?(job_|queue_info)[^>]*>[ \n]*)+} $temp # temp
	regsub -all {<([^ >]*)[^>]*>([^>]*)</([^>]+)>} $temp \\1\t\\2 temp
	set list [split $temp #]
	set result {}
	foreach el $list {
		set el [string trim $el]
		if {$el eq ""} continue
		set a(tasks) ""
		array set a [split $el \t\n]
		set resultline $a(JB_job_number),$a(tasks)
		set a(run) ?
		set a(runversion) ?
		regexp {^j([^.]+)\.([0-9_-]+)\.} $a(JB_name) temp a(run) a(runversion)
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
