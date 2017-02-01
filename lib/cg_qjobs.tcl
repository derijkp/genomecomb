proc cg_qjobs {args} {
	set basedir [file_absolute [pwd]]
	set options {}
	cg_options qjobs args {
	} {} 0 0	
	set xml [exec qstat -xml]
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
		foreach field {JB_job_number tasks state JB_name} {
			lappend resultline $a($field)
		}
		lappend result $resultline
	}
	set result [lsort -dictionary -index 0 $result]
	foreach line $result {
		puts [join [lrange $line 1 end] \t]
	}
}
