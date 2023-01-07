proc cg_qjobs {args} {
	set basedir [file_absolute [pwd]]
	set options {}
	set user {}
	cg_options qjobs args {
		-u {set user $value}
	} {} 0 0 {
		returns running and waiting jobs on a grid engine cluster in tsv format (so they can by analysed using cg select)
	}
	set result {}
	if {![catch {
		exec which qstat
	}]} {
		if {$user ne ""} {
			lappend options -u $value
		}
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
			} elseif {[regexp {^j([^#]+)\.([0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]_[0-9][0-9-]+)\.\.\.\.} $a(JB_name) temp a(run) a(runversion)]} {
			} elseif {[regexp {^j(.+)\.([0-9_-]+)$} [file root $a(JB_name)] temp a(run) a(runversion)]} {
			} elseif {[regexp {^j([^.]+)\.([0-9_-]+)\.} $a(JB_name) temp a(run) a(runversion)]} {
			} else {
				regexp {^j(.+)\.([0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]_[0-9][0-9]-[0-9][0-9]-[0-9][0-9]-[0-9][0-9][0-9])\.} $a(JB_name) temp a(run) a(runversion)
			}
			foreach field {JB_job_number tasks state JB_submission_time JAT_start_time JB_priority JAT_prio JB_owner queue_name slots run runversion JB_name} {
				lappend resultline [get a($field) .]
			}
			lappend result $resultline
		}
		set result [bsort -index 0 $result]
		puts [join {id tasks state submissiontime starttime priority JAT_prio owner queue slots run runversion name} \t]
		foreach line $result {
			puts [join [lrange $line 1 end] \t]
		}
	} else {
package require json
		if {$user ne ""} {
			lappend options -u $value
		}
		set json [exec squeue --json {*}$options]
		set d [json::json2dict $json]
		set list [dict get $d jobs]
		foreach el $list {
			unset -nocomplain a
			set a(tasks) ""
			set a(run) ?
			set a(runversion) ?
			set a(JAT_prio) ?
			array set a $el
			if {$a(submit_time) eq ""} {
				set a(submit_time) ""
			} else {
				set a(submit_time) [clock format $a(submit_time) -format "%Y-%m-%d %H:%M:%S"]
			}
			if {$a(start_time) eq ""} {
				set a(start_time) ""
			} else {
				set a(start_time) [clock format $a(start_time) -format "%Y-%m-%d %H:%M:%S"]
			}
			set resultline $a(job_id),$a(array_task_id)
			if {[regexp {^j([^#]+)\.([0-9_-]+)\#} $a(name) temp a(run) a(runversion)]} {
			} elseif {[regexp {^j([^#]+)\.([0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]_[0-9][0-9-]+)\.\.\.\.} $a(name) temp a(run) a(runversion)]} {
			} elseif {[regexp {^j(.+)\.([0-9_-]+)$} [file root $a(name)] temp a(run) a(runversion)]} {
			} elseif {[regexp {^j([^.]+)\.([0-9_-]+)\.} $a(name) temp a(run) a(runversion)]} {
			} else {
				regexp {^j(.+)\.([0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]_[0-9][0-9]-[0-9][0-9]-[0-9][0-9]-[0-9][0-9][0-9])\.} $a(name) temp a(run) a(runversion)
			}
			foreach field {job_id tasks job_state submit_time start_time priority JAT_prio user_name partition cpus run runversion name standard_error standard_output} {
				lappend resultline [get a($field) .]
			}
			lappend result $resultline
		}
		set result [bsort -index 0 $result]
		puts [join {id tasks state submissiontime starttime priority JAT_prio owner queue slots run runversion name standard_error standard_output} \t]
		foreach line $result {
			puts [join [lrange $line 1 end] \t]
		}
	}
}
