proc cg_qjobs {args} {
	set basedir [file_absolute [pwd]]
	set options {}
	set user {}
	set summary 0
	cg_options qjobs args {
		-u {set user $value}
		-s - -summary {
			set summary [true $value]
		}
	} {} 0 0 {
		returns running and waiting jobs on a grid engine cluster in tsv format (so they can by analysed using cg select)
	}
	set result {}
	set type [job_distribute]
	if {$type ni "sge slurm"} {
		if {![catch {
			exec which qstat
		}]} {
			set type sge
		} elseif {![catch {
			exec which squeue
		}]} {
			set type slurm
		} else {
			error "could not determine job manager (supported: sge, slurm)"
		}
	}
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
		if {!$summary} {
			puts [join {id tasks state submissiontime starttime priority JAT_prio owner queue slots run runversion name} \t]
			foreach line $result {
				puts [join [lrange $line 1 end] \t]
			}
		} else {
			unset -nocomplain a
			foreach line $result {
				foreach {run state} [list_sub [lrange $line 1 end] {10 2}] break
				incr a([list $run $state])
			}
		}
	} else {
		package require json
		if {$user ne ""} {
			lappend options -u $value
		}
		set json [exec squeue --json {*}$options]
		set d [json::json2dict $json]
		set result {}
		set header {id tasks state submissiontime starttime priority JAT_prio owner queue slots run runversion name standard_error standard_output	}
		set list [dict get $d jobs]
		foreach el $list {
			unset -nocomplain a
			set a(tasks) ""
			set a(run) ?
			set a(runversion) ?
			set a(JAT_prio) ?
			array set a $el
			if {$a(submit_time) eq ""} {
				set a(submissiontime) ""
			} else {
				set a(submissiontime) [clock format $a(submit_time) -format "%Y-%m-%d %H:%M:%S"]
			}
			if {$a(start_time) eq ""} {
				set a(starttime) ""
			} else {
				set a(starttime) [clock format $a(start_time) -format "%Y-%m-%d %H:%M:%S"]
			}
			set a(id) $a(job_id),[lindex $a(array_task_id) end]
			set a(tasks) [lindex $a(tasks) end]
			set a(state) $a(job_state)
			set a(owner) $a(user_name)
			set a(queue) $a(partition)
			set a(slots) [lindex $a(cpus) end]
			set a(priority) [lindex $a(priority) end]
			if {[regexp {^j([^#]+)\.([0-9_-]+)\#} $a(name) temp a(run) a(runversion)]} {
			} elseif {[regexp {^j([^#]+)\.([0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]_[0-9][0-9-]+)\.\.\.\.} $a(name) temp a(run) a(runversion)]} {
			} elseif {[regexp {^j(.+)\.([0-9_-]+)$} [file root $a(name)] temp a(run) a(runversion)]} {
			} elseif {[regexp {^j([^.]+)\.([0-9_-]+)\.} $a(name) temp a(run) a(runversion)]} {
			} else {
				regexp {^j(.+)\.([0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]_[0-9][0-9]-[0-9][0-9]-[0-9][0-9]-[0-9][0-9][0-9])\.} $a(name) temp a(run) a(runversion)
			}
			set resultline {}
			foreach field $header {
				lappend resultline [get a($field) .]
			}
			lappend result $resultline
		}
		set result [bsort -index 0 $result]
		if {!$summary} {
			puts [join $header \t]
			foreach line $result {
				puts [join $line \t]
			}
		} else {
			set poss [list_cor $header {run state}]
			unset -nocomplain a
			foreach line $result {
				foreach {run state} [list_sub $line $poss] break
				incr a([list $run $state])
			}
		}
	}
	if {$summary} {
		puts [join {run state count} \t]
		foreach key [bsort [array names a]] {
			puts [join [list {*}$key $a($key)] \t]
		}
	}
}
