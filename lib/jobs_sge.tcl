proc job_process_sge_init {} {
	# we will use par (parallel) code with some specifics for sge
	if {[info commands job_process_par] eq ""} {auto_load job_process_par}
	set ::job_method_info {}
	interp alias {} job_process {} job_process_par
	interp alias {} job_running {} job_running_sge
	interp alias {} job_wait {} job_process_sge_wait
	interp alias {} job_process_par_submit {} job_process_sge_submit
}

proc job_running_sge {jobid} {
	expr {[isint $jobid] && ![catch {exec qstat -j $jobid}]} 
}

proc job_status_sge {job {jobloginfo {}}} {
	global cgjob_distr_running
	if {$jobloginfo eq ""} {
		if {![file exists $job.log]} {return unkown}
		set jobloginfo [job_parse_log $job $totalduration]
	}
	foreach {failed starttime endtime run duration totalduration} $jobloginfo break
	if {$failed} {
		return error
	} elseif {$endtime ne ""} {
		return finished
	} else {
		set jobnum [job_process_par_jobid $job]
		if {[isint $jobnum]} {
			if {$starttime ne ""} {
				return running
			} else {
				return submitted
			}
		} else {
			return error
		}
	}
}

proc job_process_sge_submit {job runfile args} {
	global cgjob
	set options {}
	set soft {}
	set hard {}
	set priority [get cgjob(priority) 0]
	set pos 0
	foreach {opt value} $args {
		switch -- $opt {
			-deps {
				set value [list_remove $value {}]
				if {[llength $value]} {
					lappend options -hold_jid [join $value ,]
				}
				incr pos 2
			}
			-priority {
				set priority $value
			}
			-cores {
				if {![info exists cgjob_distr(no_local_pe)]} {
					set cgjob_distr(no_local_pe) [catch {exec qconf -sp local}]
				}
				if {!$cgjob_distr(no_local_pe)} {
					lappend hard -pe local $value
				}
				incr pos 2
			}
			-hard {
				lappend hard -l $value
				incr pos 2
			}
			-soft {
				lappend soft -l $value
				incr pos 2
			}
			-host {
				lappend hard -l hostname=$value
				incr pos 2
			}
			-io {
				lappend soft -l io=$value
				incr pos 2
			}
			-mem {
				lappend hard -l mem_free=$value
				incr pos 2
			}
			-- {
				incr pos 1
				break
			}
			default {
				break
			}
		}
	}
	if {[llength $soft]} {
		lappend options -soft {*}$soft
	}
	if {[llength $hard]} {
		lappend options -hard {*}$hard
	}
	set name "[file tail $job] $job"
	regsub -all {[^A-Za-z0-9_.-]} $name __ name
	if {[string length $name] > 200} {
		set name [string range $name 0 100]....[string range $name end-100 end]
	}
	set dir [file dir $job]
	catch {file delete $job.finished}
	catch {file delete $job.out}
	catch {file delete $job.err}
	set jnum [exec qsub -N j$name -q all.q -o $job.out -e $job.err -p $priority {*}$options $runfile]
	regexp {[0-9]+} $jnum jobnum
	return $jobnum
}

if 0 {
# debug proc
	set curjobnum 1
	proc job_process_sge_submit {job runfile args} {
		global curjobnum
		return [incr curjobnum]
	}
	interp alias {} job_process_par_submit {} job_process_sge_submit
}

proc job_process_sge_wait {} {
	global cgjob
	job_logfile_par_close
}
