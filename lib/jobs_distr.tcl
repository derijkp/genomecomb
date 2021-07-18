proc job_process_init_distr {} {
	global cgjob cgjob_distr cgjob_distr_running cgjob_pid cgjob_exit
	unset -nocomplain cgjob_distr
	unset -nocomplain cgjob_distr_running
	set cgjob_distr(queue) {}
	set cgjob_distr(num) 0
	set cgjob(pid) [pid]
	# we will use par (parallel) code with some specifics for distr
	if {[info commands job_process_par] eq ""} {auto_load job_process_par}
	set ::job_method_info [list pid [pid]]
	interp alias {} job_process {} job_process_par
	interp alias {} job_running {} job_running_distr
	interp alias {} job_wait {} job_wait_distr
	interp alias {} job_process_submit_par {} job_process_submit_distr
}

proc job_running_distr {jobnum} {
	global cgjob cgjob_distr_running cgjob_distr_queue
	if {![info exists cgjob_distr_running($jobnum)]} {return 0}
	set job [lindex $cgjob_distr_running($jobnum) 1]
	if {![catch {file_read $job.pid} pid]} {
		if {![catch {exec ps $pid}]} {
			return 1
		}
		if {!$cgjob(silent)} {puts "   -=- ending $job ($jobnum)"}
		file delete $job.pid
	}
	unset -nocomplain cgjob_distr_running($jobnum)
	unset -nocomplain cgjob_distr_queue($jobnum)
	return 0
}

proc job_status_distr {job {jobloginfo {}}} {
	global cgjob cgjob_distr_running
	set totalduration {0 0}
	if {$jobloginfo eq ""} {
		if {![file exists $job.log]} {return unkown}
		set jobloginfo [job_parse_log $job]
	}
	foreach {status starttime endtime run duration} $jobloginfo break
	if {$status ni {submitted running}} {return $status}
	if {![info exists cgjob(pid)] || [catch {exec ps $cgjob(pid)}]} {return error}
	if {$status eq "submitted"} {return $status}
	if {![catch {file_read $job.pid} pid] && ![catch {exec ps $pid}]} {
		return running
	} else {
		return error
	}
}

proc job_process_distr_jobmanager {} {
	global cgjob cgjob_distr cgjob_distr_running cgjob_exit cgjob_distr_queue
	after cancel job_process_distr_jobmanager
	foreach jobnum [bsort [array names cgjob_distr_running]] {
		# job_running_distr checks whether the job is running
		# if it is not, it will (a.o.) clear the entry in cgjob_distr_running
		job_running_distr $jobnum
	}
	set running [array names cgjob_distr_running]
	set countrunning [llength $running]
	set maxrunning [get cgjob(distribute) 4]
	if {$countrunning >= $maxrunning} {
		if {!$cgjob(silent)} {puts -nonewline stderr .}
		after 1000 job_process_distr_jobmanager
		return
	}
	set torun [expr {$maxrunning - $countrunning}]
	set pos -1
	set added {}
	foreach line $cgjob_distr(queue) {
		incr pos
		foreach {jobnum deps name job runfile options} $line break
		set do 1
		foreach dep $deps {
			if {[info exists cgjob_distr_queue($dep)]} {set do 0 ; break}
		}
		if {!$do} continue
		if {[llength [list_common $deps $running]]} continue
		if {!$cgjob(silent)} {puts "   -=- starting $job"}
		set cgjob_pid [lindex [exec $runfile > $job.out 2> $job.err &] end]
		file_write $job.pid $cgjob_pid
		set cgjob_distr_running($jobnum) [list $cgjob_pid $job]
		incr torun -1
		lappend added $pos
		if {$torun == 0} break
	}
	set running [array names cgjob_distr_running]
	set countrunning [llength $running]
	if {[llength $added]} {
		set cgjob_distr(queue) [list_sub $cgjob_distr(queue) -exclude $added]
		set countqueue [llength cgjob_distr(queue)]
		if {!$cgjob(silent)} {puts "   -=- [llength [array names cgjob_distr_running]] running $countqueue in queue"}
		after 200 job_process_distr_jobmanager
	} else {
		set countqueue [llength $cgjob_distr(queue)]
		if {!$countqueue && !$countrunning} {
			set cgjob_exit 1
			return
		}
		if {!$cgjob(silent)} {puts -nonewline stderr .}
		after 1000 job_process_distr_jobmanager
	}
}

proc job_process_submit_distr {job runfile args} {
#set jobnum [incr cgjob_distr(num)]
#return $jobnum
	global cgjob_distr cgjob_distr_queue
	set options {}
	set deps {}
	set cores 1
	set io 1
	set pos 0
	foreach {opt value} $args {
		switch -- $opt {
			-deps {
				set deps [list_remove $value {}]
				incr pos 2
			}
			-cores {
				# not used yet
				set cores $value
				incr pos 2
			}
			-hard {
				incr pos 2
			}
			-soft {
				incr pos 2
			}
			-host {
				incr pos 2
			}
			-io {
				# not used yet
				set io $value
				incr pos 2
			}
			-mem {
				# not used yet
				set mem $value
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
	set name "[file tail $job] $job"
	# replace all invalid chars and replace invisible chars 
	regsub -all {[^A-Za-z0-9_.-]} $name __ name
	set dir [file dir $job]
	catch {file delete $job.finished}
	catch {file delete $job.ok}
	if {![info exists cgjob_distr(num)]} {set cgjob_distr(num) 0}
	set jobnum [incr cgjob_distr(num)]
	catch {file delete $job.out}
	catch {file delete $job.err}
	lappend cgjob_distr(queue) [list $jobnum $deps $name $job $runfile $options]
	set cgjob_distr_queue($jobnum) 1
	job_process_distr_jobmanager
	return $jobnum
}

proc job_logfile_distr_close {} {
	global cgjob
	if {[file exists $cgjob(logfile).running]} {
		job_update $cgjob(logfile).running $cgjob(cleanup) 1 $cgjob(removeold) 1
	}
	set statusok [file exists $cgjob(logfile).finished]
	if {$cgjob(cleanup) eq "allways" || ($cgjob(cleanup) eq "success" && $statusok)} {
		job_cleanup
		set result [glob $cgjob(logfile).finished $cgjob(logfile).running $cgjob(logfile).error]
		job_cleanlogs $result
		# only keep result logfile if -d option was given explicitely
		if {!$cgjob(hasargs)} {file delete $result}
	}
}

proc job_wait_distr {} {
	global cgjob cgjob_exit cgjob_running
	update
	job_logfile_par_close
	unset -nocomplain cgjob_exit
	after cancel job_process_distr_jobmanager
	after 1000 job_process_distr_jobmanager
	if {[catch {vwait cgjob_exit} e]} {
		puts "job_wait warning: $e"
	}
	# puts "All jobs done"
	after cancel job_process_distr_jobmanager
	unset -nocomplain cgjob_exit
	update
	job_logfile_distr_close
}
