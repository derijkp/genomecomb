proc job_process_distr_init {} {
	global cgjob cgjob_distr cgjob_distr_running cgjob_pid cgjob_exit
	unset -nocomplain cgjob_distr
	unset -nocomplain cgjob_distr_running
	set cgjob_distr(queue) {}
	set cgjob_distr(num) 0
	# we will use par (parallel) code with some specifics for distr
	if {[info commands job_process_par] eq ""} {auto_load job_process_par}
	interp alias {} job_process {} job_process_par
	interp alias {} job_running {} job_running_distr
	interp alias {} job_wait {} job_process_distr_wait
	interp alias {} job_process_par_submit {} job_process_distr_submit
}

proc job_running_distr {jobnum} {
	global cgjob_distr_running
	info exists cgjob_distr_running($jobnum)
}

proc job_process_distr_jobmanager {} {
	global cgjob cgjob_distr cgjob_distr_running cgjob_exit cgjob_distr_queue
	after cancel job_process_distr_jobmanager
	set pids [exec ps -eo pid]
	foreach jobnum [ssort -natural [array names cgjob_distr_running]] {
		if {![inlist $pids [lindex $cgjob_distr_running($jobnum) 0]]} {
			if {!$cgjob(silent)} {puts "   -=- ending [lindex $cgjob_distr_running($jobnum) 1] ($jobnum)"}
			unset cgjob_distr_running($jobnum)
			unset cgjob_distr_queue($jobnum)
		}
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
		if {!$cgjob(silent)} {puts "   -=- starting $name"}
		set cgjob_pid [lindex [exec $runfile > $job.out 2> $job.err &] end]
		set cgjob_distr_running($jobnum) [list $cgjob_pid $name]
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

proc job_process_distr_submit {job runfile args} {
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
	if {![info exists cgjob_distr(num)]} {set cgjob_distr(num) 0}
	set jobnum [incr cgjob_distr(num)]
	catch {file delete $job.out}
	catch {file delete $job.err}
	lappend cgjob_distr(queue) [list $jobnum $deps $name $job $runfile $options]
	set cgjob_distr_queue($jobnum) 1
	job_process_distr_jobmanager
	return $jobnum
}

proc job_process_distr_wait {} {
	global cgjob cgjob_exit cgjob_running
	update
	unset -nocomplain cgjob_exit
	after cancel job_process_distr_jobmanager
	after 1000 job_process_distr_jobmanager
	if {[catch {vwait cgjob_exit} e]} {
		puts "job_wait warning: $e"
	}
	# puts "All jobs done"
	unset -nocomplain cgjob_exit
}

