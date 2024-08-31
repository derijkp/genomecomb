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
	interp alias {} job_runall {} job_runall_par
	interp alias {} job_running {} job_running_distr
	interp alias {} job_wait {} job_wait_distr
	interp alias {} job_process_submit_par {} job_process_submit_distr
}

proc job_running_distr {jobnum} {
	global cgjob cgjob_distr_running cgjob_distr_queue
	if {![info exists cgjob_distr_running($jobnum)]} {return 0}
	set job [lindex $cgjob_distr_running($jobnum) 1]
	if {$cgjob(dsubmit) eq ""} {
		if {![catch {file_read [job.file pid $job]} pid]} {
			if {![catch {exec kill -0 $pid} msg]} {
				return 1
			}
			file delete [job.file pid $job]
		}
	} elseif {$cgjob(dsubmit) eq "sge"} {
		set jnum [lindex $cgjob_distr_running($jobnum) 0]
		if {[job_running_sge $jnum]} {
			return 1
		}
	} elseif {$cgjob(dsubmit) eq "slurm"} {
		set jnum [lindex $cgjob_distr_running($jobnum) 0]
		if {[job_running_slurm $jnum]} {
			return 1
		}
	} else {
		set jnum [lindex $cgjob_distr_running($jobnum) 0]
		if {[exec $cgjob(dsubmit) running $jnum]} {
			return 1
		}
	}
	if {!$cgjob(silent)} {puts "   -=- [timestamp] ending $job ($jobnum)"}
	unset -nocomplain cgjob_distr_running($jobnum)
	unset -nocomplain cgjob_distr_queue($jobnum)
	return 0
}

proc job_status_distr {job {jobloginfo {}}} {
	global cgjob cgjob_distr_running
	set totalduration {0 0}
	if {$jobloginfo eq ""} {
		if {![file exists [job.file log $job]]} {return unkown}
		set jobloginfo [job_parse_log $job]
	}
	foreach {status starttime endtime run duration} $jobloginfo break
	if {$status ni {submitted running}} {return $status}
	if {![info exists cgjob(pid)] || [catch {exec kill -0 $cgjob(pid)}]} {return error}
	if {$status eq "submitted"} {return $status}
	if {![catch {file_read [job.file pid $job]} pid] && ![catch {exec kill -0 $pid}]} {
		return running
	} else {
		return error
	}
}


proc job_process_distr_progress {running} {
	global cgjob
	if {[get cgjob(distr_count) 0] >= 79} {
		global cgjob_distr_running
		puts .
		foreach job $running {
			puts "   -=- Running $cgjob_distr_running($job)"
		}
		puts "   -=- in queue: [llength $::cgjob_distr(queue)] jobs"
		set cgjob(distr_count) 0
	} else {
		incr cgjob(distr_count)
		puts -nonewline stderr .
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
		if {!$cgjob(silent)} {
			job_process_distr_progress $running
		}
		after 1000 job_process_distr_jobmanager
		return
	}
	set torun [expr {$maxrunning - $countrunning}]
	set pos -1
	set added {}
	if {$cgjob(dmaxmem) ne ""} {
		set currentmem 0
		foreach jobnum [bsort [array names cgjob_distr_running]] {
			set currentmem [expr {$currentmem + [lindex $cgjob_distr_running($jobnum) 2]}]
		}
	}
	foreach line $cgjob_distr(queue) {
		incr pos
		foreach {jobnum deps name job runfile options mem} $line break
		# foreach {jobid jobname job_logdir pwd deps ftargetvars ftargets fskip checkcompressed code submitopts frmtargets precode jobforce optional cores} $line break
		set do 1
		foreach dep $deps {
			if {[info exists cgjob_distr_queue($dep)]} {set do 0 ; break}
		}
		if {!$do} continue
		if {[llength [list_common $deps $running]]} continue
		if {$cgjob(dsubmit) eq ""} {
			set waitingformem 0
			set waitingmaxmem 0
			if {$cgjob(dmaxmem) ne ""} {
				set testmem [expr {$currentmem + $mem}]
				if {$testmem > $cgjob(dmaxmem)} {
					# if {!$cgjob(silent)} {puts "   -=- not enough memory free to start $job: needs [display_memory $mem]"}
					incr waitingformem
					if {$mem > $waitingmaxmem} {set waitingmaxmem $mem}
					continue
				}
				set currentmem $testmem
			}
			if {!$cgjob(silent)} {puts "   -=- [timestamp] starting $job"}
			set cgjob_pid [lindex [exec $runfile > [job.file out $job] 2> [job.file err $job] &] end]
			file_write [job.file pid $job] $cgjob_pid
			set cgjob_distr_running($jobnum) [list $cgjob_pid $job $mem]
		} elseif {$cgjob(dsubmit) eq "sge"} {
			set jnum [job_process_submit_sge $job $cmd -mem $mem]
			set cgjob_distr_running($jobnum) [list $jnum $job $mem]
		} elseif {$cgjob(dsubmit) eq "slurm"} {
			set jnum [job_process_submit_slurm $job $cmd -mem $mem]
			set cgjob_distr_running($jobnum) [list $jnum $job $mem]
		} else {
			set jnum [exec $cgjob(dsubmit) submit -mem $mem -o [job.file out $job] -e [job.file err $job] $runfile]
			set cgjob_distr_running($jobnum) [list $jnum $job $mem]
		}
		incr torun -1
		lappend added $pos
		if {$torun == 0} break
	}
	set running [array names cgjob_distr_running]
	set countrunning [llength $running]
	if {[llength $added]} {
		set cgjob_distr(queue) [list_sub $cgjob_distr(queue) -exclude $added]
		set countqueue [llength cgjob_distr(queue)]
		if {$waitingformem} {set temp " ($waitingformem waiting for memory, max needed [display_memory $waitingmaxmem])"} else {set temp {}}
		if {!$cgjob(silent)} {puts "   -=- [llength [array names cgjob_distr_running]] running $countqueue in queue$temp"}
		after 200 job_process_distr_jobmanager
	} else {
		set countqueue [llength $cgjob_distr(queue)]
		if {!$countqueue && !$countrunning} {
			set cgjob_exit 1
			return
		}
		if {!$cgjob(silent)} {
			job_process_distr_progress $running
		}
		after 1000 job_process_distr_jobmanager
	}
}

proc job_process_submit_distr {job cmd args} {
#set jobnum [incr cgjob_distr(num)]
#return $jobnum
	global cgjob cgjob_distr cgjob_distr_queue
	set options {}
	set deps {}
	set cores 1
	set io 1
	set pos 0
	set mem 500m
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
				set mem $value
				incr pos 2
			}
			-time {
				# not used yet
				set time $value
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
	catch {file delete [job.file finished $job]}
	catch {file delete [job.file ok $job]}
	if {![info exists cgjob_distr(num)]} {set cgjob_distr(num) 0}
	set jobnum [incr cgjob_distr(num)]
	catch {file delete [job.file out $job]}
	catch {file delete [job.file err $job]}
	# make runscript
	set runcmd {}
	append runcmd {#!/bin/sh}
	append runcmd \n
	append runcmd {#$ -S /bin/bash} \n
	append runcmd {#$ -V} \n
	append runcmd {#$ -cwd} \n
	append runcmd $cmd
	set runfile [job.file run $job]
	file_write $runfile $runcmd
	file attributes $runfile -permissions u+x
	# submit
	if {$cgjob(nosubmit) || $cgjob(dry)} {
		putslog "nosubmit run, would be distr_submit: added to queue [list jobnum $jobnum deps $deps name $name job $job runfile $runfile options $options]"
	} else {
		putslog "distr_submit: added to queue [list jobnum $jobnum deps $deps name $name job $job runfile $runfile options $options mem $mem]"
		if {$mem eq ""} {set mem 500m}
		set mem [entry_memory $mem]
		lappend cgjob_distr(queue) [list $jobnum $deps $name $job $runfile $options $mem]
		set cgjob_distr_queue($jobnum) 1
		job_process_distr_jobmanager
	}
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
	if {$statusok} {
		putslog "all jobs finished"
	} else {
		set nrerrors [lindex [cg select -g all -q {$status eq "error"} $cgjob(logfile).error] end]
		puts stderr "\nAnalysis ended with $nrerrors errors"
		puts stderr "You can use the following command for an overview of failed jobs and their error messages:"
		puts stderr "cg error_report $cgjob(logfile).error"
	}
}

proc job_wait_distr {} {
	global cgjob cgjob_exit
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
