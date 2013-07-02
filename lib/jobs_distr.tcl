proc job_process_distr_running {jobs} {
	global cgjob_running
	foreach job $jobs {
		if {$job eq ""} continue
		if {$job eq "q"}  {
			return 1
		}
		if {[get cgjob_running($job) 0]} {
			return 1
		}
	}
	return 0
}

proc job_process_distr_progress {job args} {
	global cgjob_info
	if {![info exists cgjob_info($job,out)]} {
		set cgjob_info($job,out) [open $job.out w]
	}
	set f $cgjob_info($job,out)
	foreach line $args {
		puts $f $line
		# job_log $job "$line"
	}
	flush $f
}

proc tail {file {num 1}} {
	set f [open $file]
	set result {}
	while {1} {
		set line [gets $f]
		lappend result $line
		if {[llength $result] > $num} {set result [lrange $result end-$num end]}
		if {[eof $f]} break
	}
	close $f
	incr num -1
	return [lrange $result end-$num end]
}

proc job_process_distr_watchdog {} {
	global cgjob cgjob_running cgjob_info
	if {!$cgjob(silent)} {puts stderr "   -=- [llength [array names cgjob_running]] running, [llength $cgjob(queue)] in queue"}
	set pids [exec ps -eo pid]
	foreach job [array names cgjob_running] {
		if {![inlist $pids $cgjob_running($job)]} {
			puts "[job_timestamp] watchdog stopped $job"
			job_process_distr_done $job
		}
	}
	after 1000 job_process_distr_watchdog
}

proc job_process_distr_finishjob {job} {
	global cgjob cgjob_running cgjob_info cgjob_id cgjob_ptargets
	set targets [get cgjob_info($job,targets) ""]
	set ptargets [get cgjob_info($job,ptargets) ""]
	# unset cgjob_id($target) to indicate they are no longer running
	# we no longer write results of the check to the log, as the job has already done that
	foreach target $targets {
		unset -nocomplain cgjob_id($target)
	}
	if {[llength $ptargets]} {
		foreach ptarget $ptargets {
			unset -nocomplain cgjob_ptargets($ptarget)
		}
	}
	if {[info exists cgjob_info($job,out)]} {
		catch {close $cgjob_info($job,out)}
	}
	job_logclose $job
	unset -nocomplain cgjob_running($job)
	unset -nocomplain cgjob_info($job,chan)
}

proc job_process_distr_done {job args} {
	global cgjob cgjob_running cgjob_info
	# job has ended
	set jobname [file tail $job]
	set line [tail $job.log 3]
	if {![file exists $job.finished]} {
		job_log $job "-----> job $jobname failed: did not finish\n"
	} elseif {[regexp failed $line]} {
		file delete $job.finished
		job_log $job "-----> job $job failed: not all targets made\n"
	} else {
		job_log $job "-------------------- end [file tail $job] --------------------"
		job_log $job "-----> job $jobname finished successfully\n"
	}
	if {!$cgjob(silent)} {puts stderr "   -=- [llength [array names cgjob_running]] running, [llength $cgjob(queue)] in queue"}
	job_process_distr_finishjob $job
	# other jobs to do, if not signal end to vwait by setting var cgjob_exit
#	if {![llength [array names cgjob_running]] && ![llength $cgjob(queue)]} {
#	}
	update
	if {[llength $cgjob(queue)]} {
		after idle job_process_distr
	} elseif {![llength [array names cgjob_running]]} {
		set ::cgjob_exit 1
	}
	update
}

proc job_status {} {
	global cgjob_id
	set result ""
	foreach name [array names cgjob_id] {
		lappend result "$name: $cgjob_id($name)"
	}
	return $result
}

proc job_process_distr {} {
	global cgjob cgjob_id cgjob_running cgjob_pid cgjob_ptargets cgjob_info

	update
	set queue $cgjob(queue)
	if {![llength $queue]} return
	set cgjob(queue) {}
	set jobroot [pwd]
	set initlen [llength $queue]
	# join [list_subindex $queue {0 1 2 3 4 5 6}] \n
	while {[llength $queue]} {
		set countrunning 0
		foreach job [array names cgjob_running] {
			incr countrunning [get cgjob_info($job,cores) 1]
		}
		if {$countrunning >= $cgjob(distribute)} {
			break
		}
		set line [list_shift queue]
		foreach {jobid jobname job_logdir pwd deps foreach ftargetvars ftargets fptargets fskip code submitopts} $line break
		cd $pwd
		set job [job_logname $job_logdir $jobname]
		# check foreach deps, skip if not fullfilled
		# check for foreach patterns, expand into one ore more entries in the queue
		if {[llength $foreach]} {
			if {[catch {job_finddeps $job $foreach ftargetvars 1 fids time} fadeps]} {
				if {[regexp {^missing dependency} $fadeps]} {
					job_log $job "$fadeps"
				} elseif {[regexp {^ptargets hit} $fadeps]} {
					# if one of the deps hit a ptarget, we cannot continue:
					# future jobs deps may depend on this lines targets,
					# which depend on the outcome of the ptarget job
					# we wait for the ptarget job to finish, by breaking the loop (will reenter after next job is finished)
					job_logclear $job
					job_log $job "blocking at $jobname: $fadeps"
					lappend cgjob(queue) $line
					job_logclose $job
					break
				} else {
					job_log $job "error in foreach dependencies for $jobname: $fadeps"
				}
				job_log $job "-----> job $jobname skipped: dependencies not found"
				job_process_distr_finishjob $job
				continue
			}
			set temp {}
			# make foreach empty
			lset line 5 {}
			foreach fdep $fadeps ftargetvar $ftargetvars {
				lset line 1 $jobname-$fdep
				lset line 4 [list $fdep {*}$deps]
				lset line 6 $ftargetvar
				lappend temp $line
			}
			set queue [list_concat $temp $queue]
			job_logclose $job
			continue
		}
		set opos [lsearch $submitopts -cores]
		if {$opos != -1} {
			set value [lindex $submitopts [incr opos]]
			set cgjob_info($job,cores) $value
			if {$value <= $cgjob(distribute) && $value > ($cgjob(distribute)-$countrunning)} continue
		}
		job_lognf $job "==================== $jobname ===================="
		cd $pwd
		# check deps, skip if not fullfilled
		if {[catch {job_finddeps $job $deps newtargetvars 0 ids time $ftargetvars} adeps]} {
			# dependencies not found (or error) -> really skip job
			if {[regexp {^missing dependency} $adeps]} {
				job_log $job "$adeps"
			} elseif {[regexp {^ptargets hit} $adeps]} {
				# if one of the deps hit a ptarget, we wil not continue:
				# future jobs deps may depend on this lines targets,
				# which depend on the outcome of the ptarget job
				# we wait for the ptarget job to finish, by breaking the loop (will reenter after next job is finished)
				job_logclear $job
				job_log $job "blocking at $jobname: $adeps"
				lappend cgjob(queue) $line
				job_logclose $job
				break
			} else {
				job_log $job "error in dependencies for $jobname: $adeps"
			}
			job_log $job "-----> job $jobname skipped: dependencies not found"
			job_process_distr_finishjob $job
			continue
		}
		if {[job_process_distr_running $ids]} {
			set depsrunning 1
		} else {
			set depsrunning 0
		}
		set targetvars $ftargetvars
		lappend targetvars {*}$newtargetvars
		# check skip targets, if already done or running, skip
		if {!$cgjob(force) && [llength $fskip]} {
			set skip [job_targetsreplace $fskip $targetvars]
			if {[llength $skip] && [job_checktargets $job $skip $time running]} {
				job_log $job "skipping $jobname: skip targets already completed or running"
				job_process_distr_finishjob $job
				continue
			}
		}
		# check targets, if already done or running, skip
		set targets [job_targetsreplace $ftargets $targetvars]
		set newtargets 0
		if {![job_checktargets $job $targets $time targetsrunning]} {
			set newtargets 1
		}
		set ptargets [job_targetsreplace $fptargets $targetvars]
		if {[llength $ptargets] && ![llength [job_findptargets $ptargets]]} {
			set newtargets 1
		}
		# indicate targets are in the queue, so job_finddeps will find them
		foreach target $targets {
			if {![info exists cgjob_id($target)]} {
				set cgjob_id($target) q
			}
		}
		set cgjob_info($job,targets) $targets
		# if deps or overlapping targets are not finished yet, put line back in the queue and skip running
		if {$depsrunning || [llength $targetsrunning]} {
			job_logclear $job
			lappend cgjob(queue) $line
			job_logclose $job
			continue
		}
		job_log $job
		if {$cgjob(force)} {
			foreach target $targets {
				# job_backup $target 1
				file delete $target
			}
			set newtargets 1
		}
		if {!$newtargets} {
			job_log $job "skipping $jobname: targets already completed or running"
			foreach ptarget $ptargets {
				unset -nocomplain cgjob_ptargets($ptarget)
			}
			job_process_distr_finishjob $job
			continue
		}
		# put ptargets in the queue only if we will actually be running the job
		# Otherwise they would keep dependend jobs from being started
		# (normal targets do not cause this, so these can stay)
		foreach ptarget $ptargets {
			if {![info exists cgjob_id($ptarget)]} {
				set cgjob_ptargets($ptarget) q
			}
		}
		set cgjob_info($job,ptargets) $ptargets
		job_log $job "-------------------- running $jobname --------------------"
		# run code
		set cmd {#!/bin/sh}
		append cmd "\n\# the next line restarts using cgsh \\\n"
		append cmd {exec cg source "$0" "$@"} \n
		append cmd [job_generate_code $job $pwd $adeps $targetvars $targets $ptargets $code]\n
		catch {file delete $job.finished}
		set runfile $job.run
		file_write $runfile $cmd
		file attributes $runfile -permissions u+x
		Extral::bgexec -no_error_redir -pidvar cgjob_pid  \
			-command [list job_process_distr_done $job] \
			-progresscommand [list job_process_distr_progress $job] \
			$runfile 2> $job.err
		foreach target $targets {
			set cgjob_id($target) $job
		}
		foreach ptarget $ptargets {
			set cgjob_ptargets($ptarget) $job
		}
		set cgjob_running($job) $cgjob_pid
		lappend running $job
	}

	if {[llength $queue]} {
		lappend cgjob(queue) {*}$queue
	}
	set running [array names cgjob_running]
	cd $jobroot
	if {![llength $running]} {
		if {[llength $cgjob(queue)]} {
			puts stderr "[llength $cgjob(queue)] jobs in queue cannot be executed: [list_subindex $cgjob(queue) 1]"
		}
		set ::cgjob_exit 1
	}
	after idle update
}

proc job_process_distr_wait {} {
	global cgjob cgjob_exit
	update
	unset -nocomplain cgjob_exit
	after cancel job_process_distr_watchdog
	set running [array names cgjob_running]
	if {![llength $running]} {
		if {[llength $cgjob(queue)]} {
			after idle job_process_distr
		} else {
			return
		}
	}
	after 1000 job_process_distr_watchdog
	if {[catch {vwait cgjob_exit} e]} {
		puts "job_wait warning: $e"
	}
	puts "All jobs done"
	unset -nocomplain cgjob_exit
}

