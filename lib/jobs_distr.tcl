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
	foreach line $args {
		job_log $job "$line"
	}
}

proc job_process_distr_done {job targets ptargets args} {
	global cgjob cgjob_id cgjob_running cgjob_ptargets
	job_log $job "-------------------- end [file tail $job] --------------------"
	# job has ended
	unset -nocomplain cgjob_running($job)
	# unset cgjob_id($target) to indicate they are nolonger running
	# we no longer write results of the check to the log, as the job has already done that
	foreach target $targets {
		unset -nocomplain cgjob_id($target)
	}
	if {[llength $ptargets]} {
		foreach ptarget $ptargets {
			unset cgjob_ptargets($ptarget)
		}
	}
	job_logclose $job
	# other jobs to do, if not signal end to vwait by setting var cgjob_exit
	update
	if {[llength $cgjob(queue)]} {
		after idle job_process_distr
	} elseif {![llength [array names cgjob_running]]} {
		set ::cgjob_exit 1
	}
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
	global cgjob cgjob_id cgjob_running cgjob_pid cgjob_ptargets

	set queue $cgjob(queue)
	update
	if {![llength $queue]} return
	set cgjob(queue) {}
	set jobroot [pwd]
	set running [array names cgjob_running]
	set initlen [llength $queue]
	# join [list_subindex $queue {0 1 2 3 4 5 6}] \n
	while {[llength $queue]} \
 {
		if {[llength $running] > $cgjob(distribute)} {
			break
		}
		set line [list_shift queue]
		foreach {jobid jobname pwd deps foreach ftargetvars ftargets fptargets fskip code submitopts} $line break
puts $jobname

		cd $pwd
		set job [job_logname $jobname]
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
					break
				} else {
					job_log $job "error in foreach dependencies for $jobname: $fadeps"
				}
				job_log $job "job $jobname skipped: dependencies not found"
				continue
			}
			set temp {}
			# make foreach empty
			lset line 4 {}
			foreach fdep $fadeps ftargetvar $ftargetvars {
				lset line 1 $jobname-$fdep
				lset line 3 [list $fdep {*}$deps]
				lset line 5 $ftargetvar
				lappend temp $line
			}
			set queue [list_concat $temp $queue]
			continue
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
				break
			} else {
				job_log $job "error in dependencies for $jobname: $adeps"
			}
			job_log $job "job $jobname skipped: dependencies not found"
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
		# if deps or overlapping targets are not finished yet, put line back in the queue and skip running
		if {$depsrunning || [llength $targetsrunning]} {
			job_logclear $job
			lappend cgjob(queue) $line
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
		job_log $job "-------------------- running $jobname --------------------"
		# run code
		set cmd {#!/bin/sh}
		append cmd "\n\# the next line restarts using cgsh \\\n"
		append cmd {exec cg source "$0" "$@"} \n
		append cmd [job_generate_code $job $pwd $adeps $targetvars $targets $ptargets $code]\n
		set runfile $job.run
		file_write $runfile $cmd
		file attributes $runfile -permissions u+x
		Extral::bgexec -no_error_redir -pidvar cgjob_pid \
			-command [list job_process_distr_done $job $targets $ptargets] \
			-progresscommand [list job_process_distr_progress $job] \
			$runfile 2> $runfile.stderr
		foreach target $targets {
			set cgjob_id($target) $job
		}
		foreach ptarget $ptargets {
			set cgjob_ptargets($ptarget) $job
		}
		set cgjob_running($job) $cgjob_pid
		lappend running $job
	}

	set running [array names cgjob_running]
	lappend cgjob(queue) {*}$queue
	cd $jobroot
	after idle update
}

proc job_process_distr_wait {} {
	global cgjob cgjob_exit
	unset -nocomplain cgjob_exit
	set running [array names cgjob_running]
	if {![llength cgjob(queue)] && ![llength $running]} return
	if {[catch {vwait cgjob_exit} e]} {
		puts "job_wait warning: $e"
	}
	unset -nocomplain cgjob_exit
}

