proc job_process_sge_runninginit {} {
	global cgjob_running_sge
	set cgjob_running_sge() 1
	catch {exec qstat -xml} jobxml
	set jobs [regexp -all -inline {<job_list.+?</job_list>} $jobxml]
	unset -nocomplain cgjob_running_sge
	foreach job $jobs {
		set task {} ; set name {}
		regexp {<tasks>(.*?)</tasks>} $job temp task
		regexp {<JB_name>(.*?)</JB_name>} $job temp name
		set name [string range $name 1 end]
		if {[string is int $task]} {
			set cgjob_running_sge([list $name $task]) 1
		} else {
			set tasks {}
			regexp {([0-9]+)-([0-9]+)} $task temp start end
			for {} {$start <= $end} {incr start} {
				set cgjob_running_sge([list $name $start]) 1
			}
		}
	}
}

proc job_process_sge_running {jobs} {
	global cgjob_running cgjob_running_sge
	if {![info exists cgjob_running_sge()]} job_process_sge_runninginit
	foreach job $jobs {
		if {$job eq ""} continue
		if {$job eq "q"}  {
			return 1
		}
		if {[info exists cgjob_running($job)]} {
			set jobid $cgjob_running($job)
			if {[catch {exec qstat -j $jobid}]} {
				set cgjob_running($job) ""
			} else {
				return 1
			}
		}
	}
	return 0
}

proc job_process_sge_jobid {job} {
	global cgjob_running cgjob_running_sge
	if {![file exists $job.jobid]} {
		return ""
	}
	set jobid [file_read $job.jobid]
	if {[catch {exec qstat -j $jobid}]} {
		return ""
	} else {
		return $jobid
	}
}

proc job_process_sge_marktargets {targets ptargets id} {
	global cgjob_id
	foreach target $targets {
		if {![info exists cgjob_id($target)]} {
			set cgjob_id($target) $id
		}
	}
	foreach ptarget $ptargets {
		if {![info exists cgjob_id($ptarget)]} {
			set cgjob_ptargets($ptarget) $id
		}
	}
}

proc job_process_sge {} {
	global cgjob cgjob_id cgjob_running cgjob_pid cgjob_ptargets
	set queue $cgjob(queue)
	update
	if {![llength $queue]} return
	set cgjob(queue) {}
	set jobroot [pwd]
	while {[llength $queue]} {
		set line [list_shift queue]
		foreach {jobid jobname pwd deps foreach ftargetvars ftargets fptargets fskip code submitopts} $line break
		cd $pwd
		set job [job_logname $jobname]
		file mkdir [file dir $job]
		# check foreach deps, skip if not fullfilled
		# check for foreach patterns, expand into one ore more entries in the queue
		if {[llength $foreach]} {
			if {[catch {job_finddeps $job $foreach ftargetvars fids} fadeps]} {
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
				job_log $job "job $jobname failed"
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
		# check if job is already running, if so, mark targets with jobid
		set jobnum [job_process_sge_jobid $job]
		if {[isint $jobnum]} {
			job_process_sge_marktargets [file_read $job.targets] [file_read $job.$ptargets] $jobnum
		}
		# check deps, skip if not fullfilled
		if {[catch {job_finddeps $job $deps newtargetvars ids $ftargetvars} adeps]} {
			# dependencies not found (or error) -> really skip job
			if {[regexp {^missing dependency} $adeps]} {
				job_log $job "$adeps"
			} elseif {[regexp {^ptargets hit} $adeps]} {
				error "ptargets not supported by sge (yet)"
				# if one of the deps hit a ptarget, we cannot continue:
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
			job_log $job "job $jobname failed"
			continue
		}
		file_write $job.deps $adeps
#		if {[job_process_sge_running $ids]} {
#			set depsrunning 1
#		} else {
#			set depsrunning 0
#		}
		set targetvars $ftargetvars
		lappend targetvars {*}$newtargetvars
		# check skip targets, if already done or running, skip
		if {!$cgjob(force) && [llength $fskip]} {
			set skip [job_targetsreplace $fskip $targetvars]
			if {[llength $skip] && [job_checktargets $job $skip running]} {
				job_log $job "skipping $jobname: skip targets already completed or running"
				continue
			}
		}
		# check targets, if already done or running, skip
		set targets [job_targetsreplace $ftargets $targetvars]
		file_write $job.targets $targets
		set newtargets 0
		if {![job_checktargets $job $targets targetsrunning]} {
			set newtargets 1
		}
		set ptargets [job_targetsreplace $fptargets $targetvars]
		file_write $job.ptargets $ptargets
		if {[llength $ptargets] && ![llength [job_findptargets $ptargets]]} {
			set newtargets 1
		}
		# indicate targets are in the queue, so job_finddeps will find them
		job_process_sge_marktargets $targets $ptargets q
#		# if deps or overlapping targets are not finished yet, put line back in the queue and skip running
#		if {$depsrunning || [llength $targetsrunning]} {
#			job_logclear $job
#			lappend cgjob(queue) $line
#			continue
#		}
		job_log $job
		if {$cgjob(force)} {
			foreach target [list_concat $targets $ptargets] {
				job_backup $target 1
			}
			set newtargets 1
		}
		if {!$newtargets} {
			job_log $job "skipping $jobname: targets already completed or running"
			continue
		}
		job_log $job "-------------------- submitting $jobname --------------------"
		# run code
		set cmd {#!/bin/sh}
		append cmd \n
		append cmd {#$ -S /bin/bash} \n
		append cmd {#$ -V} \n
		append cmd {#$ -cwd} \n
		append cmd "\n\# the next line restarts using cgsh \\\n"
		append cmd {exec cg source "$0" "$@"} \n
		append cmd "file_add \{$job.log\} \"[job_timestamp]\\tstarting $jobname\"\n"
		append cmd [job_generate_code $job $pwd $adeps $targetvars $targets $code]\n
		append cmd "file_add \{$job.log\} \"[timestamp]\\tending $jobname\"\n"
		set runfile $job.run
		file_write $runfile $cmd
		file attributes $runfile -permissions u+x
		# submit job
		set jobnum [job_process_sge_submit $job $runfile -deps $ids]
		file_write $job.jobid $jobnum
		job_process_sge_marktargets $targets $ptargets $jobnum
		set cgjob_running($job) $jobnum
	}
	lappend cgjob(queue) {*}$queue
	cd $jobroot
	after idle update
}

proc job_process_sge_submit {job runfile args} {
	set options {}
	set soft {}
	set hard {}
	set pos 0
	foreach {opt value} $args {
		switch -- $opt {
			-deps {
				if {[llength $value]} {
					lappend options -hold_jid $value
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
	regsub -all / $name __ name
	if {[string length $name] > 200} {
		set name [string range $name 0 100]....[string range $name end-100 end]
	}
	set dir [file dir $job]
	set jnum [exec qsub -N j$name -q all.q -o $job.out -e $job.err {*}$options $runfile]
	regexp {[0-9]+} $jnum jobnum
	return $jobnum
}

proc job_process_sge_wait {} {
}

