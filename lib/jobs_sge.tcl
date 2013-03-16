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
		file delete $job.jobid
		return ""
	} else {
		return $jobid
	}
}

proc job_process_sge_marktargets {targets ptargets id} {
	global cgjob_id cgjob_ptargets
	foreach target $targets {
		if {[get cgjob_id($target) q] eq "q"} {
			set cgjob_id($target) $id
		}
	}
	foreach ptarget $ptargets {
		if {[get cgjob_id($ptarget) q] eq "q"} {
			set cgjob_ptargets($ptarget) $id
		}
	}
}

proc job_process_checkptargetsfinished {} {
	global cgjob_ptargets
	set done 0
	foreach ptarget [array names cgjob_ptargets] {
		set jobid $cgjob_ptargets($ptarget)
		if {[isint $jobid] && [catch {exec qstat -j $jobid}]} {
			unset cgjob_ptargets($ptarget)
			incr done
		}
	}
	return $done
}

proc job_process_sge {} {
	global cgjob
	while 1 {
		job_process_sge_onepass
		if {![llength $cgjob(queue)]} break
		# The queue may be not empty if jobs with deps on ptargets were present
		# only if (a) original job with ptargets is finished can we proceed
		after 1000
		while {![job_process_checkptargetsfinished]} {
			after 1000
		}
	}
}

proc job_process_sge_onepass {} {
	global cgjob cgjob_id cgjob_running cgjob_pid cgjob_ptargets cgjob_blocked
	set queue $cgjob(queue)
	update
	if {![llength $queue]} return
	set cgjob(queue) {}
	set jobroot [pwd]
	while {[llength $queue]} {
		set line [list_shift queue]
		foreach {jobid jobname job_logdir pwd deps foreach ftargetvars ftargets fptargets fskip code submitopts} $line break
		cd $pwd
		set job [job_logname $job_logdir $jobname]
		file mkdir [file dir $job]
		# If this job was previously blocked because of ptargets deps,
		# the ptargets set to stop further processing are cleared here
		# (They can still be reapplied later if they depend on ptargets that are not finished yet)
		if {[info exists cgjob_blocked($job)]} {
			foreach ptarget [job_targetsreplace $ftargets {}] {
				unset cgjob_ptargets($ptarget)
			}
			unset cgjob_blocked($job)
		}
		# check foreach deps, skip if not fullfilled
		# check for foreach patterns, expand into one ore more entries in the queue
		if {[llength $foreach]} {
			if {[catch {job_finddeps $job $foreach ftargetvars 1 fids time} fadeps]} {
				if {[regexp {^missing dependency} $fadeps]} {
					job_log $job "$fadeps"
				} elseif {[regexp {^ptargets hit} $fadeps]} {
					# ptarget dependency means the job cannot yet be submitted
					# nor any of the jobs that may be dependent on it
					# we wil re-add it to the queue for later processing
					# and create ptargets from its targets
					# job_logclear $job
					job_log $job "blocking at $jobname: $fadeps"
					lappend cgjob(queue) $line
					set cgjob_blocked($job) 1
					job_process_sge_marktargets {} [job_targetsreplace $ftargets {}] q
					continue
				} else {
					job_log $job "error in foreach dependencies for $jobname: $fadeps"
				}
				job_log $job "job $jobname skipped: dependencies not found"
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
			continue
		}
		job_lognf $job "==================== $jobname ===================="
		cd $pwd
		# check if job is already running, if so, mark targets with jobid
		set jobnum [job_process_sge_jobid $job]
		if {[isint $jobnum]} {
			job_process_sge_marktargets [file_read $job.targets] [file_read $job.ptargets] $jobnum
			job_log $job "job $jobname is already running, skip"
			continue
		}
		# check deps, skip if not fullfilled
		if {[catch {job_finddeps $job $deps newtargetvars 0 ids time $ftargetvars} adeps]} {
			# dependencies not found (or error) -> really skip job
			if {[regexp {^missing dependency} $adeps]} {
				job_log $job "$adeps"
			} elseif {[regexp {^ptargets hit} $adeps]} {
				# ptarget dependency means the job cannot yet be submitted
				# nor any of the jobs that may be dependent on it
				# we wil re-add it to the queue for later processing
				# and create ptargets from its targets
				# job_logclear $job
				job_log $job "blocking at $jobname: $adeps"
				lappend cgjob(queue) $line
				set cgjob_blocked($job) 1
				job_process_sge_marktargets {} [job_targetsreplace $ftargets {}] q
				continue
			} else {
				job_log $job "error in dependencies for $jobname: $adeps"
			}
			job_log $job "job $jobname skipped: dependencies not found"
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
			if {[llength $skip] && [job_checktargets $job $skip $time running]} {
				job_log $job "skipping $jobname: skip targets already completed or running"
				continue
			}
		}
		# check targets, if already done or running, skip
		set targets [job_targetsreplace $ftargets $targetvars]
		file_write $job.targets $targets
		set newtargets 0
		if {![job_checktargets $job $targets $time targetsrunning]} {
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
			foreach target $targets {
				# job_backup $target 1
				file delete $target
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
		append cmd [job_generate_code $job $pwd $adeps $targetvars $targets $ptargets $code]\n
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
	catch {file delete $job.finished $job.failed}
	set jnum [exec qsub -N j$name -q all.q -o $job.out -e $job.err {*}$options $runfile]
	regexp {[0-9]+} $jnum jobnum
	return $jobnum
}

proc job_process_sge_wait {} {
}

