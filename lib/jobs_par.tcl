proc job_process_par_jobid {job} {
	if {![job_file_or_link_exists [job.file jid $job]]} {
		return ""
	}
	set jobid [file_read [job.file jid $job]]
	if {[catch {exec qstat -j $jobid}]} {
		file delete [job.file jid $job]
		return ""
	} else {
		return $jobid
	}
}

proc job_process_par_marktargets {targets rmtargets id} {
	global cgjob_id cgjob_rm job_getinfo cgjob_getinfo
	foreach target $targets {
		if {[get cgjob_id($target) q] eq "q"} {
			set gzfile [gzfiles $target]
			if {$gzfile ne ""} {set target $gzfile}
			if {[get job_getinfo 0] && ![info exists cgjob_id($target)]} {lappend cgjob_getinfo(id) $target}
			set cgjob_id($target) $id
		}
	}
	foreach rmtarget $rmtargets {
		if {[get job_getinfo 0] && ![info exists cgjob_rm($rmtarget)]} {lappend cgjob_getinfo(rm) $rmtarget}
		set cgjob_rm($rmtarget) $id
	}
}

proc job_timecmd {job} {
	global job_timecmd
	if {![info exists job_timecmd]} {
		if {[catch {
			set job_timecmd [exec which time]
		}]} {
			set job_timecmd {}
		}
	}
	return "\"$job_timecmd\" -o \"[job.file tm $job]\" -v"
}

proc job_process_par_onepass {} {
	global cgjob cgjob_id cgjob_running
	set currentrun [file tail [get cgjob(logfile) ""]]
	set queue $cgjob(queue)
	# join [list_subindex $queue 1] \n
	update
	if {![llength $queue]} return
	set cgjob(queue) {}
	set jobroot [pwd]
	while {[llength $queue]} {
		set joberror {}
		set line [list_shift queue]
		foreach {jobid jobname job_logdir pwd deps ftargetvars ftargets fskip checkcompressed code submitopts frmtargets precode jobforce optional cores} $line break
		cd $pwd
		set job [job_logname $job_logdir $jobname]
		file mkdir [file dir $job]
		set time 0
		set timefile {}
		set submittime [timestamp]
		job_lognf $job "==================== $jobname ===================="
		cd $pwd
		# check if job is already running, if so, mark targets with jobid
		set jobnum [job_process_par_jobid $job]
		if {[isint $jobnum]} {
			set temptargets [job_getfromlog $job cgjobinfo_targets]
			set temprmtargets [job_getfromlog $job cgjobinfo_rmtargets]
			job_process_par_marktargets $temptargets $temprmtargets $jobnum
			job_log $job "job $jobname is already running, skip"
			job_logfile_add $job $jobnum running $ftargets $cores
			job_logclose $job
			continue
		}
		# check deps, skip if not fullfilled
		if {[catch {job_finddeps $job $deps newtargetvars 0 ids time timefile $checkcompressed $ftargetvars} adeps]} {
			# dependencies not found (or error) -> really skip job
			if {[regexp {^missing dependency} $adeps]} {
				job_log $job "$adeps"
				if {!$optional && !$cgjob(skipjoberrors)} {
					set joberror "error trying to run job $jobname:\n$adeps"
				} else {
					job_log $job "-----> job $jobname skipped: dependencies not found"
					job_logfile_add $job . skipped $ftargets $cores "dependencies not found" $submittime
					job_logclose $job
					continue
				}
			} else {
				if {!$optional && !$cgjob(skipjoberrors)} {
					set joberror "error in dependencies for job $jobname:\n$adeps"
				} else {
					job_log $job "error in dependencies for $jobname: $adeps"
					job_log $job "-----> job $jobname skipped: dependencies not found"
					job_logfile_add $job . skipped $ftargets $cores "error in dependencies: $adeps" $submittime
					job_logclose $job
					continue
				}
			}
		}
		job_log $job "cgjobinfo_deps: [list $adeps]"
		set targetvars $ftargetvars
		lappend targetvars {*}$newtargetvars
		if {$jobforce} {job_log $job "forcing $jobname"}
		if {$cgjob(force)} {set time force}
		# check skip targets, if already done or running, skip
		if {!$jobforce && !$cgjob(force) && [llength $fskip]} {
			set doskip 0
			foreach skip $fskip {
				set skip [job_targetsreplace $skip $targetvars]
				if {[llength $skip] && [job_checktargets $job $skip 1 $time $timefile $checkcompressed running]} {
					set doskip 1
					break
				}
			}
			if {$doskip} {
				job_log $job "skipping $jobname: skip targets already completed or running"
				job_logfile_add $job . skipped $ftargets $cores "skip targets already completed or running" $submittime
				job_logclose $job
				continue
			}
		}
		# check targets, if already done or running, skip
		if {$ftargets ne ""} {
			set targets [job_targetsreplace $ftargets $targetvars]
			job_log $job "cgjobinfo_targets: [list $targets]"
			set newtargets 0
			if {$jobforce || ![job_checktargets $job $targets 0 $time $timefile $checkcompressed targetsrunning]} {
				set newtargets 1
			}
		} else {
			set targets {}
			job_log $job "cgjobinfo_targets: {}"
			set targetsrunning {}
			set newtargets 1
		}
		if {$frmtargets ne ""} {
			set rmtargets [job_targetsreplace $frmtargets $targetvars]
			job_log $job "cgjobinfo_rmtargets: [list $rmtargets]"
		} else {
			set rmtargets {}
		}
		# indicate targets are in the queue, so job_finddeps will find them
		job_process_par_marktargets $targets $rmtargets q
		job_log $job
		if {$cgjob(force)} {
			foreach target $targets {
				# job_backup $target 1
				file delete $target
			}
			set newtargets 1
		}
		if {!$jobforce && !$newtargets} {
			job_log $job "skipping $jobname: targets already completed or running"
			job_logfile_add $job . skipped $targets $cores "targets already completed or running" $submittime
			job_logclose $job
			continue
		}
		if {$joberror ne ""} {error $joberror}
		# job_log $job "-------------------- submitting $jobname --------------------"
		# run code
		set cmd {#!/bin/sh}
		append cmd \n
		append cmd {#$ -S /bin/bash} \n
		append cmd {#$ -V} \n
		append cmd {#$ -cwd} \n
		append cmd "\n\# the next line restarts using runcmd (specialised tclsh) \\\n"
		append cmd "exec [job_timecmd $job] $cgjob(runcmd) \"\$0\" \"\$@\"\n"
		append cmd [job_generate_code $job $pwd $adeps $targetvars $targets $checkcompressed $code]\n
		append cmd "file_add \{[job.file log $job]\} \"\[job_timestamp\]\\tending $jobname\"\n"
		set runfile [job.file run $job]
		file_write $runfile $cmd
		file attributes $runfile -permissions u+x
		# submit job
		set ids [list_remove [list_remdup $ids] {}]
		set jobnum [job_process_submit_par $job $runfile -deps $ids {*}$submitopts]
		job_log $job "-------------------- submitted $jobname ($jobnum <- $ids) (run $currentrun) --------------------"
		job_logfile_add $job $jobnum submitted $targets $cores "" $submittime
		file_write [job.file jid $job] $jobnum
		job_log $job "cgjobinfo_jobid: [list $jobid]"
		job_logclose $job
		job_process_par_marktargets $targets $rmtargets $jobnum
		set cgjob_running($job) $jobnum
	}
	lappend cgjob(queue) {*}$queue
	cd $jobroot
}

proc job_process_par {} {
	global cgjob
	while 1 {
		job_process_par_onepass
		if {![llength $cgjob(queue)]} break
		after 1000
	}
}

proc job_process_par_checkrunning {} {
	global cgjob_running
	set running 0
	foreach job [array names cgjob_running] {
		set jobid $cgjob_running($job)
		if {![job_running $jobid]} {
			unset -nocomplain cgjob_running($job)
		} else {
			incr running
		}
	}
	return $running
}

proc job_runall_par {} {
	global cgjob
	while 1 {
		job_process_par_onepass
		# empty queue
		after 1000
		while {[job_process_par_checkrunning]} {
			after 1000
		}
		if {![llength $cgjob(queue)]} break
	}
}

proc job_logfile_par_close {} {
	global cgjob
	if {![info exists cgjob(f_logfile)]} return
	puts $cgjob(f_logfile) [join [list total . running $cgjob(starttime) "" "" "" "" "" "" "" ""] \t]
	close $cgjob(f_logfile)
	file rename $cgjob(logfile).submitting $cgjob(logfile).running
	# job_update must be managed by the actual close (e.g. job_logfile_distr_close)
	# job_update $cgjob(logfile).running $cgjob(cleanup) 1 $cgjob(removeold)
}

