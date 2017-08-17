proc job_process_par_jobid {job} {
	if {![job_file_exists $job.jobid]} {
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

proc job_process_par_marktargets {targets ptargets rmtargets id} {
	global cgjob_id cgjob_ptargets cgjob_rm
	foreach target $targets {
		if {[get cgjob_id($target) q] eq "q"} {
			set gzfile [gzfiles $target]
			if {$gzfile ne ""} {set target $gzfile}
			set cgjob_id($target) $id
		}
	}
	foreach ptarget $ptargets {
		if {[get cgjob_ptargets($ptarget) q] eq "q"} {
			set cgjob_ptargets($ptarget) $id
		}
	}
	foreach rmtarget $rmtargets {
		set cgjob_rm($rmtarget) $id
	}
}

proc job_process_par_checkptargetsfinished {} {
	global cgjob_ptargets
	set done 0
	foreach ptarget [array names cgjob_ptargets] {
		set jobid $cgjob_ptargets($ptarget)
		if {![job_running $jobid]} {
			unset -nocomplain cgjob_ptargets($ptarget)
			incr done
		}
	}
	return $done
}

proc job_process_par_onepass {} {
	global cgjob cgjob_id cgjob_running cgjob_ptargets cgjob_blocked
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
		foreach {jobid jobname job_logdir pwd deps foreach ftargetvars ftargets fptargets fskip checkcompressed code submitopts frmtargets precode jobforce optional} $line break
		cd $pwd
		set job [job_logname $job_logdir $jobname]
		file mkdir [file dir $job]
		set time 0
		set timefile {}
		# If this job was previously blocked because of ptargets deps,
		# the ptargets set to stop further processing are cleared here
		# (They can still be reapplied later if they depend on ptargets that are not finished yet)
		if {[info exists cgjob_blocked($job)]} {
			foreach ptarget [job_targetsreplace $ftargets {}] {
				unset -nocomplain cgjob_ptargets($ptarget)
			}
			unset cgjob_blocked($job)
		}
		# check foreach deps, skip if not fullfilled
		# check for foreach patterns, expand into one ore more entries in the queue
		set submittime [timestamp]
		if {[llength $foreach]} {
			if {[catch {job_finddeps $job $foreach ftargetvars 1 fids time timefile $checkcompressed} fadeps]} {
				if {[regexp {^missing dependency} $fadeps]} {
					job_log $job "$fadeps"
					job_logfile_add $job . error $ftargets $errormsg $submittime
				} elseif {[regexp {^ptargets hit} $fadeps]} {
					# ptarget dependency means the job cannot yet be submitted
					# nor any of the jobs that may be dependent on it
					# we wil re-add it to the queue for later processing
					# and create ptargets from its targets
					# job_logclear $job
					job_log $job "blocking at $jobname: $fadeps"
					lappend cgjob(queue) $line
					set cgjob_blocked($job) 1
					job_process_par_marktargets {} [job_targetsreplace $ftargets {}] {} q
					job_logclose $job
					continue
				} else {
					set errormsg "error in foreach dependencies for $jobname: $fadeps"
					job_log $job $errormsg
					job_logfile_add $job . error $ftargets $errormsg $submittime
				}
				job_log $job "-----> job $jobname skipped: dependencies not found"
				job_logclose $job
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
			job_logclear $job
			continue
		}
		job_lognf $job "==================== $jobname ===================="
		cd $pwd
		# check if job is already running, if so, mark targets with jobid
		set jobnum [job_process_par_jobid $job]
		if {[isint $jobnum]} {
			if {[job_file_exists $job.targets]} {
				set temptargets [file_read $job.targets]
			} else {
				set temptargets {}
			}
			if {[job_file_exists $job.rmtargets]} {
				set temprmtargets [file_read $job.rmtargets]
			} else {
				set temprmtargets {}
			}
			if {![job_file_exists $job.ptargets]} {
				error "job $job ($jobnum) seems to be running, but there is no $job.ptargets"
			}
			if {[job_file_exists $job.ptargets]} {
				set tempptargets [file_read $job.ptargets]
			} else {
				set tempptargets {}
			}
			job_process_par_marktargets $temptargets $tempptargets $temprmtargets $jobnum
			job_log $job "job $jobname is already running, skip"
			job_logfile_add $job $jobnum running $ftargets
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
					job_logfile_add $job . skipped $ftargets "dependencies not found" $submittime
					job_logclose $job
					continue
				}
			} elseif {[regexp {^ptargets hit} $adeps]} {
				# ptarget dependency means the job cannot yet be submitted
				# nor any of the jobs that may be dependent on it
				# we wil re-add it to the queue for later processing
				# and create ptargets from its targets
				# job_logclear $job
				job_log $job "blocking at $jobname: $adeps"
				lappend cgjob(queue) $line
				set cgjob_blocked($job) 1
				job_process_par_marktargets {} [job_targetsreplace $ftargets {}] {} q
				job_logclose $job
				continue
			} else {
				if {!$optional && !$cgjob(skipjoberrors)} {
					set joberror "error in dependencies for job $jobname:\n$adeps"
				} else {
					job_log $job "error in dependencies for $jobname: $adeps"
					job_log $job "-----> job $jobname skipped: dependencies not found"
					job_logfile_add $job . skipped $ftargets "error in dependencies: $adeps" $submittime
					job_logclose $job
					continue
				}
			}
		}
		file_write $job.deps $adeps
		set targetvars $ftargetvars
		lappend targetvars {*}$newtargetvars
		if {$jobforce} {job_log $job "forcing $jobname"}
		if {$cgjob(force)} {set time force}
		# check skip targets, if already done or running, skip
		if {!$jobforce && !$cgjob(force) && [llength $fskip]} {
			set doskip 0
			foreach skip $fskip {
				set skip [job_targetsreplace $skip $targetvars]
				if {[llength $skip] && [job_checktargets $job $skip $time $timefile $checkcompressed running]} {
					set doskip 1
					break
				}
			}
			if {$doskip} {
				job_log $job "skipping $jobname: skip targets already completed or running"
				job_logfile_add $job . skipped $ftargets "skip targets already completed or running" $submittime
				job_logclose $job
				continue
			}
		}
		# check targets, if already done or running, skip
		if {$ftargets ne ""} {
			set targets [job_targetsreplace $ftargets $targetvars]
			file_write $job.targets $targets
			set newtargets 0
			if {$jobforce || ![job_checktargets $job $targets $time $timefile $checkcompressed targetsrunning]} {
				set newtargets 1
			}
		} else {
			set targets {}
			file_write $job.targets {}
			set targetsrunning {}
			set newtargets 1
		}
		if {$frmtargets ne ""} {
			set rmtargets [job_targetsreplace $frmtargets $targetvars]
			file_write $job.rmtargets $rmtargets
		} else {
			set rmtargets {}
		}
		set ptargets [job_targetsreplace $fptargets $targetvars]
		file_write $job.ptargets $ptargets
		if {[llength $ptargets] && ![llength [job_findptargets $ptargets $checkcompressed]]} {
			set newtargets 1
		}
		# indicate targets are in the queue, so job_finddeps will find them
		job_process_par_marktargets $targets $ptargets $rmtargets q
#		# if deps or overlapping targets are not finished yet, put line back in the queue and skip running
#		if {$depsrunning || [llength $targetsrunning]} {
#			job_logclear $job
#			lappend cgjob(queue) $line
#			job_logclose $job
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
		if {!$jobforce && !$newtargets} {
			job_log $job "skipping $jobname: targets already completed or running"
			job_logfile_add $job . skipped $targets "targets already completed or running" $submittime
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
		append cmd "exec $cgjob(runcmd) \"\$0\" \"\$@\"\n"
		append cmd [job_generate_code $job $pwd $adeps $targetvars $targets $ptargets $checkcompressed $code]\n
		append cmd "file_add \{$job.log\} \"\[job_timestamp\]\\tending $jobname\"\n"
		set runfile $job.run
		file_write $runfile $cmd
		file attributes $runfile -permissions u+x
		# submit job
		set ids [list_remove $ids {}]
		set jobnum [job_process_par_submit $job $runfile -deps $ids {*}$submitopts]
		job_log $job "-------------------- submitted $jobname ($jobnum <- $ids) (run $currentrun) --------------------"
		job_logfile_add $job $jobnum submitted $targets "" $submittime
		job_logclose $job
		file_write $job.jobid $jobnum
		job_process_par_marktargets $targets $ptargets $rmtargets $jobnum
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
		# The queue may be not empty if jobs with deps on ptargets were present
		# only if (a) original job with ptargets is finished can we proceed
		after 1000
		while {![job_process_par_checkptargetsfinished]} {
			after 1000
		}
	}
}

proc job_logfile_par_close {} {
	global cgjob
	if {![info exists cgjob(f_logfile)]} return
	puts $cgjob(f_logfile) [join [list total . running $cgjob(starttime) "" "" "" "" "" ""] \t]
	close $cgjob(f_logfile)
	file rename $cgjob(logfile).submitting $cgjob(logfile).running
	job_update $cgjob(logfile).running $cgjob(cleanup) 0 $cgjob(removeold)
}

