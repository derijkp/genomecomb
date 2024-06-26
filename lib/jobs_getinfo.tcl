# allows you to get deps and targets of a command (set of commands) without actually running/submitting:
# job_getinfo 1
# ... commands ...
# foreach {deps targets} [job_getinfo 0] break
proc job_getinfo {{value {}}} {
	if {$value eq ""} {
		get ::job_getinfo 0
	} elseif {$value eq "1"} {
		global cgjob_getinfo
		set ::job_getinfo $value
		set cgjob_getinfo(deps) {}
		set cgjob_getinfo(targets) {}
	} elseif {$value eq "0"} {
		global cgjob_id cgjob_rm job_getinfo cgjob_getinfo
		foreach target [get cgjob_getinfo(id) ""] {unset -nocomplain cgjob_id($target)} 
		foreach target [get cgjob_getinfo(rm) ""] {unset -nocomplain cgjob_rm($target)} 
		foreach job [get cgjob_getinfo(jobs) ""] {
			catch {close $cgjob(f,$job)}
			unset -nocomplain cgjob(f,$job)
			unset -nocomplain cgjob(buffer,$job)
		}
		set ::job_getinfo $value
		set result [list [list_remove [list_remdup $cgjob_getinfo(deps)] {}] [list_remove [list_remdup $cgjob_getinfo(targets)] {}]]
		unset cgjob_getinfo
		return $result
	} else {
		error "job_getinfo error: value $value not supported"
	}
}

proc job_process_getinfo {jobid jobname job_logdir pwd deps ftargetvars ftargets fskip checkcompressed code submitopts frmtargets precode jobforce optional cores} {
# putsvars jobid jobname job_logdir pwd deps ftargetvars ftargets fskip checkcompressed code submitopts frmtargets precode jobforce optional
	global cgjob job_deps
	global cgjob_getinfo
	set jobroot [pwd]
	set joberror {}
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
			return
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
					return
				}
			} else {
				if {!$optional && !$cgjob(skipjoberrors)} {
					set joberror "error in dependencies for job $jobname:\n$adeps"
				} else {
					job_log $job "error in dependencies for $jobname: $adeps"
					job_log $job "-----> job $jobname skipped: dependencies not found"
					job_logfile_add $job . skipped $ftargets $cores "error in dependencies: $adeps" $submittime
					job_logclose $job
					return
				}
			}
		}
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
				job_log $job "skipping $jobname: skip targets ($skip) already completed or running"
				job_logfile_add $job . skipped $ftargets $cores "skip targets ($skip) already completed or running" $submittime
				job_logclose $job
				return
			}
		}
		# check targets, if already done or running, skip
		if {$ftargets ne ""} {
			set targets [job_targetsreplace $ftargets $targetvars]
			set newtargets 0
			if {$jobforce || ![job_checktargets $job $targets 0 $time $timefile $checkcompressed targetsrunning]} {
				set newtargets 1
			}
		} else {
			set targets {}
			set targetsrunning {}
			set newtargets 1
		}
		if {$frmtargets ne ""} {
			set rmtargets [job_targetsreplace $frmtargets $targetvars]
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
			return
		}
		if {$joberror ne ""} {
			job_logclose $job
			error $joberror
		}
		job_logclose $job
	# compile info
#set o [open ~/tmp/temp a]
#puts $o "===== [file tail $job] ====="
#foreach var {job adeps targets rmtargets cgjob_getinfo(deps) cgjob_getinfo(targets)} {
#puts $o "----- $var -----\n[join [get $var] \n]\n"
#}
#close $o
	lappend cgjob_getinfo(jobs) $job
	list_addnew cgjob_getinfo(deps) {*}[list_lremove $adeps $cgjob_getinfo(targets)]
	if {$rmtargets ne ""} {
		set cgjob_getinfo(targets) [list_lremove [list_remdup $cgjob_getinfo(targets)] $rmtargets]
	}
	list_addnew cgjob_getinfo(targets) {*}$targets
	cd $jobroot
}

