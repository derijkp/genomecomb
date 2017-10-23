proc job_getinfo {{value {}}} {
	if {$value eq ""} {
		get ::job_getinfo 0
	} elseif {$value eq "1"} {
		global cgjob_getinfo
		set ::job_getinfo $value
		set cgjob_getinfo(deps) {}
		set cgjob_getinfo(targets) {}
	} elseif {$value eq "0"} {
		global cgjob_id cgjob_ptargets cgjob_rm job_getinfo cgjob_getinfo
		foreach target [get cgjob_getinfo(id) ""] {unset -nocomplain cgjob_id($target)} 
		foreach target [get cgjob_getinfo(ptargets) ""] {unset -nocomplain cgjob_ptargets($target)} 
		foreach target [get cgjob_getinfo(rm) ""] {unset -nocomplain cgjob_rm($target)} 
		foreach job [get cgjob_getinfo(jobs) ""] {
			unset -nocomplain cgjob(f,$job)
			unset -nocomplain cgjob(buffer,$job)
		}
		set ::job_getinfo $value
		set result [list $cgjob_getinfo(deps) $cgjob_getinfo(targets)]
		unset cgjob_getinfo
		return $result
	} else {
		error "job_getinfo error: value $value not supported"
	}
}

proc job_process_getinfo {jobid jobname job_logdir pwd deps foreach ftargetvars ftargets fptargets fskip checkcompressed code submitopts frmtargets precode jobforce optional} {
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
					job_logfile_add $job . skipped $ftargets "dependencies not found" $submittime
					job_logclose $job
					return
				}
			} else {
				if {!$optional && !$cgjob(skipjoberrors)} {
					set joberror "error in dependencies for job $jobname:\n$adeps"
				} else {
					job_log $job "error in dependencies for $jobname: $adeps"
					job_log $job "-----> job $jobname skipped: dependencies not found"
					job_logfile_add $job . skipped $ftargets "error in dependencies: $adeps" $submittime
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
				if {[llength $skip] && [job_checktargets $job $skip $time $timefile $checkcompressed running]} {
					set doskip 1
					break
				}
			}
			if {$doskip} {
				job_log $job "skipping $jobname: skip targets already completed or running"
				job_logfile_add $job . skipped $ftargets "skip targets already completed or running" $submittime
				job_logclose $job
				return
			}
		}
		# check targets, if already done or running, skip
		if {$ftargets ne ""} {
			set targets [job_targetsreplace $ftargets $targetvars]
			set newtargets 0
			if {$jobforce || ![job_checktargets $job $targets $time $timefile $checkcompressed targetsrunning]} {
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
		set ptargets [job_targetsreplace $fptargets $targetvars]
		if {[llength $ptargets] && ![llength [job_findptargets $ptargets $checkcompressed]]} {
			set newtargets 1
		}
		# indicate targets are in the queue, so job_finddeps will find them
		job_process_par_marktargets $targets $ptargets $rmtargets q
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
			list_addnew cgjob_getinfo(deps) {*}[list_lremove $adeps $cgjob_getinfo(targets)]
			list_addnew cgjob_getinfo(targets) {*}$targets
			return
		}
		if {$joberror ne ""} {error $joberror}
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

