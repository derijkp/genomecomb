proc job_process_init_status {} {
	global curjobnum jobsrunning graph curgraphid graphid totalduration
	# we will use par (parallel) code with some specifics for status
	if {[info commands job_process_par] eq ""} {auto_load job_process_par}
	set ::job_method_info {}
	interp alias {} job_process {} job_process_parstatus
	interp alias {} job_runall {} job_runall_par
	interp alias {} job_running {} job_running_status
	interp alias {} job_wait {} job_wait_status
	interp alias {} job_process_submit_par {} job_process_submit_status
	# interp alias {} job_process_par_checkjobid {} job_process_checkjobid_status
	set jobsrunning {}
	set totalduration {0 0}
	set graph {}
	set curjobnum 0
	set curgraphid 0
	unset -nocomplain graphid
	puts "status\tjobname\tjobnum\tduration\tmsg\tjob"
}

proc job_running_status {job} {
	return 1
}

proc job_process_pargraph {job jobname status duration checkcompressed adeps ids targets} {
	global graph job_name graphid curgraphid graphdone
#	set adeps [list_remdup $adeps]
#	set ids [list_remdup $ids]
#	set targets [list_remdup $targets]
	if {![info exists graphid($job)]} {
		set args {}
		set graphid($job) [incr curgraphid]
		set label [file tail $jobname]
		if {$duration ne "" && ![inlist {skipped error} $status]} {
			append label "\\n$duration"
		}
		if {$status eq "running"} {
			lappend args "color=green"
		} elseif {$status eq "skipped"} {
			lappend args "color=orange"
		} elseif {$status ne "ok"} {
			lappend args "color=red"
		}
		lappend args "label=\"$label\""
		lappend graph "$graphid($job) \[[join $args ,]\]"
	}
	set jobid $graphid($job)
	foreach file $adeps {
		if {$file eq ""} continue
		if {![info exists graphid($file)]} {
			set graphid($file) [incr curgraphid]
			if {![gzexists $file $checkcompressed]} {set color ",color=red"} else {set color {}}
			lappend graph "$graphid($file) \[shape=box$color,label=\"[file tail $file]\"\]"
		}
		lappend graph "$graphid($file) -> $jobid"
	}
	foreach id $ids {
		if {$id eq ""} continue
		set id $job_name($id)
		if {![info exists graphdone($graphid($id),$jobid)]} {
			lappend graph "$graphid($id) -> $jobid \[color=grey\]"
			set graphdone($graphid($id),$jobid) 1
		}
	}
	foreach file $targets {
		if {$file eq ""} continue
		if {![info exists graphid($file)]} {
			set graphid($file) [incr curgraphid]
			if {![gzexists $file $checkcompressed]} {set color ",color=red"} else {set color {}}
			lappend graph "$graphid($file) \[shape=box$color,label=\"[file tail $file]\"\]"
		}
		lappend graph "$jobid -> $graphid($file)"
	}
}

proc patternglob {pattern checkcompressed} {
	if {[string index $pattern 0] eq "^" && [string index $pattern end] eq "\$"} {
		set pattern [string range $pattern 1 end-1]
		set pattern [file_absolute $pattern]
		set glob [regexp2glob $pattern]
		set files {}
		if {$checkcompressed} {
			set list [bsort [gzfiles $glob]]
		} else {
			set list [bsort [glob $glob]]
		}
		foreach file $list {
			if {[regexp ^$pattern\$ $file]} {
				lappend files $file
			}
		}
	} else {
		set pattern [file_absolute $pattern]
		if {$checkcompressed} {
			set files [bsort [gzfiles $pattern]]
		} else {
			set files [bsort [glob $pattern]]
		}
	}
	return $files
}

proc job_process_parstatus {} {
	global cgjob cgjob_id cgjob_running curjobid curjobnum jobsrunning graph job_name totalduration
	set queue $cgjob(queue)
	# join [list_subindex $queue 1] \n
	update
	if {![llength $queue]} return
	set cgjob(queue) {}
	set jobroot [pwd]
	while {[llength $queue]} {
		set line [list_shift queue]
		foreach {jobid jobname job_logdir pwd deps ftargetvars ftargets fskip checkcompressed code submitopts frmtargets precode jobforce optional cores} $line break
		cd $pwd
		set job [job_logname $job_logdir $jobname]
		set submittime ""
		# get job log information -> duration
		set duration {}
		if {[file_or_link_exists [job.file log $job]]} {
			set jobloginfo [job_parse_log $job]
			foreach {status time endtime run duration submittime} $jobloginfo break
			if {$time eq ""} {unset time}
			if {$endtime eq ""} {unset endtime}
		}
		# check if job is already running, if so, mark targets with jobid
		set jobnum [job_process_par_jobid $job]
		if {[isint $jobnum]} {
			set temptargets [job_getfromlog $job cgjobinfo_targets]
			set temprmtargets [job_getfromlog $job cgjobinfo_rmtargets]
			job_process_par_marktargets $temptargets $temprmtargets $jobnum
			puts "running\t$jobname\t$jobnum\t$duration\t\t$job"
			job_logfile_add $job $jobnum running $ftargets $cores
			lappend jobsrunning $jobnum
			catch {job_finddeps $job $deps newtargetvars 0 ids time timefile $checkcompressed $ftargetvars} adeps
			set job_name($jobnum) $job
			job_process_pargraph $job $jobname running $duration $checkcompressed $adeps $ids $temptargets
			continue
		} else {
			set jobnum [incr curjobnum]
			set job_name($jobnum) $job
		}
		# check deps, skip if not fullfilled
		if {[catch {job_finddeps $job $deps newtargetvars 0 ids time timefile $checkcompressed $ftargetvars} adeps]} {
			# dependencies not found (or error) -> really skip job
			if {[regexp {^missing dependency} $adeps]} {
				if {!$optional} {
					job_log "could not run job $jobname:\n$adeps"
				}
#				job_log $job "$adeps"
			} else {
				puts "joberror\t$jobname\t\terror in dependencies: $adeps\t$job"
			}
			puts "skipped\t$jobname\t\tdependencies not found\t$job"
			job_logfile_add $job . skipped $ftargets $cores "dependencies not found" $submittime
			job_logclose $job
			continue
		}
		set targetvars $ftargetvars
		lappend targetvars {*}$newtargetvars
		if {$jobforce} {job_log $job "forced $jobname"}
		if {$jobforce || $cgjob(force)} {set time force}
		# check targets, if already done or running, skip
		if {$ftargets ne ""} {
			set targets [job_targetsreplace $ftargets $targetvars]
			job_log $job "cgjobinfo_targets: [list $targets]"
			set newtargets 0
			set time 0
			set timefile {}
			if {!$jobforce && ![job_checktargets $job $targets 0 $time $timefile $checkcompressed targetsrunning]} {
				set newtargets 1
			}
		} else {
			set targets {}
			job_log $job "cgjobinfo_targets: [list $targets]"
			set targetsrunning {}
			set newtargets 1
		}
		if {$frmtargets ne ""} {
			set rmtargets [job_targetsreplace $frmtargets $targetvars]
			job_log $job "cgjobinfo_rmtargets: [list $rmtargets]"
		} else {
			set rmtargets {}
		}
		# indicate targets
		set job_err [job.file err $job]
		job_process_par_marktargets $targets $rmtargets $jobnum
		if {[file_or_link_exists $job_err]} {
			puts "error\t$jobname\t$jobnum\t\terror file available\t$job_err"
			job_logfile_add $job $jobnum error $targets $cores [job_cleanmsg [file_read $job_err]] $submittime
			job_process_pargraph $job $jobname error $duration $checkcompressed $adeps $ids $targets
		} elseif {!$newtargets} {
			puts "ok\t$jobname\t$jobnum\t$duration\ttargets found\t$job"
			job_logfile_add $job $jobnum finished $targets $cores "" $submittime
			job_process_pargraph $job $jobname ok $duration $checkcompressed $adeps $ids $targets
			job_logclose $job
			continue
		} else {
			puts "wrong\t$jobname\t$jobnum\t\ttargets not ok, no error file\t$job"
			job_logfile_add $job $jobnum error $targets $cores "some error, no errorfile was left" $submittime
			job_process_pargraph $job $jobname wrong $duration $checkcompressed $adeps $ids $targets
		}
		set cgjob_running($job) $jobnum
		job_logclose $job
	}
	cd $jobroot
}

proc job_wait_status {} {
	global cgjob jobsrunning graph totalduration
	job_logfile_par_close
	if {[file exists $cgjob(logfile).running]} {
		job_update $cgjob(logfile).running $cgjob(cleanup) 0 $cgjob(removeold)
	}
	puts \n
	set tduration [time_format [list 0 [lindex $totalduration end]] "%H %M %S %s"]
	regexp {[1-9]*[0-9]$} [lindex $tduration 0] hours
	set duration "[expr {24*[lindex $totalduration 0]+$hours}]:[lindex $tduration 1]:[lindex $tduration 2].[lindex $tduration 3]"
	puts "total single thread duration: $duration"
	if {[llength $jobsrunning]} {puts "running [join $jobsrunning ,]"}
	file_write jobdeps.dot "digraph deps \{\n[join $graph \;\n]\;\n\}\n"
	exec dot -Tps jobdeps.dot -o jobdeps.ps
	puts "Written jobdeps.dot and jobdeps.ps"
}
