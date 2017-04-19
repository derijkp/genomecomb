proc job_process_status_init {} {
	global curjobnum jobsrunning graph curgraphid graphid totalduration
	# we will use par (parallel) code with some specifics for status
	if {[info commands job_process_par] eq ""} {auto_load job_process_par}
	set ::job_method_info {}
	interp alias {} job_process {} job_process_parstatus
	interp alias {} job_running {} job_running_status
	interp alias {} job_wait {} job_process_status_wait
	interp alias {} job_process_par_submit {} job_process_status_submit
	interp alias {} job_process_par_checkjobid {} job_process_status_checkjobid
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

proc job_process_pargraph {job jobname status duration checkcompressed adeps ids targets ptargets} {
	global graph job_name graphid curgraphid graphdone
#	set adeps [list_remdup $adeps]
#	set ids [list_remdup $ids]
#	set targets [list_remdup $targets]
#	set ptargets [list_remdup $ptargets]
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
	foreach pattern $ptargets {
		if {$pattern eq ""} continue
		foreach file [patternglob $pattern $checkcompressed] {
			if {![info exists graphid($file)]} {
				set graphid($file) [incr curgraphid]
				if {![gzexists $file $checkcompressed]} {set color ",color=red"} else {set color {}}
				lappend graph "$graphid($file) \[shape=box$color,label=\"[file tail $file]\"\]"
			}
			lappend graph "$jobid -> $graphid($file) \[style=dashed\]"
		}
	}
}

proc patternglob {pattern checkcompressed} {
	if {[string index $pattern 0] eq "^" && [string index $pattern end] eq "\$"} {
		set pattern [string range $pattern 1 end-1]
		set pattern [file_absolute $pattern]
		set glob [regexp2glob $pattern]
		set files {}
		if {$checkcompressed} {
			set list [lsort -dict [gzfiles $glob]]
		} else {
			set list [lsort -dict [glob $glob]]
		}
		foreach file $list {
			if {[regexp ^$pattern\$ $file]} {
				lappend files $file
			}
		}
	} else {
		set pattern [file_absolute $pattern]
		if {$checkcompressed} {
			set files [lsort -dict [gzfiles $pattern]]
		} else {
			set files [lsort -dict [glob $pattern]]
		}
	}
	return $files
}

proc job_process_parstatus {} {
	global cgjob cgjob_id cgjob_running cgjob_ptargets cgjob_blocked curjobid curjobnum jobsrunning graph job_name totalduration
	set queue $cgjob(queue)
	# join [list_subindex $queue 1] \n
	update
	if {![llength $queue]} return
	set cgjob(queue) {}
	set jobroot [pwd]
	while {[llength $queue]} {
		set line [list_shift queue]
		foreach {jobid jobname job_logdir pwd deps foreach ftargetvars ftargets fptargets fskip checkcompressed code submitopts frmtargets precode jobforce optional} $line break
		cd $pwd
		set job [job_logname $job_logdir $jobname]
#		# If this job was previously blocked because of ptargets deps,
#		# the ptargets set to stop further processing are cleared here
#		# (They can still be reapplied later if they depend on ptargets that are not finished yet)
#		if {[info exists cgjob_blocked($job)]} {
#			foreach ptarget [job_targetsreplace $ftargets {}] {
#				unset -nocomplain cgjob_ptargets($ptarget)
#			}
#			unset cgjob_blocked($job)
#		}
		# check foreach deps, skip if not fullfilled
		# check for foreach patterns, expand into one ore more entries in the queue
		if {[llength $foreach]} {
			# we assume ptarget locks are resolved
			unset -nocomplain cgjob_ptargets
			if {[catch {job_finddeps $job $foreach ftargetvars 1 fids time timefile $checkcompressed} fadeps]} {
				if {[regexp {^missing dependency} $fadeps]} {
#					job_log $job "$fadeps"
				} elseif {[regexp {^ptargets hit} $fadeps]} {
					error "cannot get status on unfinished ptargets hit"
				} else {
					puts "joberror\t$jobname\t\terror in foreach dependencies: $fadeps\t$job"
				}
				puts "skipped\t$jobname\t\tdependencies not found\t$job"
				lappend graph $jobname
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
		cd $pwd
		# get job log information -> duration
		set duration {}
		if {[job_file_exists $job.log]} {
			set jobloginfo [job_parse_log $job $totalduration]
			foreach {failed time endtime run duration totalduration submittime} $jobloginfo break
			if {$time eq ""} {unset time}
			if {$endtime eq ""} {unset endtime}
		}
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
			puts "running\t$jobname\t$jobnum\t$duration\t\t$job"
			lappend jobsrunning $jobnum
			catch {job_finddeps $job $deps newtargetvars 0 ids time timefile $checkcompressed $ftargetvars} adeps
			set job_name($jobnum) $job
			job_process_pargraph $job $jobname running $duration $checkcompressed $adeps $ids $temptargets $tempptargets
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
			} elseif {[regexp {^ptargets hit} $adeps]} {
			} else {
				puts "joberror\t$jobname\t\terror in dependencies: $adeps\t$job"
			}
			puts "skipped\t$jobname\t\tdependencies not found\t$job"
			continue
		}
		set targetvars $ftargetvars
		lappend targetvars {*}$newtargetvars
		if {$jobforce} {job_log $job "forced $jobname"}
		if {$jobforce || $cgjob(force)} {set time force}
		# check targets, if already done or running, skip
		if {$ftargets ne ""} {
			set targets [job_targetsreplace $ftargets $targetvars]
			file_write $job.targets $targets
			set newtargets 0
			set time 0
			set timefile {}
			if {!$jobforce && ![job_checktargets $job $targets $time $timefile $checkcompressed targetsrunning]} {
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
		# expand ptargets to available files
		set temp {}
		foreach pattern $ptargets {
			if {$pattern eq ""} continue
			foreach file [patternglob $pattern $checkcompressed] {
				lappend temp $file
			}
		}
		# indicate targets
		job_process_par_marktargets [list_concat $targets $temp] $ptargets $rmtargets $jobnum
		if {[job_file_exists $job.err]} {
			puts "error\t$jobname\t$jobnum\t\terror file available\t$job.err"
			job_process_pargraph $job $jobname error $duration $checkcompressed $adeps $ids $targets $ptargets
		} elseif {!$newtargets} {
			puts "ok\t$jobname\t$jobnum\t$duration\ttargets found\t$job"
			job_process_pargraph $job $jobname ok $duration $checkcompressed $adeps $ids $targets $ptargets
			continue
		} else {
			puts "wrong\t$jobname\t$jobnum\t\ttargets not ok, no error file\t$job"
			job_process_pargraph $job $jobname wrong $duration $checkcompressed $adeps $ids $targets $ptargets
		}
		set cgjob_running($job) $jobnum
	}
	cd $jobroot
}

proc job_process_status_wait {} {
	global jobsrunning graph totalduration
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
