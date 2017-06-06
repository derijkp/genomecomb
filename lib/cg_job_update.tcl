proc job_cleanlogs {logfile} {
	putslog "cleaning up log_jobs for $logfile"
	set f [open $logfile]
	while 1 {
		if {[gets $f line] == -1} break
		if {[string index $line 0] ne {#}} break
		if {[regexp {^# ([^:]+): (.*)$} $line temp key value]} {
			set cgjob($key) $value
		}
	}
	if {[split $line \t] ne {job jobid status submittime starttime endtime duration targets msg run}} {
		close $f
		error "file $logfile is not a proper logfile (must be tsv with fields: job jobid status submittime starttime endtime duration targets msg run)"
	}
	while 1 {
		if {[gets $f line] == -1} break
		set sline [split $line \t]
		foreach {jobo jobid status submittime starttime endtime duration targets msg run} $sline break
		if {$jobo eq "total"} continue
		if {[get cgjob(basedir) ""] ne "" && [file pathtype $jobo] ne "absolute"} {
			set job $cgjob(basedir)/$jobo
		} else {
			set job $jobo
		}
		set files [glob -nocomplain $job.*]
		if {![llength $files]} continue
		file delete {*}$files
		set dirsa([file dir $job]) 1
	}
	foreach dir [array names dirsa] {
		if {[catch {glob $dir/*}]} {
			catch {file delete $dir}
		}
	}
}

# status can be
# submitting: while submitting; cannot be updated
# running: still running, can be updated
# finished: after completely successful run
# error: some jobs had errors
proc job_update {logfile {cleanup success} {force 0}} {
	global cgjob
	if {![file exists $logfile]} {
		puts stderr "logfile $logfile not found, checking for finished logfile"
		set root [file root $logfile]
		foreach ext {finished error} {
			if {[file exists $root.$ext]} {
				puts stderr "run already finished with status $ext: $root.$ext"
				return
			}
		}
		if {[file exists $root.submitting]} {
			error "cannot update, run still submitting: $logfile"
		}
		error "logfile $logfile not found"
	}
	if {[file extension $logfile] eq ".submitting"} {
		error "cannot update logfile while still submitting: $logfile"
	}
	catch {close $f} ; catch {close $o}
	set tempfile [filetemp $logfile]
	set f [open $tempfile]
	hardcopy $logfile $tempfile
	# init needs to be here, because otherwise variables like cgpub(pid) will 
	# be reset on first use of functions in jobs.tcl
	set tempresult [filetemp $logfile]
	set o [open $tempresult w]
	while 1 {
		if {[gets $f line] == -1} break
		if {[string index $line 0] ne {#}} break
		if {[regexp {^# ([^:]+): (.*)$} $line temp key value]} {
			set cgjob($key) $value
		}
		puts $o $line
	}
	if {[isint $cgjob(distribute)]} {
		if {$cgjob(distribute) <= 1} {
			set target direct
		} else {
			set target distr
		}
	} else {
		set target $cgjob(distribute)
	}
	if {[split $line \t] ne {job jobid status submittime starttime endtime duration targets msg run}} {
		close $o ; close $f
		error "file $logfile is not a proper logfile (must be tsv with fields: job jobid status submittime starttime endtime duration targets msg run)"
	}
	puts $o $line
	set endstartcode {}
	set endstarttime {}
	set endendcode {}
	set endendtime {}
	set endstatus finished
	set totalduration {0 0}
	while 1 {
		if {[gets $f line] == -1} break
		set sline [split $line \t]
		foreach {jobo jobid status submittime starttime endtime duration targets msg run} $sline break
		set startcode [timescan $starttime]
		set endcode [timescan $endtime]
		if {$jobo eq "total"} {
			puts $o [join [list total $jobid $endstatus $submittime $endstarttime $endendtime [timediff2duration $totalduration] $targets [job_cleanmsg $msg] $run] \t]
			break
		}
		if {[get cgjob(basedir) ""] ne "" && [file pathtype $jobo] ne "absolute"} {
			set job $cgjob(basedir)/$jobo
		} else {
			set job $jobo
		}
		if {$status in {submitted running}} {set starttime {} ; set endtime {} ; set duration {}}
		if {($starttime eq "" || $endtime eq "" | $duration eq "" | $force) && [job_file_exists $job.log]} {
			set jobloginfo [job_parse_log $job $totalduration]
			foreach {status starttime endtime run duration totalduration} $jobloginfo break
			if {$status eq "error"} {
				if {[catch {set msg [file_read $job.err]}]} {set msg ""}
			} elseif {$status in {submitted running}} {
				# still in queue, running or hang/error?
				set status [job_status_$target $job $jobloginfo]
				if {[string range $duration end-2 end] eq "..."} {set endtime "" ; set duration ""}
				if {$status eq "error"} {
					if {[catch {set msg [file_read $job.err]}]} {
						set msg "job no longer running, but no error message found"
					}
				}
			}
		}
		if {$status in "running submitted"} {
			set endstatus running
		} elseif {$status eq "error" && $endstatus ne "running"} {
			set endstatus error
		}
		if {$startcode ne "" && $endcode ne ""} {
			set diff [lmath_calc $endcode - $startcode]
			set totalduration [lmath_calc $totalduration + $diff]
		}
		puts $o [join [list $jobo $jobid $status $submittime $starttime $endtime $duration $targets [job_cleanmsg $msg] $run] \t]
		if {$endstartcode eq "" || [time_comp $startcode $endstartcode] > 0} {set endstartcode $startcode ; set endstarttime $starttime}
		if {$endendcode eq "" || ($endcode ne "" && [time_comp $endendcode $endcode] > 0)} {set endendcode $endcode ; set endendtime $endtime}
	}
	close $o
	close $f
	# file rename -force $logfile $logfile.old
	if {$endstatus eq "running"} {
		set result [file root $logfile].running
	} elseif {$endstatus eq "error"} {
		set result [file root $logfile].error
	} else {
		set result [file root $logfile].finished
	}
	file delete -force $tempfile
	file rename -force $tempresult $result
	if {$cleanup eq "allways" || ($cleanup eq "success" && $endstatus eq "finished")} {
		job_cleanlogs $result
	}
	if {$logfile ne $result} {file delete $logfile}
}

proc cg_job_update args {
	job_init
	set cleanup success
	set force 0
	cg_options job_update args {
		-dcleanup - -cleanup - -c {
			if {$value ni {success never allways}} {error "$value not a valid option for -cleanup, should be one of: success, never, allways"}
			set cleanup $value
		}
		-f - -force {
			set force $value
		}
	}
	foreach logfile [lsort -dict -decreasing $args] {
		if {[file extension $logfile] in ".running .finished .error .submitting"} break
	}
	set logfile [file_absolute $logfile]
	job_update $logfile $cleanup $force
}
