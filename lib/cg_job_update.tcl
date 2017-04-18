proc job_update {logfile} {
	global cgjob
	if {![file exists $logfile]} {
		error "cannot update non(no longer)-existing logfile: $logfile"
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
		set target $type
	}
	if {[split $line \t] ne {job jobid status submittime starttime endtime duration targets msg}} {
		close $o ; close $f
		error "file $logfile is not a proper logfile (must be tsv with fields: job jobid status submittime starttime endtime duration targets msg)"
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
		foreach {job jobid status submittime starttime endtime duration targets msg} $sline break
		set startcode [timescan $starttime]
		set endcode [timescan $endtime]
		if {$job eq "total"} {
			puts $o [join [list $job $jobid $endstatus $submittime $endstarttime $endendtime [timediff2duration $totalduration] $targets [job_cleanmsg [job_cleanmsg $msg]]] \t]
			break
		}
		switch $status {
			submitted - running {
				set duration {}
				set starttime {} ; set endtime {}
				if {![job_file_exists $job.log]} {
					puts $o $line
					continue
				}
				set jobloginfo [job_parse_log $job $totalduration]
				foreach {failed starttime endtime duration totalduration} $jobloginfo break
				if {$failed} {
					if {[catch {set errormsg [file_read $job.err]}]} {set errormsg ""}
					puts $o [join [list $job $jobid error $submittime $starttime $endtime $duration $targets $errormsg] \t]
					if {$endstatus ne "running"} {set endstatus error}
				} elseif {$endtime ne ""} {
					puts $o [join [list $job $jobid finished $submittime $starttime $endtime $duration $targets ""] \t]
				} else {
					# still in queue, running or hang/error?
					set status [job_status_$target $job $jobloginfo]
					set msg ""
					if {$status in "running submitted"} {
						set endstatus running
					} elseif {$status eq "error"} {
						if {$endstatus ne "running"} {set endstatus error}
						if {[catch {set msg [file_read $job.err]}]} {
							set msg "job no longer running, but no error message found"
						}
					}
					puts $o [join [list $job $jobid $status $submittime $starttime $endtime $duration $targets [job_cleanmsg $msg]] \t]
				}
			}
			error {
				if {$endstatus ne "running"} {set endstatus error}
				puts $o $line
			}
			default {
				if {$startcode ne "" && $endcode ne ""} {
					set diff [lmath_calc $endcode - $startcode]
					set totalduration [lmath_calc $totalduration + $diff]
				}
				puts $o $line
			}
		}
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
	if {$logfile ne $result} {file delete $logfile}
}

proc cg_job_update args {
	job_init
	cg_options job_update args {
	} logfile 1 1
	set logfile [file_absolute $logfile]
#	if {[regexp {^(.*)\.[0-9_-]+\.([a-z]+)$} $logfile temp base filestatus]} {
#	} elseif {[regexp {^(.*)\.[0-9_-]+$} $logfile temp base]} {
#	} else {
#		set base $logfile
#	}
#	set logfile [lindex [lsort -dict [glob -nocomplain $base.*.running]] end]
	job_update $logfile
}
