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
	set header [split $line \t]
	if {[list_remove $header time_seconds cores] ne {job jobid status submittime starttime endtime duration targets msg run}} {
		close $f
		error "file $logfile is not a proper logfile (must be tsv with fields: job jobid status submittime starttime endtime duration time_seconds targets msg run cores)"
	}
	unset -nocomplain donea
	while 1 {
		if {[gets $f line] == -1} break
		set sline [split $line \t]
		foreach {jobo jobid status submittime starttime endtime duration time_seconds targets msg run cores} $sline break
		if {$jobo eq "total"} continue
		if {[get cgjob(basedir) ""] ne "" && [file pathtype $jobo] ne "absolute"} {
			set job $cgjob(basedir)/$jobo
		} else {
			set job $jobo
		}
		if {[info exists donea($job)]} continue
		foreach ext {log run out err ok} {
			file delete -force [job.file $ext $job]
		}
		catch {file delete [job.file finished $job]}
		set dirsa([file dir $job]) 1
		set donea($job) 1
	}
	foreach dir [lsort -decreasing [array names dirsa]] {
		job_delete_ifempty $dir 1
	}
	putslog "cleanup finished"
}

# correct error in values from older version of timestamp
proc correct_time_ms {timeVar} {
	upvar $timeVar time
	if {![regexp {^(.*\.)([0-9]?[0-9]?)$} $time temp pre ms]} {
		return 0
	}
	set time $pre[format %03d $ms]
	return 1
}

# status can be
# submitting: while submitting; cannot be updated
# running: still running, can be updated
# finished: after completely successful run
# error: some jobs had errors
# set cleanup never ; set force 0 ; set removeold 0 ; set rundone 0
proc job_update {logfile {cleanup success} {force 0} {removeold 0} {rundone 0}} {
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
	unset -nocomplain joblogcachea
	# get data from old log files
	set logroot [file root [file root $logfile]]
	set oldlogfiles [bsort [glob -nocomplain $logroot.*.finished $logroot.*.running $logroot.*.error $logroot.*.submitting]]
	set pos [lsearch $oldlogfiles $logfile]
	incr pos -1
	set oldlogfiles [lrange $oldlogfiles 0 $pos]
	unset -nocomplain oldlogsa
	set expectedheader "job jobid status submittime starttime endtime duration time_seconds targets msg run cores"
	foreach oldlogfile $oldlogfiles {
		set f [gzopen $oldlogfile]
		set header [tsv_open $f]
		set addseconds 0
		set addcores 0
		if {$header eq "job jobid status submittime starttime endtime duration targets msg run"} {
			set addseconds 1
		} elseif {$header eq "job jobid status submittime starttime endtime duration time_seconds targets msg run"} {
			set addcores 1
		} elseif {$header ne $expectedheader} {
			error "error in format of logfile $oldlogfile"
		}
		while {[gets $f line] != -1} {
			set line [split $line \t]
			foreach {job jobid status submittime starttime endtime} $line break
			if {$addseconds} {
				set time_seconds [timebetween_inseconds $starttime $endtime]
				set line [linsert $line 7 $time_seconds]
			}
			if {$addcores} {
				lappend line {}
			}
			if {![info exists oldlogsa($job)] || $status ne "skipped" || [lindex $oldlogsa($job) 2] ne "finished"} {
				set oldlogsa($job) $line
			}
		}
		gzclose $f
	}
	# start
	catch {close $f} ; catch {close $o}
	set tempfile [filetemp $logfile]
	hardcopy $logfile $tempfile
	set tempresult [filetemp $logfile]
	set f [open $tempfile]
	set o [open $tempresult w]
	# copy comments
	while 1 {
		if {[gets $f line] == -1} break
		if {[string index $line 0] ne {#}} break
		if {[regexp {^# ([^:]+): (.*)$} $line temp key value]} {
			set cgjob($key) $value
		}
		puts $o $line
	}
	if {$rundone} {unset -nocomplain cgjob(pid)}
	if {[isint $cgjob(distribute)]} {
		if {$cgjob(distribute) <= 1} {
			set target direct
		} else {
			set target distr
		}
	} else {
		set target $cgjob(distribute)
	}
	set header [split $line \t]
	set addseconds 0
	if {$header eq "job jobid status submittime starttime endtime duration targets msg run"} {
		set addseconds 1
	} elseif {$header eq "job jobid status submittime starttime endtime duration time_seconds targets msg run"} {
		set addcores 1
	} elseif {$header ne $expectedheader} {
		close $o ; close $f
		error "error in format of logfile $oldlogfile (must be tsv with fields: job jobid status submittime starttime endtime duration time_seconds targets msg run)"
	}
	puts $o [join $expectedheader \t]
	set totalstartcode {}
	set totalstarttime {}
	set totalendcode {}
	set totalendtime {}
	set endstatus finished
	set totalseconds 0
	while 1 {
		if {[gets $f line] == -1} break
		set sline [split $line \t]
		if {$addseconds} {
			foreach {jobo jobid status submittime starttime endtime duration targets msg run} $sline break
			set cs [correct_time_ms starttime] ; set ce [correct_time_ms endtime]
			set startcode [timescan starttime "Could not interpret starttime $starttime correctly in $logfile, line $line"]
			set endcode [timescan endtime "Could not interpret endtime $endtime correctly in $logfile, line $line"]
			if {$cs || $ce} {
				set duration [timediff2duration [lmath_calc $endcode - $startcode]]
			}
			set time_seconds [timebetween_inseconds $starttime $endtime]
		} else {
			set cores {}
			foreach {jobo jobid status submittime starttime endtime duration time_seconds targets msg run cores} $sline break
			set cs [correct_time_ms starttime] ; set ce [correct_time_ms endtime]
			set startcode [timescan starttime "Could not interpret starttime $starttime correctly in $logfile, line $line"]
			set endcode [timescan endtime "Could not interpret endtime $endtime correctly in $logfile, line $line"]
			if {$cs || $ce} {
				set duration [timediff2duration [lmath_calc $endcode - $startcode]]
				set time_seconds [timebetween_inseconds $starttime $endtime]
			}
		}
		if {$jobo eq "total"} {
			set walltime [timediff2duration [lmath_calc $totalendcode - $totalstartcode]]
			puts $o [join [list total $jobid $endstatus $submittime $totalstarttime $totalendtime "$walltime (wall) [time_seconds2duration $totalseconds] (total)" "$totalseconds (total)" $targets [job_cleanmsg $msg] $run {}] \t]
			break
		}
		if {[get cgjob(basedir) ""] ne "" && [file pathtype $jobo] ne "absolute"} {
			set job $cgjob(basedir)/$jobo
		} else {
			set job $jobo
		}
		if {[file exists [job.file jid $job]] && [job_running [file_read [job.file jid $job]]]} {
			set status submitted
		}
		if {$status in {submitted running}} {set endtime {} ; set duration {}; set time_seconds {}}
		if {$status eq "skipped" && [info exists oldlogsa($jobo)]} {
			foreach {jobo jobid status submittime starttime endtime duration time_seconds targets msg run} $oldlogsa($jobo) break
			if {$status in "error skipped" && [job_file_or_link_exists [job.file log $job]]} {
				set jobloginfo [job_parse_log $job]
				foreach {status starttime endtime run duration submittime time_seconds} $jobloginfo break
			}
			correct_time_ms starttime
			correct_time_ms endtime
			set startcode [timescan starttime "Could not interpret starttime $starttime correctly in $logfile, line $line"]
			set endcode [timescan endtime "Could not interpret endtime $endtime correctly in $logfile, line $line"]
			set duration [timediff2duration [lmath_calc $endcode - $startcode]]
			set time_seconds [timebetween_inseconds $starttime $endtime]
		}
		if {$starttime eq "" || $endtime eq "" | $duration eq "" | $force | ($status eq "error" && [job_file_or_link_exists [job.file log $job]])} {
			if {[job_file_or_link_exists [job.file log $job]]} {
				set jobloginfo [job_parse_log $job]
				foreach {status starttime endtime run duration submittime time_seconds} $jobloginfo break
				if {($starttime eq "" || $endtime eq "") && [info exists joblogcachea($job)]} {
					set jobloginfo $joblogcachea($job)
					foreach {status starttime endtime run duration submittime time_seconds} $jobloginfo break
				}
				if {$status eq "skipped"} {
					set joblogcachea($job) $jobloginfo
				}
				set cs [correct_time_ms starttime] ; set ce [correct_time_ms endtime]
				set startcode [timescan starttime "Could not interpret starttime $starttime correctly in $logfile, line $line"]
				set endcode [timescan endtime "Could not interpret endtime $endtime correctly in $logfile, line $line"]
				if {$cs || $ce} {
					set duration [timediff2duration [lmath_calc $endcode - $startcode]]
					set time_seconds [timebetween_inseconds $starttime $endtime]
				}
				set msg {}
			} else {
				set jobloginfo [list $status $starttime $endtime $run $duration $submittime $time_seconds]
			}
			if {$status eq "error"} {
				if {[catch {set msg [file_read [job.file err $job]]}]} {set msg ""}
			} elseif {$status in {submitted running}} {
				# still in queue, running or hang/error?
				set status [job_status_$target $job $jobloginfo]
				if {[string range $duration end-2 end] eq "..."} {set endtime "" ; set duration ""}
				if {$status eq "error"} {
					if {[catch {set msg [file_read [job.file err $job]]}]} {
						set msg "job $job no longer running, but no error message found"
					}
				}
			}
		} else {
			set duration [timediff2duration [lmath_calc $endcode - $startcode]]
		}
		if {$status in "running submitted"} {
			set endstatus running
		} elseif {$status eq "error" && $endstatus ne "running"} {
			set endstatus error
		}
		if {$startcode ne "" && $endcode ne "" && $time_seconds ne ""} {
			if {![isint $cores]} {set cores 1}
			set totalseconds [expr {$totalseconds + $cores * $time_seconds}]
		}
		puts $o [join [list $jobo $jobid $status $submittime $starttime $endtime $duration $time_seconds $targets [job_cleanmsg $msg] $run $cores] \t]
		if {$totalstartcode eq "" || [time_comp $startcode $totalstartcode] < 0} {
			set totalstartcode $startcode ; set totalstarttime $starttime
		}
		if {$totalendcode eq "" || ($endcode ne "" && [time_comp $endcode $totalendcode] > 0)} {
			set totalendcode $endcode ; set totalendtime $endtime
		}
	}
	close $o
	close $f
	# file rename -force -- $logfile $logfile.old
	if {$endstatus eq "running"} {
		set result [file root $logfile].running
	} elseif {$endstatus eq "error"} {
		set result [file root $logfile].error
	} else {
		set result [file root $logfile].finished
	}
	file delete -force $tempfile
	file rename -force -- $tempresult $result
	if {$cleanup eq "allways" || ($cleanup eq "success" && $endstatus eq "finished")} {
		job_cleanlogs $result
	}
	if {$logfile ne $result} {file delete $logfile}
	if {$removeold} {
		foreach oldlogfile $oldlogfiles {
			if {$cleanup eq "allways" || ($cleanup eq "success" && $endstatus eq "finished")} {
				job_cleanlogs $oldlogfile
			}
			file delete $oldlogfile
		}
	}
}

proc cg_job_update args {
	# init needs to be here, because otherwise variables like cgpub(pid) will 
	# be reset on first use of functions in jobs.tcl
	set cleanup success
	set force 0
	set removeold 0
	cg_options job_update args {
		-dcleanup - -cleanup - -c {
			if {$value ni {success never allways}} {error "$value not a valid option for -cleanup, should be one of: success, never, allways"}
			set cleanup $value
		}
		-dremoveold - -removeold - -r {
			set removeold $value
		}
		-f - -force {
			set force $value
		}
	}
	foreach logfile [bsort -decreasing $args] {
		if {[file extension $logfile] in ".running .finished .error .submitting"} break
	}
	set logfile [file_absolute $logfile]
	set distribute sge
	set f [open $logfile]
	while 1 {
		if {[gets $f line] == -1} break
		if {[string index $line 0] ne "\#"} break
		if {[regexp {# *distribute: *(.+)} $line temp distribute]} break
	}
	job_init -d $distribute
	job_update $logfile $cleanup $force $removeold
}
