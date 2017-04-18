proc job_process_direct_init {} {
	set ::job_method_info {}
	interp alias {} job_process {} job_process_direct
	interp alias {} job_running {} job_running_direct
	interp alias {} job_wait {} job_process_direct_wait
}

# not actually used for direct jobs
proc job_running_direct {job} {
	return 0
}

# not actually used for direct jobs
proc job_status_direct {job {jobloginfo {}}} {
	if {$jobloginfo eq ""} {
		if {![file exists $job.log]} {return unkown}
		set jobloginfo [job_parse_log $job $totalduration]
	}
	foreach {failed starttime endtime duration totalduration} $jobloginfo break
	if {$failed} {
		return error
	} elseif {$endtime ne ""} {
		return finished
	} else {
		return error
	}
}

proc stderr2file {fileout {fileerr {}}} {
	if {$fileout ne ""} {
		set ::stderr_redirect $fileerr
		if {![llength [info command __puts] ]} {
			rename puts __puts
		}
		proc puts {args} [string_change {
			if {[llength $args] > 1} {
				set chan [lindex $args 0]
				if {$chan eq "stderr"} {
					set f [open {@FILEERR@} a+]
					catch {__puts {*}[lreplace $args 0 0 $f]}
					close $f
				} elseif {$chan eq "stdout"} {
					set f [open {@FILEOUT@} a+]
					catch {__puts {*}[lreplace $args 0 0 $f]}
					close $f
				}
			} else {
				set f [open {@FILEOUT@} a+]
				catch {__puts $f [lindex $args 0]}
				close $f
			}
			if {[catch {
				__puts {*}$args
			} err]} {
				return -code error $err
			}
		} [list @FILEOUT@ [file_absolute $fileout] @FILEERR@ [file_absolute $fileerr]]]
	} else {
		unset -nocomplain ::stderr_redirect
		if {[llength [info command __puts] ]} {
			rename puts {}
			rename __puts puts
		}
	}
}

proc job_process_direct {} {
	global cgjob job_deps
	set jobroot [pwd]
	while {[llength $cgjob(queue)]} {
		set joberror {}
		set line [list_shift cgjob(queue)]
		foreach {jobid jobname job_logdir pwd deps foreach ftargetvars ftargets fptargets fskip checkcompressed code submitopts frmtargets precode jobforce optional} $line break
		cd $pwd
		set timefile {}
		set job [job_logname $job_logdir $jobname]
		# check foreach deps, skip if not fullfilled
		# add all resulting (foreach) jobs in front of the queue, and go back to running the queue
		set submittime [timestamp]
		if {[llength $foreach]} {
			if {[catch {job_finddeps $job $foreach ftargetvars 1 fids time timefile $checkcompressed} fadeps]} {
				if {![regexp {^missing dependency} $fadeps]} {
					set errormsg "error in foreach dependencies for $jobname: $fadeps"
					job_log $job $errormsg
					job_logfile_add $job . error $ftargets $errormsg $submittime
				} else {
					job_log $job "$fadeps"
					job_logfile_add $job . error $ftargets $fadeps $submittime
				}
				job_log $job "-----> job $jobname failed"
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
			set cgjob(queue) [list_concat $temp $cgjob(queue)]
			job_logclose $job
			continue
		}
		cd $pwd
		job_log $job "==================== $jobname ===================="
		# check deps, skip if not fullfilled
		if {[catch {job_finddeps $job $deps newtargetvars 0 ids time timefile $checkcompressed $ftargetvars} adeps]} {
			if {![regexp {^missing dependency} $adeps]} {
				job_log $job "error in dependencies for $jobname: $adeps"
			} else {
				job_log $job "$adeps"
			}
			if {$optional || $cgjob(skipjoberrors)} {
				job_log $job "-----> job $jobname skipped: dependencies not found"
				job_logfile_add $job . skipped $ftargets "dependencies not found" $submittime
				job_logclose $job
				continue
			} else {
				set joberror "error trying to run job $jobname:\n$adeps"
			}
		}
		set targetvars $ftargetvars
		lappend targetvars {*}$newtargetvars
		if {$cgjob(force)} {set time force}
		# check skip targets, if already done or running, skip
		if {$jobforce} {job_log $job "forcing $jobname"}
		if {!$cgjob(force) && !$jobforce && [llength $fskip]} {
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
				job_logfile_add $job . skipped $ftargets "skip targets ($skip) already completed or running" $submittime
				continue
			}
		}
		# check targets, if already done or running, skip
		# =============
		set run 0
		if {$ftargets ne ""} {
			set targets [job_targetsreplace $ftargets $targetvars]
			if {!$jobforce && ![job_checktargets $job $targets $time $timefile $checkcompressed running]} {
				set run 1
			}
		} else {
			set targets {}
			set running {}
			set run 1
		}
		set ptargets [job_targetsreplace $fptargets $targetvars]
		if {$jobforce || ([llength $ptargets] && ![llength [job_findptargets $ptargets $checkcompressed]])} {
			set run 1
		}
		if {$cgjob(force)} {
			if {[llength $running]} {
				error "cannot force job with still running tasks ([join $running ,])"
			}
			foreach target $targets {
				# job_backup $target 1
				file delete $target
			}
			set run 1
		}
		if {!$run} {
			job_log $job "skipping $jobname: targets already completed or running"
			job_logfile_add $job . skipped $targets "targets already completed or running" $submittime
			job_logclose $job
			continue
		}
		if {$joberror ne ""} {
			job_logfile_add $job . error $ftargets $joberror $submittime
			error $joberror
		}
		job_log $job "-------------------- running $jobname --------------------"
		set starttime [timestamp]
		# run code
		# ========
		set cmd "proc job_run {} \{\n"
		append cmd [job_generate_code $job $pwd $adeps $targetvars $targets $ptargets $checkcompressed $code]
		append cmd \}
		set ok 1
		if {[catch {eval $cmd} result]} {
			set ok 0
			job_log $job "error creating $jobname: $result"
		}
		catch {file delete $job.finished}
		set f [open $job.err w]; close $f
		stderr2file $job.out $job.err
		set error [catch {job_run} result]
		stderr2file {}
		if {$error} {
			set errormessage $result\n$::errorInfo
			file_add $job.err $errormessage
			if ($cgjob(skipjoberrors)) {
				puts stderr $result
			} else {
				error $result $::errorInfo
			}
		}
		set endtime [timestamp]
		# log results
		# ===========
		if {![job_file_exists $job.finished]} {
			job_log $job "$jobname failed: did not finish\nerror:\n$result\n"
			job_logfile_add $job . error $targets "did not finish\nerror:\n$result" $submittime $starttime $endtime
		} elseif {$error} {
			file delete $job.finished
			job_logfile_add $job . error $targets $errormessage $submittime $starttime $endtime
		} else {
			job_log $job "-------------------- end $jobname --------------------"
			job_log $job "$jobname finished successfully\n"
			job_logfile_add $job . finished $targets "" $submittime $starttime $endtime
		}
		job_logclose $job
	}
	cd $jobroot
}

proc job_logfile_direct_close {} {
	global cgjob
	if {![info exists cgjob(f_logfile)]} return
	puts $cgjob(f_logfile) [join [list total . finished $cgjob(starttime) "" $cgjob(endtime) [timediff2duration $cgjob(totalduration)] "" ""] \t]
	close $cgjob(f_logfile)
	if {$cgjob(status) eq "error"} {
		file rename $cgjob(logfile).submitting $cgjob(logfile).error
	} else {
		file rename $cgjob(logfile).submitting $cgjob(logfile).finished
	}
}

proc job_process_direct_wait {} {
	job_logfile_direct_close
}
