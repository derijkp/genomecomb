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
		} [list @FILEOUT@ [file normalize $fileout] @FILEERR@ [file normalize $fileerr]]]
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
		set line [list_shift cgjob(queue)]
		foreach {jobid jobname job_logdir pwd deps foreach ftargetvars ftargets fptargets fskip code submitopts} $line break
		cd $pwd
		set job [job_logname $job_logdir $jobname]
		# check foreach deps, skip if not fullfilled
		# add all resulting (foreach) jobs in front of the queue, and go back to running the queue
		if {[llength $foreach]} {
			if {[catch {job_finddeps $job $foreach ftargetvars 1 fids time} fadeps]} {
				if {![regexp {^missing dependency} $fadeps]} {
					job_log $job "error in foreach dependencies for $jobname: $fadeps"
				} else {
					job_log $job "$fadeps"
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
		if {[catch {job_finddeps $job $deps newtargetvars 0 ids time $ftargetvars} adeps]} {
			if {![regexp {^missing dependency} $adeps]} {
				job_log $job "error in dependencies for $jobname: $adeps"
			} else {
				job_log $job "$adeps"
			}
			job_log $job "-----> job $jobname skipped: dependencies not found"
			job_logclose $job
			continue
		}
		set targetvars $ftargetvars
		lappend targetvars {*}$newtargetvars
		# check skip targets, if already done or running, skip
		if {!$cgjob(force) && [llength $fskip]} {
			set doskip 0
			foreach skip $fskip {
				set skip [job_targetsreplace $skip $targetvars]
				if {[llength $skip] && [job_checktargets $job $skip $time running]} {
					set doskip 1
					break
				}
			}
			if {$doskip} {
				job_log $job "skipping $jobname: skip targets already completed or running"
				continue
			}
		}
		# check targets, if already done or running, skip
		set run 0
		set targets [job_targetsreplace $ftargets $targetvars]
		if {![job_checktargets $job $targets $time running]} {
			set run 1
		}
		set ptargets [job_targetsreplace $fptargets $targetvars]
		if {[llength $ptargets] && ![llength [job_findptargets $ptargets]]} {
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
			job_logclose $job
			continue
		}
		job_log $job "-------------------- running $jobname --------------------"
		# run code
		set cmd "proc job_run {} \{\n"
		append cmd [job_generate_code $job $pwd $adeps $targetvars $targets $ptargets $code]
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
			puts stderr $result
			file_add $job.err $result
		}
		if {![file exists $job.finished]} {
			job_log $job "-----> job $jobname failed: did not finish\nerror:\n$result\n"
		} elseif {$error} {
			file delete $job.finished
			job_log $job "-----> job $job failed: error processing job\nerror:\n$result\n"
		} else {
			job_log $job "-------------------- end $jobname --------------------"
			job_log $job "-----> job $jobname finished successfully\n"
		}
		job_logclose $job
	}
	cd $jobroot
}

proc job_process_direct_wait {} {
}
