proc job_process_sge_init {} {
	# we will use par (parallel) code with some specifics for sge
	if {[info commands job_process_par] eq ""} {auto_load job_process_par}
	interp alias {} job_process {} job_process_par
	interp alias {} job_running {} job_running_sge
	interp alias {} job_wait {} job_process_sge_wait
	interp alias {} job_process_par_submit {} job_process_sge_submit
}

proc job_running_sge {jobid} {
	expr {[isint $jobid] && ![catch {exec qstat -j $jobid}]} 
}

proc job_process_sge_submit {job runfile args} {
	set options {}
	set soft {}
	set hard {}
	set pos 0
	foreach {opt value} $args {
		switch -- $opt {
			-deps {
				set value [list_remove $value {}]
				if {[llength $value]} {
					lappend options -hold_jid [join $value ,]
				}
				incr pos 2
			}
			-cores {
				if {![info exists cgjob_distr(no_local_pe)]} {
					set cgjob_distr(no_local_pe) [catch {exec qconf -sp local}]
				}
				if {!$cgjob_distr(no_local_pe)} {
					lappend hard -pe local $value
				}
				incr pos 2
			}
			-hard {
				lappend hard -l $value
				incr pos 2
			}
			-soft {
				lappend soft -l $value
				incr pos 2
			}
			-host {
				lappend hard -l hostname=$value
				incr pos 2
			}
			-io {
				lappend soft -l io=$value
				incr pos 2
			}
			-- {
				incr pos 1
				break
			}
			default {
				break
			}
		}
	}
	if {[llength $soft]} {
		lappend options -soft {*}$soft
	}
	if {[llength $hard]} {
		lappend options -hard {*}$hard
	}
	set name "[file tail $job] $job"
	regsub -all / $name __ name
	if {[string length $name] > 200} {
		set name [string range $name 0 100]....[string range $name end-100 end]
	}
	set dir [file dir $job]
	catch {file delete $job.finished}
	catch {file delete $job.out}
	catch {file delete $job.err}
	set jnum [exec qsub -N j$name -q all.q -o $job.out -e $job.err {*}$options $runfile]
	regexp {[0-9]+} $jnum jobnum
	return $jobnum
}

if 0 {
# debug proc
	set curjobnum 1
	proc job_process_sge_submit {job runfile args} {
		global curjobnum
		return [incr curjobnum]
	}
	interp alias {} job_process_par_submit {} job_process_sge_submit
}

proc job_process_sge_wait {} {
}

