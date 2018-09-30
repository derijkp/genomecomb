proc job_process_sge_init {} {
	# we will use par (parallel) code with some specifics for sge
	if {[info commands job_process_par] eq ""} {auto_load job_process_par}
	set ::cgjob(alljobids) {}
	set ::job_method_info {}
	interp alias {} job_process {} job_process_par
	interp alias {} job_running {} job_running_sge
	interp alias {} job_wait {} job_process_sge_wait
	interp alias {} job_process_par_submit {} job_process_sge_submit
}

proc job_running_sge {jobid} {
	expr {[isint $jobid] && ![catch {exec qstat -j $jobid}]} 
}

proc job_status_sge {job {jobloginfo {}}} {
	global cgjob_distr_running
	if {$jobloginfo eq ""} {
		if {![file exists $job.log]} {return unkown}
		set jobloginfo [job_parse_log $job]
	}
	foreach {status starttime endtime run duration totalduration submittime} $jobloginfo break
	if {$status ni {submitted running}} {return $status}
	set jobnum [job_process_par_jobid $job]
	if {[job_running_sge $jobnum]} {
		return $status
	} else {
		return error
	}
}

proc job_process_sge_submit {job runfile args} {
	global cgjob
	set options {}
	set soft {}
	set hard {}
	set priority [get cgjob(priority) 0]
	set pos 0
	set dqueue [get cgjob(dqueue) all.q]
	foreach {opt value} $args {
		switch -- $opt {
			-deps {
				set value [list_remove $value {}]
				if {[llength $value]} {
					lappend options -hold_jid [join $value ,]
				}
				incr pos 2
			}
			-priority {
				set priority $value
			}
			-cores {
				# we will use a PE named local, which must be present/made on the cluster
				# -l slots=$value mentioned in http://www.softpanorama.org/HPC/Grid_engine/Reference/qsub.shtml
				# would be easier, but does not seem to work (if there is any PE defined?)
				# lappend hard -l slots=$value
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
			-mem {
				lappend hard -l mem_free=$value
				incr pos 2
			}
			-dqueue {
				set dqueue $value
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
	set prefix [file tail [get cgjob(prefix) "j"]]
	regsub -all {[^A-Za-z0-9_.-]} $name __ name
	regsub -all {[^A-Za-z0-9_.-]} $prefix __ prefix
	set name ${prefix}#$name
	if {[string length $name] > 200} {
		set name [string range $name 0 100]....[string range $name end-100 end]
	}
	set dir [file dir $job]
	catch {file delete $job.finished}
	catch {file delete $job.out}
	catch {file delete $job.err}
	set jnum [exec qsub -N j$name -q $dqueue -o $job.out -e $job.err -p $priority {*}$options $runfile]
	regexp {[0-9]+} $jnum jobnum
	lappend cgjob(alljobids) $jobnum
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
	global cgjob cgjob_id
	set priority [get cgjob(priority) 0]
	job_logfile_par_close
	set logfile $cgjob(logfile).running
	set name [file tail $logfile]
	regsub -all {[^A-Za-z0-9_.-]} $name __ name
	set runfile $logfile.update
	set outfile $logfile.update.out
	set errfile $logfile.update.err
	set cmd {#!/bin/sh}
	append cmd \n
	append cmd {#$ -S /bin/bash} \n
	append cmd {#$ -V} \n
	append cmd {#$ -cwd} \n
	append cmd "\n\# the next line restarts using runcmd (specialised tclsh) \\\n"
	append cmd "exec $cgjob(runcmd) \"\$0\" \"\$@\"\n\n"
	append cmd job_init\n
	append cmd [list job_update $logfile $cgjob(cleanup)]\n
	append cmd [list file delete -force $runfile]\n
	append cmd [list file delete -force $outfile]\n
	append cmd [list file delete -force $errfile]\n
	file_write $runfile $cmd
	file attributes $runfile -permissions u+x
	set options {}
	if {[llength $cgjob(alljobids)]} {
		lappend options -hold_jid [join $cgjob(alljobids) ,]
	}
	exec qsub -N j$name -q [get cgjob(dqueue) all.q] -o $outfile -e $errfile -p $priority {*}$options $runfile
}
