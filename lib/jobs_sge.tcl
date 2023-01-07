proc job_process_init_sge {} {
	# we will use par (parallel) code with some specifics for sge
	if {[info commands job_process_par] eq ""} {auto_load job_process_par}
	set ::cgjob(endjobids) {}
	set ::job_method_info {}
	interp alias {} job_process {} job_process_par
	interp alias {} job_runall {} job_runall_par
	interp alias {} job_running {} job_running_sge
	interp alias {} job_wait {} job_wait_sge
	interp alias {} job_process_submit_par {} job_process_submit_sge
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
	foreach {status starttime endtime run duration submittime} $jobloginfo break
	if {$status ni {submitted running}} {return $status}
	set jobnum [job_process_par_jobid $job]
	if {[job_running_sge $jobnum]} {
		return $status
	} else {
		return error
	}
}

proc job_process_submit_sge {job runfile args} {
	global cgjob
	set options {}
	set soft {}
	set hard {}
	set priority [get cgjob(priority) 0]
	set cores 1
	set mem {}
	set time {}
	set pos 0
	set dqueue [get cgjob(dqueue) all.q]
	foreach {opt value} $args {
		switch -- $opt {
			-deps {
				set value [list_remove $value {} q]
				if {[llength $value]} {
					lappend options -hold_jid [join $value ,]
					set cgjob(endjobids) [list_lremove $cgjob(endjobids) $value]
				}
				incr pos 2
			}
			-priority {
				set priority $value
			}
			-cores {
				# we will use a PE named local, which must be present/made on the cluster
				# if PE local does not exists, PE smp is tried
				# -l slots=$value mentioned in http://www.softpanorama.org/HPC/Grid_engine/Reference/qsub.shtml
				# would be easier, but does not seem to work (if there is any PE defined?)
				# lappend hard -l slots=$value
				set cores $value
				if {![info exists cgjob_distr(local_pe)]} {
					set error [catch {exec qconf -sp local}]
					if {!$error} {
						set cgjob_distr(local_pe) local
					} else {
						set error [catch {exec qconf -sp smp}]
						if {!$error} {
							set cgjob_distr(local_pe) smp
						} else {
							set cgjob_distr(local_pe) {}
						}
					}
				}
				if {$cgjob_distr(local_pe) ne ""} {
					lappend hard -pe $cgjob_distr(local_pe) $cores -R y
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
				set mem $value
				incr pos 2
			}
			-time {
				set time $value
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
	if {[llength $cgjob(dmem)] == 1} {set cgjob(dmem) [list * $cgjob(dmem)]}
	foreach {pattern value} $cgjob(dmem) {
		if {$pattern eq "*" || [regexp $pattern $job]} {
			if {[job_memgt $value $mem]} {
				set mem $value
			}
		}
	}
	if {[llength $cgjob(dtime)] == 1} {set cgjob(dtime) [list * $cgjob(dtime)]}
	foreach {pattern value} $cgjob(dtime) {
		if {$pattern eq "*" || [regexp $pattern $job]} {
			if {[bsort [list $time $value]] eq [list $time $value]} {
				set time $value
			}
		}
	}
	if {$mem ne ""} {
		# mem_free: only start job if the given amount of memory is free on the node
		#  -> often not sufficient if running jobs increase memory use during run
		# virtual_free: reserves this much memory, can be defined as a consumable resource
		#	-> if all memory is reserved for running kjobs, no new jobs are started on node
		# not using h_vmem, because that would kill any job going (even a bit) above reserved memory
		set temp [job_mempercore $mem $cores]
		lappend hard -l mem_free=$temp,virtual_free=$temp
	}
	if {$time ne ""} {
		# (time format is hh:mm:ss)
		lappend soft -l s_rt=$time
	}
	if {[llength $soft]} {
		lappend options -soft {*}$soft
	}
	if {[llength $hard]} {
		lappend options -hard {*}$hard
	}
	set name "[file tail $job] $job"
	set prefix [file tail [get cgjob(prefix) "j"]]
	set name [sge_safename $name $prefix]
	set dir [file dir $job]
	catch {file delete $job.finished}
	catch {file delete $job.ok}
	catch {file delete $job.out}
	catch {file delete $job.err}
	if {[regexp , $job]} {
		error "Cannot submit job to sge: it has a comma in the output file $job.out, which grid engine sometimes has problems with"
	}
	putslog "sge_submit: [list qsub -N j$name -q $dqueue -o $job.out -e $job.err -p $priority {*}$options $runfile]"
	set jnum [exec qsub -N j$name -q $dqueue -o $job.out -e $job.err -p $priority {*}$options $runfile]
	regexp {[0-9]+} $jnum jobnum
	lappend cgjob(endjobids) $jobnum
	return $jobnum
}

proc sge_safename {name {prefix {}}} {
	regsub -all {[^A-Za-z0-9_.-]} $name __ name
	if {$prefix ne ""} {
		regsub -all {[^A-Za-z0-9_.-]} $prefix __ prefix
		set plen [string length $prefix]
		if {$plen > 99} {
			set prefix [string range $prefix 0 47]..[string range $prefix end-48 end]#
		} else {
			set prefix ${prefix}#
		}
	}
	set name ${prefix}$name
	if {[string length $name] > 200} {
		set name [string range $name 0 100]....[string range $name end-100 end]
	}
	return $name
}

if 0 {
# debug proc
	set curjobnum 1
	proc job_process_submit_sge {job runfile args} {
		global curjobnum
		return [incr curjobnum]
	}
	interp alias {} job_process_submit_par {} job_process_submit_sge
}

proc job_logfile_sge_close {} {
	global cgjob
	set logfile $cgjob(logfile).running
	set runfile $logfile.update
	set outfile $logfile.update.out
	set errfile $logfile.update.err
	job_update $logfile $cgjob(cleanup) 1
	set statusok [file exists $cgjob(logfile).finished]
	if {$cgjob(cleanup) eq "allways" || ($cgjob(cleanup) eq "success" && $statusok)} {
		set result [glob $cgjob(logfile).finished $cgjob(logfile).running $cgjob(logfile).error]
		job_cleanup
		job_cleanlogs $result
		# only keep result logfile if -d option was given explicitely
		if {!$cgjob(hasargs)} {file delete $result}
	}
	file delete -force $runfile
	file delete -force $outfile
	file delete -force $errfile
}

proc job_wait_sge {} {
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
	append cmd "job_init -d sge\n"
	append cmd [list array set cgjob [array get cgjob]]\n
	append cmd job_logfile_sge_close\n
	file_write $runfile $cmd
	file attributes $runfile -permissions u+x
	set options {}
	if {[llength $cgjob(endjobids)]} {
		lappend options -hold_jid [join $cgjob(endjobids) ,]
	}
	exec qsub -N j$name -q [get cgjob(dqueue) all.q] -o $outfile -e $errfile -p $priority {*}$options $runfile
}
