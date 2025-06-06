proc job_process_init_slurm {} {
	# we will use par (parallel) code with some specifics for slurm
	if {[info commands job_process_par] eq ""} {auto_load job_process_par}
	set ::cgjob(endjobids) {}
	set ::job_method_info {}
	interp alias {} job_process {} job_process_par
	interp alias {} job_runall {} job_runall_par
	interp alias {} job_running {} job_running_slurm
	interp alias {} job_wait {} job_wait_slurm
	interp alias {} job_process_submit_par {} job_process_submit_slurm
}

proc job_running_slurm {jobid} {
	if {![isint $jobid]} {return 0}
	set error [catch {exec squeue -o %t -j $jobid} msg]
	if {$error} {return 0}
	set state [lindex $msg end]
	if {$state in "CF CG PD R RD RF RH RQ RS RV SI SO"} {return 1} else {return 0}
}

proc job_status_slurm {job {jobloginfo {}}} {
	global cgjob_distr_running
	if {$jobloginfo eq ""} {
		if {![file exists [job.file log $job]]} {return unkown}
		set jobloginfo [job_parse_log $job]
	}
	foreach {status starttime endtime run duration submittime} $jobloginfo break
	if {$status ni {submitted running}} {return $status}
	set jobnum [job_process_par_jobid $job]
	if {[job_running_slurm $jobnum]} {
		return $status
	} else {
		return error
	}
}

proc job_process_submit_slurm {job cmd args} {
	global cgjob cgjob_id cgjob_running
	set options {}
	set soft {}
	set hard {}
	set priority [get cgjob(priority) 0]
	set cores 1
	set mem {}
	set time {}
	set pos 0
	set dqueue {}
	foreach {opt value} $args {
		switch -- $opt {
			-deps {
				set value [list_remove $value {} q]
				if {[llength $value]} {
					lappend options --dependency=afterany:[join $value :]
					set cgjob(endjobids) [list_lremove $cgjob(endjobids) $value]
				}
				incr pos 2
			}
			-priority {
				set priority $value
				lappend options --nice=$priority
			}
			-cores {
				set cores $value
				lappend options --ntasks=1 --cpus-per-task=$cores
				incr pos 2
			}
			-hard {
				# not used on slurm
				incr pos 2
			}
			-soft {
				# not used on slurm
				incr pos 2
			}
			-host {
				lappend options --nodelist=[join $host ,]
				incr pos 2
			}
			-io {
				# not used on slurm
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
		set temp [string tolower [job_mempercore $mem $cores]]
		lappend options --mem-per-cpu=$temp
	}
	if {$time ne ""} {
		lappend options --time=$time
	}
	set name "[file tail $job] $job"
	set prefix [file tail [get cgjob(prefix) "j"]]
	set name [slurm_safename $name $prefix]
	set dir [file dir $job]
	set job_out [job.file out $job]
	set job_err [job.file err $job]
	catch {file delete [job.file finished $job]}
	catch {file delete [job.file ok $job]}
	catch {file delete $job_out}
	catch {file delete $job_err}
	if {[regexp , $job]} {
		error "Cannot submit job to slurm: it has a comma in the output file $job_out"
	}
	# make runscript
	set runcmd {}
	append runcmd {#!/bin/sh}
	append runcmd \n
	append runcmd {#$ -S /bin/bash} \n
	append runcmd {#$ -V} \n
	append runcmd {#$ -cwd} \n
	append runcmd $cmd
	set runfile [job.file run $job]
	file_write $runfile $runcmd
	file attributes $runfile -permissions u+x
	if {$dqueue eq "" && [info exists cgjob(dqueue)]} {
		set dqueue $cgjob(dqueue)
	}
	if {$dqueue ne ""} {
		lappend options --partition=[join $dqueue ,]
	}
	if {$priority != 0} {
		lappend options --nice=$priority
	}
	if {[info exists cgjob(submitoptions)] && $cgjob(submitoptions) ne ""} {
		lappend options {*}$cgjob(submitoptions)
	}
	if {!$cgjob(nosubmit) && !$cgjob(dry)} {
		putslog "slurm_submit: [list sbatch --job-name=j$name --output=$job_out --error=$job_err --export=ALL {*}$options $runfile]"
		set jnum [exec sbatch --job-name=j$name --output=$job_out --error=$job_err --export=ALL {*}$options $runfile]
	} else {
		putslog "nosubmit run, would be slurm_submit: [list sbatch --job-name=j$name --output=$job_out --error=$job_err --export=ALL {*}$options $runfile]"
		set jnum [incr cgjob(nosubmit)]
	}
	regexp {[0-9]+} $jnum jobnum
	lappend cgjob(endjobids) $jobnum
	return $jobnum
}

proc slurm_safename {name {prefix {}}} {
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
	proc job_process_submit_slurm {job runfile args} {
		global curjobnum
		return [incr curjobnum]
	}
	interp alias {} job_process_submit_par {} job_process_submit_slurm
}

proc job_logfile_slurm_close {} {
	global cgjob
	set logfile $cgjob(logfile).running
	set runfile $logfile.update
	set outfile $logfile.update.out
	set errfile $logfile.update.err
#puts *****close*****
#putsvars cgjob(logfile) cgjob(cleanup) cgjob(cleanupfiles) cgjob(cleanupifemptyfiles)
	job_update $logfile $cgjob(cleanup) 1
#puts *****close_afterupdate*****
#putsvars cgjob(logfile) cgjob(cleanupfiles)
	set statusok [file exists $cgjob(logfile).finished]
#putsvars statusok
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

proc job_wait_slurm {} {
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
	append cmd "exec [job_timecmd $name] $cgjob(runcmd) \"\$0\" \"\$@\"\n\n"
	append cmd "job_init -d slurm\n"
	append cmd [list array set cgjob [array get cgjob]]\n
	# append cmd [list putsvars wait cgjob(cleanupfiles) cgjob(logfile)]\n
	append cmd job_logfile_slurm_close\n
	file_write $runfile $cmd
	file attributes $runfile -permissions u+x
	set options {}
	if {[llength $cgjob(endjobids)]} {
		lappend options --dependency=afterany:[join $cgjob(endjobids) :]
	}
	if {[info exists cgjob(dqueue)]} {
		lappend options --partition=[join $cgjob(dqueue) ,]
	}
	if {$priority != 0} {
		lappend options --nice=$priority
	}
	if {[info exists cgjob(submitoptions)] && $cgjob(submitoptions) ne ""} {
		lappend options {*}$cgjob(submitoptions)
	}
	if {!$cgjob(nosubmit) && !$cgjob(dry)} {
		putslog "slurm_submit: [list sbatch --job-name=j$name -o $outfile -e $errfile --export=ALL {*}$options $runfile]"
		exec sbatch --job-name=j$name -o $outfile -e $errfile --export=ALL {*}$options $runfile
	} else {
		putslog "nosubmit run, would be slurm_submit: [list sbatch --job-name=j$name -o $outfile -e $errfile --export=ALL {*}$options $runfile]"
		set jnum [incr cgjob(nosubmit)]
	}
}

proc grid_wait_slurm {} {
	while 1 {
		after 500
		puts -nonewline .
		flush stdout
		if {[llength [split [string trim [exec squeue]] \n]] == 1} break
	}
	puts ""
}

# https://slurm.schedmd.com/sbatch.html
# https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/antwerp/SLURM_UAntwerp.html
# https://curc.readthedocs.io/en/latest/running-jobs/slurm-commands.html