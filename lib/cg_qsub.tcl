proc cg_qsub {args} {
	global cgjob_distr
	set basedir [file_absolute [pwd]]
	set options {}
	set dqueue all.q
	set soft {}
	set hard {}
	set priority 0
	set cores 1
	set mem {}
	set memlimit {}
	set time {}
	set lang {}
	set submitoptions {}
	cg_options qsub args {
		-deps {
			lappend options -hold_jid [join $value ,]
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
		}
		-hard {
			lappend hard -l $value
		}
		-soft {
			lappend soft -l $value
		}
		-host {
			lappend options -l hostname=$value
		}
		-io {
			lappend options -l io=$value
		}
		-mem {
			set mem $value
		}
		-memlimit {
			set memlimit $value
		}
		-time {
			set time $value
		}
		-o - -outputfile {
			set outputfile $value
		}
		-e - -errorfile {
			set errorfile $value
		}
		-queue - -dqueue {
			set dqueue $value
		}
		-name {
			set jobname $value
		}
		-run {
			set run $value
		}
		-lang {
			set lang $value
		}
		-submitoptions {
			lappend options {*}$value
		}
	} command	
	if {$mem ne ""} {
		# mem_free: only start job if the given amount of memory is free on the node
		#  -> often not sufficient if running jobs increase memory use during run
		# virtual_free: reserves this much memory, can be defined as a consumable resource
		#	-> if all memory is reserved for running kjobs, no new jobs are started on node
		# not using h_vmem, because that would kill any job going (even a bit) above reserved memory
		set temp [job_mempercore $mem $cores]
		lappend hard -l mem_free=$temp,virtual_free=$temp
	}
	if {$memlimit ne ""} {
		set temp [job_mempercore $memlimit $cores]
		lappend hard -l h_vmem=$temp
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
	# check running jobs
	catch {exec qstat -xml} jobxml
	set jobs [regexp -all -inline {<job_list.+?</job_list>} $jobxml]
	unset -nocomplain ra
	foreach job $jobs {
		set task {} ; set name {}
		regexp {<tasks>(.*?)</tasks>} $job temp task
		regexp {<JB_name>(.*?)</JB_name>} $job temp name
		set name [string range $name 1 end]
		if {[string is int $task]} {
			set ra([list $name $task]) 1
		} else {
			set tasks {}
			regexp {([0-9]+)-([0-9]+)} $task temp start end
			for {} {$start <= $end} {incr start} {
				set ra([list $name $start]) 1
			}
		}
	}
	if {![info exists jobname]} {
		set jobname $command.[join $args .]
	}
	if {![info exists run]} {
		set run $command.[lindex $args 0].[string_change [lindex [split [timestamp] :] 0] {" " _ : - . -}]
	}
	set jobname [sge_safename $jobname $run]
	set tasknum {}
	if {[info exists ra([list $jobname $tasknum])]} {
		puts "Job $jobname.$tasknum is running, skipping"
	} else {
		job_init
		global cgjob
		set runfile job_$jobname.run
		set cmd {#!/bin/sh}
		append cmd \n
		append cmd {#$ -S /bin/bash} \n
		append cmd {#$ -V} \n
		append cmd "\n\# the next line restarts using runcmd (specialised tclsh) \\\n"
		append cmd "exec $cgjob(runcmd) \"\$0\" \"\$@\"\n"
		append cmd "cd [pwd]\n"
		if {$lang eq ""} {
			append cmd "[list exec $command {*}$args] >@ stdout 2>@ stderr\n"
		} else {
			if {$lang eq "R"} {
				set lang "R --vanilla --slave --no-restore <"
			}
			append cmd {set tempfile [tempfile]} \n
			append cmd "file_write \$tempfile [list [deindent $command]]\n"
			append cmd "exec $lang \$tempfile >@ stdout 2>@ stderr\n"
		}
		file_write $runfile $cmd
		file attributes $runfile -permissions u+x

		if {![info exists outputfile]} {set outputfile job_$jobname.out}
		if {![info exists errorfile]} {set errorfile job_$jobname.err}
		puts "run file: $runfile"
		puts "output file: $outputfile"
		puts "error file: $errorfile"
		set jnum [exec qsub -N j$jobname -cwd -q $dqueue -o $outputfile -e $errorfile -p $priority {*}$options $runfile]
		puts "$jnum $jobname"
	}
}

