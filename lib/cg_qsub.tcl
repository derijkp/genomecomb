proc cg_qsub {args} {
	set basedir [file_absolute [pwd]]
	set options {}
	set dqueue all.q
	cg_options qsub args {
		-deps {
			lappend options -hold_jid [join $value ,]
		}
		-host {
			lappend options -l hostname=$value
		}
		-io {
			lappend options -l io=$value
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
	} command	
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
	set jobname [string_change $jobname [list { } _ / __ \\ __ : _ > _ * _ ' _ ( _ ) _ , _ / __ \n __ * _ ? _]]
	if {[string length $jobname] > 200} {
		set jobname [string range $jobname 0 100]....[string range $jobname end-100 end]
	}
	set tasknum {}
	if {[info exists ra([list $jobname $tasknum])]} {
		puts "Job $jobname.$tasknum is running, skipping"
	} else {
		job_init
		global cgjob
		set cmd {#!/bin/sh}
		append cmd \n
		append cmd {#$ -S /bin/bash} \n
		append cmd {#$ -V} \n
		append cmd {#$ -cwd} \n
		append cmd "\n\# the next line restarts using runcmd (specialised tclsh) \\\n"
		append cmd "exec $cgjob(runcmd) \"\$0\" \"\$@\"\n"
		append cmd [list exec $command {*}$args]\n
		set runfile job_$jobname.run
		file_write $runfile $cmd
		file attributes $runfile -permissions u+x

		if {![info exists outputfile]} {set outputfile job_$jobname.out}
		if {![info exists errorfile]} {set errorfile job_$jobname.err}
		set jnum [exec qsub -N j$jobname -q $dqueue -o $outputfile -e $errorfile {*}$options $runfile]
		puts "$jnum $jobname"
	}
}
