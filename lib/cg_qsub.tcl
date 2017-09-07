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
	
	set name $command.[join $args .]
	regsub -all / $name __ name
	if {[string length $name] > 200} {
		set name [string range $name 0 100]....[string range $name end-100 end]
	}
	set tasknum {}
	if {[info exists ra([list $name $tasknum])]} {
		puts "Job $name.$tasknum is running, skipping"
	} else {
		if {![info exists outputfile]} {set outputfile job_$name.out}
		if {![info exists errorfile]} {set errorfile job_$name.err}
		set jnum [exec qsub -N j$name -q $dqueue -o $outputfile -e $errorfile {*}$options $::appdir/lib/repeater.sh $basedir $command {*}$args]
		puts "$jnum $name"
	}
}
