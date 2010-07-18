#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

set basedir [file normalize [pwd]]
puts $basedir

set options {}
set pos 0
foreach {opt value} $argv {
	switch -- $opt {
		-deps {
			lappend options -hold_jid $value
			incr pos 2
		}
		-host {
			lappend options -l hostname=$value
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
if {$pos} {
	set argv [lrange $argv $pos end]
}

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

set name [join $argv .]
regsub -all / $name __ name
set tasknum {}
if {[info exists ra([list $name $tasknum])]} {
	puts "Job $name.$tasknum is running, skipping"
	continue
}
puts "Submitting $name"
file mkdir osge
file mkdir esge
eval {exec qsub -N j$name -q all.q -o osge/ -e esge/} $options [file normalize ~/bin/repeater.sh] $argv

