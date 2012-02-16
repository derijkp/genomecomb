#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

set basedir [file normalize [pwd]]
#puts $basedir

if {![llength $argv]} {
	puts stderr "use\nsubmit.tcl ?-deps jobids? ?-host hostname? ?--? command ..."
	exit 1
}

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
		-io {
			lappend options -l io=$value
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
if {[string length $name] > 200} {
	set name [string range $name 0 100]....[string range $name end-100 end]
}
set tasknum {}
if {[info exists ra([list $name $tasknum])]} {
	puts "Job $name.$tasknum is running, skipping"
} else {
	file mkdir osge
	file mkdir esge
	set jnum [eval {exec qsub -N j$name -q all.q -o osge -e esge} $options [file normalize ~/bin/repeater.sh] [file normalize [pwd]] $argv]
	puts "$jnum $name"
}
