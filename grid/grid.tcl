#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

set basedir [file normalize [pwd]]
puts $basedir

if {[lindex $argv 0] eq "-deps"} {
	set dep [lindex $argv 1]
	set argv [lrange $argv 2 end]
} else {
	set dep {}
}
set command [lindex $argv 0]
set files [lrange $argv 1 end]

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

foreach file $files {
	regsub / $file __ nfile
	set name [file tail $command].[file tail $nfile]
	set tasknum {}
	if {[info exists ra([list $name $tasknum])]} {
		puts "Job $name.$tasknum is running, skipping"
		continue
	}
	puts "Submitting $name"
	file mkdir sge
	if {$dep eq ""} {
		eval {exec qsub -N j$name -q all.q -o sge/ -e sge/} $command {$file}
	} else {
		eval {exec qsub -hold_jid $dep -N j$name -q all.q -o sge/ -e sge/} $command {$file}
	}
}

