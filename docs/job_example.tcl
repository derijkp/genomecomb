#!/bin/sh
# the next line restarts using wish \
exec cg source "$0" ${1+"$@"}

proc main args {
	set num 2
	set max 10
	set pos 0
	foreach {key value} $args {
		switch -- $key {
			-num {
				set num $value
			}
			-max {
				set max $value
			}
			default break
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	set len [llength $args]
	if {$len != 1} {
		error "format is: job_example.tcl ?options? dir"
	}
	foreach {dir} $args break
	file mkdir $dir
	cd $dir
	#
	job_logdir log_jobs
	for {set i 0} {$i < $num} {incr i} {
		job initfile-$i -targets numbers-$i.txt -vars max -skip {sumtotal.txt message.txt} -code {
			# write to temp, so that when the code fails, we wont be stuck with an incomplete target file
			set f [open $target.temp w]
			for {set j 0} {$j < 10} {incr j} {
				puts $f [expr {rand()*$max}]
			}
			close $f
			# rename tempfile when finished (atomic operation)
			file rename -force $target.temp $target
		}
	}
	# separate jobs for each file following a pattern
	# these files are not necesarily created yet, they are targets from previous job
	job sum -foreach {^numbers-(.*).txt$} -targets {sum-\1.txt} -code {
		set list [file_read $dep]
		file_write $target.temp [lmath_sum $list]\n
		file rename -force $target.temp $target
	}
	job error_nonexistingdep.txt -deps {nonexistingdep.txt} -targets {error_nonexistingdep.txt} -code {
		# this depends on nonexistingdep.txt, since it does not exist, this will not be run
	}
	# multiple dependencies, multiple targets
	job sumtotal -deps {sum-*.txt} -targets {sumtotal.txt message.txt} -vars num -code {
		set total 0
		foreach dep $deps {
			set total [expr {$total + [file_read $dep]}]
		}
		file_write $target.temp $total\n
		file rename -force $target.temp $target
		file_write $target2.temp "job using $num files finished\n"
		file rename -force $target2.temp $target2
	}
}

set argv [job_init {*}$argv]
main {*}$argv
job_wait
