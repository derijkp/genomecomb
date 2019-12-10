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
	# set the job log file
	# first argument is the base name of the logfile
	# second is the base directory
	# next is the commandline (to record in the logfile)
	set cmdline [list job_example.tcl -num $num -max $max $dir]
	job_logfile $dir/example_[file tail $dir] $dir $cmdline
	#
	# set the job log directory
	# This is a directory where output, job scripts etc. will be stored
	# The combination of log dir and jobname (see further) must be unique
	# You can us different log dirs in one run (change the log dir during the program)
	job_logdir log_jobs
	for {set i 0} {$i < $num} {incr i} {
		# run a job
		job initfile-$i -targets {numbers-$i.txt} -vars max -skip {sumtotal.txt message.txt} -code {
			# write to temp, so that when the code fails, we wont be stuck with an incomplete target file
			set f [open $target.temp w]
			for {set j 0} {$j < 10} {incr j} {
				puts $f [expr {rand()*$max}]
			}
			close $f
			# rename tempfile when finished (atomic operation)
			file rename -force -- $target.temp $target
		}
	}
	# separate jobs for each file following a pattern
	# these files are not necesarily created yet, they are targets from previous job
	job sum -foreach {^numbers-(.*).txt$} -targets {sum-\1.txt} -code {
		set list [file_read $dep]
		file_write $target.temp [lmath_sum $list]\n
		file rename -force -- $target.temp $target
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
		file rename -force -- $target.temp $target
		file_write $target2.temp "job using $num files finished\n"
		file rename -force -- $target2.temp $target2
	}
}

# job_init initializes the job system.
# it accept the job options (-d, ...) and returns the rest of the arguments
set argv [job_init {*}$argv]

# run program that will submit/run jobs
main {*}$argv

# job_wait will wait (if needed) for all jobs to finish
# this is currently only useful for local distributed processing using workers (-d number)
# direct runs jobs directly, and sge submits jobs directly (with the needed dependencies)
job_wait
