#!/bin/sh
# the next line restarts using wish \
exec cg source "$0" ${1+"$@"}

# simple example of the job system in genomecomb
# You can get more info on the jobsystem with:
# cg help jobsystem
# You can get info on the extra options available to control the job system with:
# cg help joboptions

# you can run this script using:
# job_example.tcl resultdir
# or run the same program using some options:
# - using multiple (maximum 4) cores on the local machine (-d 4)
# - verbose: giving information about jobs being started (-v 1)
# - program parameter num set to 10 (-num 10): make more files than the default 5
# job_example.tcl -d 4 -v 2 -num 10 resultdir

# main body of code, called at the end of the file (after initialising the job system)
proc main args {
	set num 5
	set max 10
	set pos 0
	cg_options job_example args {
		-num {
			set num $value
		}
		-max {
			set max $value
		}
	} dir 1 1
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
	# This set_job_logdir is not really necessary, as job_logfile already sets the job logdir to $dir/log_jobs
	set_job_logdir log_jobs
	#
	# run the given number ($num) of jobs, creating files numbers-1.txt, numbers-2.txt, ...
	# if they do not exist yet
	# The code in the -code option/block may be run in parallel
	set cleanup {}
	for {set i 0} {$i < $num} {incr i} {
		# run a job
		job initfile-$i -targets {numbers-$i.txt} -vars max -code {
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
	# next runs separate jobs for each result of the previous jobs.
	# These files are not necesarily created yet (when running distributed), but are targets from previous job
	# The -code blocks will only be run when the jobs that create the dependencies (numbers-1.txt, ...) are finished
	for {set i 0} {$i < $num} {incr i} {
		job sum-$i -deps {numbers-$i.txt} -targets {sum-$i.txt} -code {
			set list [file_read $dep]
			file_write $target.temp [lmath_sum $list]\n
			file rename -force -- $target.temp $target
		}
	}
# next (commented out) would give an error, as it is a required job where a dependency does not exist, 
# nor is being made by one of the previous jobs
#	job error_nonexistingdep.txt -deps {nonexistingdep.txt} -targets {error_nonexistingdep.txt} -code {
#		# this depends on nonexistingdep.txt, since it does not exist, this will not be run
#	}
	# The -optional 1 indicates that if dependencies are not met, the job can be skipped
	job optional_nonexistingdep.txt -optional 1 -deps {nonexistingdep.txt} -targets {error_nonexistingdep.txt} -code {
		# this depends on nonexistingdep.txt, since it does not exist, this will not be run
	}
	# multiple dependencies, multiple targets
	job sumtotal -deps {sum-*.txt} -targets {sumtotal.txt message.txt} -vars {num cleanup} -code {
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

# job_wait will wait (if needed) for all jobs to finish, cleanup temporary files and update the log file
job_wait
