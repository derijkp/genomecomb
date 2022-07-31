#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

if {![info exists argv]} {set argv {}}
set pos [lsearch $argv -d]
if {$pos != -1} {
	set distribute [lindex $argv [incr pos]]
	if {$distribute eq "direct"} {
		puts "testing direct only"
		set tests {
			"direct" {uplevel job_init -skipjoberrors 1 {*}\$args}
		}
	} elseif {[string is int $distribute]} {
		puts "testing -d $distribute only"
		set tests [subst {
			"-d $distribute" {uplevel job_init -d $distribute {*}\$args}
		}]
	} elseif {$distribute eq "sge"} {
		puts "testing -d sge only"
		set testsge 1
		set tests [subst {
			"-d sge" {uplevel job_init -d sge {*}\$args}
		}]
	} elseif {$distribute eq "slurm"} {
		puts "testing -d slurm only"
		set testslurm 1
		set tests [subst {
			"-d slurm" {uplevel job_init -d slurm {*}\$args}
		}]
	} else {
		error "unknown option -d $distribute"
	}
} else {
	set tests {
		"direct" {uplevel job_init -skipjoberrors 1 -d 0 {*}$args}
		"-d 2" {uplevel job_init -d 2 {*}$args}
		"-d 4" {uplevel job_init -d 4 {*}$args}
		"-d 30" {uplevel job_init -d 30 {*}$args}
		"-d sge" {uplevel job_init -d sge {*}$args}
		"-d slurm" {uplevel job_init -d slurm {*}$args}
	}
}

source tools.tcl
set keepdir [pwd]

# use these for trying out individual tests
set testname "-d direct"
proc test_job_init {args} {
	uplevel job_init -skipjoberrors 1 -d 0 {*}$args
	uplevel job_logfile_set $::testdir/tmp/log $::testdir/tmp
}
proc gridwait {} {}

if 0 {
	set testname "-d direct"
	proc test_job_init {args} {
		uplevel job_init -skipjoberrors 1 {*}$args
		uplevel job_logfile_set $::testdir/tmp/log $::testdir/tmp
	}
	interp alias {} job_wait {} job_wait_direct
	proc gridwait {} {}

	set testname "-d 2"
	proc test_job_init {args} {
		uplevel job_init -d 2 {*}$args
		uplevel job_logfile_set $::testdir/tmp/log $::testdir/tmp
	}
	interp alias {} job_wait {} job_wait_distr
	proc gridwait {} {}

	set testname "-d 4"
	proc test_job_init {args} {
		uplevel job_init -d 4 {*}$args
		uplevel job_logfile_set $::testdir/tmp/log $::testdir/tmp
	}
	interp alias {} job_wait {} job_wait_distr
	proc gridwait {} {}

	set testname "-d 30"
	proc test_job_init {args} {
		uplevel job_init -d 30 {*}$args
		uplevel job_logfile_set $::testdir/tmp/log $::testdir/tmp
	}
	interp alias {} job_wait {} job_wait_distr
	proc gridwait {} {}

	set testname "-d sge"
	proc test_job_init {args} {
		uplevel job_init -d sge {*}$args
		uplevel job_logfile_set $::testdir/tmp/log $::testdir/tmp
	}
	interp alias {} job_wait {} job_wait_sge
	proc gridwait {} {
		while 1 {
			after 500
			puts -nonewline .
			flush stdout
			if {[exec qstat] eq ""} break
		}
	}

	set testname "-d slurm"
	proc test_job_init {args} {
		uplevel job_init -d slurm {*}$args
		uplevel job_logfile_set $::testdir/tmp/log $::testdir/tmp
	}
	interp alias {} job_wait {} job_wait_slurm
	proc gridwait {} {
		while 1 {
			after 500
			puts -nonewline .
			flush stdout
			if {[llength [split [string trim [exec squeue]] \n]] == 1} break
		}
	}
}

proc writetestfiles {args} {
	test_cleantmp
	foreach file $args {
		after 1000
		file_write $file $file
	}
}

proc jobtest {args} {
	set args [job_args $args]
	foreach {srcdir destdir header unknowntargets} $args break
	if {$unknowntargets eq "nounknowntargets"} {set unknowntargets 0} else {set unknowntargets 1}
	set srcdir [file_absolute $srcdir]
	set destdir [file_absolute $destdir]
	set job_logdir $destdir/log_jobs
	global test_names
	file mkdir $destdir
	set file [lindex [glob $srcdir/cgdat*.tsv] 0]
	set f [open $file]
	set names {}
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		lappend names [lindex $line 0]
	}
	close $f
	job sumvar -deps {
		$srcdir/cgdat*.tsv
	} -targets {
		$destdir/sum-$_names.txt
	} -vars {destdir} -code {
		for {set i 1} {$i < 5} {incr i} {
			# puts "progress $target $i"
			after 250
		}
		set targets {}
		set f [open [glob $dep]]
		while {![eof $f]} {
			set line [split [gets $f] \t]
			if {$line eq ""} continue
			# puts line=$line
			set name [lindex $line 0]
			set target $destdir/sum-$name.txt
			set o [open $target.temp w]
			set calc [join [lrange $line 1 end] +]
			puts $o $calc=[expr $calc]
			close $o
			file rename -force -- $target.temp $target
			lappend targets $target
		}
		close $f
	}
	if {!$unknowntargets} {
		set targets {$destdir/sumpattern.log $destdir/sumpattern-test1.txt $destdir/sumpattern-test2.txt $destdir/sumpattern-test3.txt}
	} else {
		set targets {$destdir/sumpattern.log}
	}
	job sumpattern -deps {$srcdir/cgdat*.tsv} \
	-targets $targets \
	-vars destdir -code {
		# var target contains target
		# var target1 contains first braced part of target (= destdir)
		# var target2 contains second braced part of target (= name)
		# var deps contains a list of all dependencies
		# var dep contains the first element of the first dependency (so you do not have to do [lindex $deps 0] to get it)
		for {set i 1} {$i < 5} {incr i} {
			# puts "progress $target $i"
			after 250
		}
		set targets {}
		set f [open [glob $dep]]
		while {![eof $f]} {
			set line [split [gets $f] \t]
			if {$line eq ""} continue
			# puts line=$line
			set name [lindex $line 0]
			set target $destdir/sumpattern-$name.txt
			set o [open $target.temp w]
			set calc [join [lrange $line 1 end] +]
			puts $o "p $calc=[expr $calc]"
			close $o
			file rename -force -- $target.temp $target
			lappend targets $target
		}
		close $f
		file_write $destdir/sumpattern.log [join $targets \n]
	}
	if {$unknowntargets} {
		job_runall
	}
	set files [jobglob $destdir/sumpattern-*.txt]
	foreach file $files {
		regexp {sumpattern-(.*).txt$} [file tail $file] temp base
		set target $destdir/sumpattern2-$base.txt
		job sumpattern2-$base -cores 2 -skip {
			$destdir/allp2.txt
		} -deps {
			$file
		} -targets {
			$target
		} -code {
			for {set i 1} {$i < 5} {incr i} {
				# puts "progress $i"
				after 250
			}
			file copy $dep $target.temp
			exec echo 2 >> $target.temp
			file rename -force -- $target.temp $target
		}
	}
	job allp2.txt -vars header -deps {$destdir/sumpattern2-*.txt} -targets {$destdir/allp2.txt} -code {
		after 500
		file_write $target.temp $header\n
		exec cat {*}$deps >> $target.temp
		file rename -force -- $target.temp $target
	}
	job allp.txt -vars header -deps {^$destdir/sumpattern-(.*)\.txt$} -targets {$destdir/allp.txt} -code {
		file_write $target.temp $header\n
		exec cat {*}$deps >> $target.temp
		file rename -force -- $target.temp $target
	}
	job test -deps {$srcdir/cgdat*.tsv} -targets {$destdir/test.txt} -code {
		exec wc $dep > $target
	}
	foreach file [jobglob $destdir/sum-*.txt] {
		regexp {sum-(.*).txt$} [file tail $file] temp base
		job sum2-$base -skip [list $destdir/all2.txt] -deps {
			$file
		} -targets {
			$destdir/sum2-$base.txt
		} -code {
			for {set i 1} {$i < 5} {incr i} {
				# puts stderr "progress $i"
				after 250
			}
			file copy $dep $target.temp
			exec echo 2 >> $target.temp
			file rename -force -- $target.temp $target
		}
	}
	job error_all.txt -optional 1 -deps {$srcdir/notpresent.txt} -targets {$destdir/all.txt} -code {
		error "This should not be executed, as the dependencies are not fullfilled, the other target is used"
	}
	job joberror -deps {^$srcdir/cgdata\.tsv$} -targets {$destdir/joberror.txt} -code {
		# this (intentionally) causes an error, and should be identified by an error in the logs
		error "Intentional job error"
	}
	job error_notdone.txt -deps {^$srcdir/cgdata\.tsv$} -targets {$destdir/error_notdone.txt} -code {
		# this (intentionally) does not produce its target, and should by identified by an error in the logs
	}
	job error_notdonedep.txt -deps {$destdir/error_notdone.txt} -targets {$destdir/error_notdonedep.txt} -code {
		# this depends on error_notdone.txt, which will fail to produce its target
		# backends should cope with this gracefully (e.g. not hanging)
	}
	job all.txt -vars header -deps {^$destdir/sum-(.*)\.txt$} -targets {$destdir/all.txt} -code {
		file_write $target.temp $header\n
		exec cat {*}[bsort $deps] >> $target.temp
		file rename -force -- $target.temp $target
	}
	set keepdir [pwd]
	cd $destdir
	job all2.txt -vars header -deps {^sum2-(.*)\.txt$} -targets {all2.txt} -code {
		file_write $target.temp $header\n
		exec cat {*}[bsort $deps] >> $target.temp
		file rename -force -- $target.temp $target
	}
	cd $keepdir
}

proc jobtestnojobs {args} {
	set args [job_args $args]
	foreach destdir $args break
	set_job_logdir $destdir/log_jobs
	job nodeps1 -optional 1 -deps {abcd} -targets {abcde} -code {
	}
	job nodeps2 -optional 1 -deps {abcde} -targets {abcdef} -code {
	}
}

proc jobtestlong {args} {
	set args [job_args $args]
	foreach destdir $args break
	file mkdir $destdir
	cd $destdir
	set_job_logdir $destdir/log_jobs
	job long -deps {} -targets {long.txt} -code {
		for {set i 0} {$i < 4} {incr i} {
			after 1000
			# puts $i
		}
		file_write long.txt done
	}
}

test job {job_expandvars} {
	unset -nocomplain ::a; unset -nocomplain ::b; unset -nocomplain ::c; unset -nocomplain ::d;
	set string {$a($b($c))}
	set a(B(C)) A(B(C))
	set b(C) B(C)
	set c C
	set d D
	job_expandvars $string
} {A(B(C))}

test job {job_expandvars} {
	unset -nocomplain ::a; unset -nocomplain ::b; unset -nocomplain ::c; unset -nocomplain ::d;
	set string {$a($b($_c))}
	set a(B(C1)) A(B(C1))
	set a(B(C2)) A(B(C2))
	set b(C1) B(C1)
	set b(C2) B(C2)
	set c {C1 C2}
	set d D
	job_expandvars $string
} {A(B(C1)) A(B(C2))}

test job {job_expandvars} {
	unset -nocomplain ::a; unset -nocomplain ::b; unset -nocomplain ::c; unset -nocomplain ::d;
	set string {$a($_b($_c))}
	set a(B1(C1)) A(B1(C1))
	set a(B2(C1)) A(B2(C1))
	set a(B(C2)) A(B(C2))
	set b(C1) {B1(C1) B2(C1)}
	set b(C2) B(C2)
	set c {C1 C2}
	set d D
	job_expandvars $string
} {A(B1(C1)) A(B2(C1)) A(B(C2))}

test job {job_expandvarslist} {
	set a {test it}
	set b now
	job_expandvarslist {$a$b $a$a}
} {{test itnow} {test ittest it}}

# ------------------------------------
# test in different "processing modes"
# ------------------------------------
foreach {testname initcode} $tests {
# start of block

if {$testname eq "-d sge"} {
	if {![get testsge 0]} {
		puts "skipping sge tests because testsge is not set (default), set testsge 1 or give testsge as a parameter to all.tcl to run"
		continue
	}
	if {[catch {exec qstat}]} {
		puts "Cannot test sge option (missing qstat; grid engine not installed?)"
		continue
	}
	interp alias {} job_wait {} job_wait_sge
	proc gridwait {} {
		while 1 {
			after 500
			puts -nonewline .
			flush stdout
			if {[exec qstat] eq ""} break
		}
	}
} elseif {$testname eq "-d slurm"} {
	if {![get testslurm 0]} {
		puts "skipping slurm tests because testslurm is not set (default), set testslurm 1 or give testslurm as a parameter to all.tcl to run"
		continue
	}
	if {[catch {exec squeue}]} {
		puts "Cannot test slurm option (missing squeue; slurm not installed?)"
		continue
	}
	interp alias {} job_wait {} job_wait_slurm
	proc gridwait {} {
		while 1 {
			after 500
			puts -nonewline .
			flush stdout
			if {[llength [split [string trim [exec squeue]] \n]] == 1} break
		}
	}
} elseif {[regexp -- {-d [0-9]} $testname]} {
	interp alias {} job_wait {} job_wait_distr
} else {
	interp alias {} job_wait {} job_wait_direct
	proc gridwait {} {}
}

# putsvars testname initcode

proc test_job_init {args} "$initcode\nuplevel job_logfile_set \$::testdir/tmp/log \$::testdir/tmp"

test job "basic chain $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	file_write test1.txt test1\n
	job job1 -deps {test1.txt} -targets {test2.txt} -code {
		set c [file_read $dep]
		file_write $target ${c}test2\n
	}
	job job2 -deps {test2.txt} -targets {test3.txt} -code {
		set c [file_read $dep]
		file_write $target ${c}test3\n
	}
	job_wait
	gridwait
	set result [list [bsort [glob *]] [file_read test3.txt]]
	cd $::testdir
	set result
} {{log.*.finished test1.txt test2.txt test3.txt} {test1
test2
test3
}} match

test job "foreach $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	file_write test1.txt test1\n
	file_write test2.txt test2\n
	foreach file [jobglob test*.txt] {
		regexp {test(.*).txt$} $file temp base
		set target rtest$base.txt
		job job1-$base -deps {$file} -targets {$target} -code {
			set c [file_read $dep]
			file_write $target r${c}
		}
	}
	job_wait
	gridwait
	set result [list [bsort [glob *]] [file_read rtest1.txt]]
	cd $::testdir
	set result
} {{log.*.finished rtest1.txt rtest2.txt test1.txt test2.txt} {rtest1
}} match

test job "chained foreach $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	job init -deps {} -targets {test1.txt test2.txt} -code {
		file_write test1.txt test1\n
		file_write test2.txt test2\n
	}
	foreach file [jobglob test*.txt] {
		job job-$file -deps {$file} -targets {r$file} -code {
			set c [file_read $dep]
			file_write $target r${c}
		}
	}
	job job2 -deps {rtest1.txt} -targets {final1.txt} -code {
		set c [file_read $dep]
		file_write $target f${c}
	}
	job_wait
	gridwait
	set result [list [bsort [glob *]] [file_read final1.txt]]
	cd $::testdir
	set result
} {{final1.txt log.*.finished rtest1.txt rtest2.txt test1.txt test2.txt} {frtest1
}} match

test job "chained foreach with glob match $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	job init -deps {} -targets {test1.txt test2.txt} -code {
		file_write test1.txt test1\n
		file_write test2.txt test2\n
	}
	foreach file [jobglob test*.txt] {
		regsub {test(.*).txt} $file {rtest\1.txt} target
		job job1-$file -deps {$file} -targets {$target} -code {
			set c [file_read $dep]
			file_write $target r${c}
		}
	}
	job job2 -deps {rtest1.txt} -targets {final1.txt} -code {
		set c [file_read $dep]
		file_write $target f${c}
	}
	job_wait
	gridwait
	set result [list [bsort [glob *]] [file_read final1.txt]]
	cd $::testdir
	set result
} {{final1.txt log.*.finished rtest1.txt rtest2.txt test1.txt test2.txt} {frtest1
}} match

test job "chained jobglob $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	job init -deps {} -targets {test1.txt test2.txt} -code {
		after 500
		file_write test1.txt test1\n
		file_write test2.txt test2\n
	}
	foreach dep [jobglob test*.txt] {
		set target r$dep
		job job1-[file tail $dep] -deps {$dep} -targets {$target} -code {
			set c [file_read $dep]
			file_write $target r${c}
		}
	}
	job job2 -deps {rtest1.txt} -targets {final1.txt} -code {
		set c [file_read $dep]
		file_write $target f${c}
	}
	job_wait
	gridwait
	set result [list [bsort [glob *]] [file_read final1.txt]]
	cd $::testdir
	set result
} {{final1.txt log.*.finished rtest1.txt rtest2.txt test1.txt test2.txt} {frtest1
}} match

test job "chained jobglob spaces in name $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	set target1 {test 1.txt}
	job init -deps {} -targets {$target1 {test 2.txt}} -code {
		after 500
		file_write {test 1.txt} test1\n
		file_write {test 2.txt} test2\n
	}
	foreach dep [jobglob test*.txt] {
		set target r$dep
		job job1-[file tail $dep] -deps {$dep} -targets {$target} -code {
			set c [file_read $dep]
			file_write $target r${c}
		}
	}
	job job2 -deps {{rtest 1.txt}} -targets {{final 1.txt}} -code {
		set c [file_read $dep]
		file_write $target f${c}
	}
	job_wait
	gridwait
	set result [list [bsort [glob *]] [file_read {final 1.txt}]]
	cd $::testdir
	set result
} {{{final 1.txt} log.*.finished {rtest 1.txt} {rtest 2.txt} {test 1.txt} {test 2.txt}} {frtest1
}} match

test job "basic chain --force 0 $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	job_args {--force 0}
	file_write test1.txt test1\n
	file_write test2.txt error2\n
	file_write test3.txt error3\n
	job job1 -deps {test1.txt} -targets {test2.txt} -code {
		set c [file_read $dep]
		file_write $target ${c}test2\n
	}
	job job2 -deps {test2.txt} -targets {test3.txt} -code {
		set c [file_read $dep]
		file_write $target ${c}test3\n
	}
	job_wait
	gridwait
	set result [list [bsort [glob *]] [file_read test3.txt]]
	cd $::testdir
	set result
} {{log.*.finished test1.txt test2.txt test3.txt} {error3
}} match

test job "basic chain --force 1 $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	job_args {--force 1}
	file_write test1.txt test1\n
	file_write test2.txt error2\n
	file_write test3.txt error3\n
	job job1 -deps {test1.txt} -targets {test2.txt} -code {
		set c [file_read $dep]
		file_write $target ${c}test2\n
	}
	job job2 -deps {test2.txt} -targets {test3.txt} -code {
		set c [file_read $dep]
		file_write $target ${c}test3\n
	}
	job_wait
	gridwait
	set result [list [bsort [glob *]] [file_read test2.txt] [file_read test3.txt]]
	cd $::testdir
	set result
} {{log.*.finished test1.txt test2.txt test3.txt} {test1
test2
} {test1
test2
test3
}} match

test job "time chain $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	file_write test1.txt test1\n
	after 1000
	file_write test3.txt error3\n
	after 1000
	file_write test2.txt error2\n
	job job1 -deps {test1.txt} -targets {test2.txt} -code {
		set c [file_read $dep]
		file_write $target ${c}test2\n
	}
	job job2 -deps {test2.txt} -targets {test3.txt} -code {
		set c [file_read $dep]
		file_write $target ${c}test3\n
	}
	job_wait
	gridwait
	set result [list [bsort [glob *]] [file_read test2.txt] [file_read test3.txt] [file_read test3.txt.old]]
	cd $::testdir
	set result
} {{log.*.finished test1.txt test2.txt test3.txt test3.txt.old} {error2
} {error2
test3
} {error3
}} match

test job "missing dep -optional 1 $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init -skipjoberrors 0
	file_write test1.txt test1\n
	job jobmissing -optional 1 -deps {missing.txt} -targets {test3.txt} -code {
		file_write $target missing\n
	}
	job job1 -deps {test1.txt} -targets {test2.txt} -code {
		set c [file_read $dep]
		file_write $target ${c}test2\n
	}
	job_wait
	gridwait
	set result [list [bsort [glob *]] [file_read test2.txt]]
	cd $::testdir
	set result
} {{log.*.finished test1.txt test2.txt} {test1
test2
}} match

test job "missing dep -optional 0 $testname" {
	test_cleantmp
	cd $::testdir/tmp
	test_job_init -skipjoberrors 0
	file_write test1.txt test1\n
	job jobmissing -optional 0 -deps {missing.txt} -targets {test3.txt} -code {
		file_write $target missing\n
	}
	job job1 -deps {test1.txt} -targets {test2.txt} -code {
		set c [file_read $dep]
		file_write $target ${c}test2\n
	}
	job_wait
	gridwait
	set result [list [bsort [glob *]] [file_read test2.txt]]
	cd $::testdir
	set result
} {error trying to run job jobmissing:
missing dependency "missing.txt"} error

test job "missing dep -optional 0 -skipjoberrors 1 $testname" {
	test_cleantmp
	cd $::testdir/tmp
	test_job_init -skipjoberrors 1
	file_write test1.txt test1\n
	job jobmissing -optional 0 -deps {missing.txt} -targets {test3.txt} -code {
		file_write $target missing\n
	}
	job job1 -deps {test1.txt} -targets {test2.txt} -code {
		set c [file_read $dep]
		file_write $target ${c}test2\n
	}
	job_wait
	gridwait
	set result [list [bsort [glob *]] [file_read test2.txt]]
	cd $::testdir
	set result
} {{log.*.finished test1.txt test2.txt} {test1
test2
}} match

test job "basic jobtest $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init -skipjoberrors 1
	jobtest ../data test testh
	job_wait
	gridwait
	set result [list \
		[bsort [glob test/*]] \
		[glob test/log_jobs/all.txt.log] \
		[file_read test/all.txt] \
		[file_read test/sum2-test3.txt] \
		[string range [file_read test/log_jobs/joberror.err] 0 20] \
	]
	cd $::testdir
	set result
} {{test/all2.txt test/all.txt test/allp2.txt test/allp.txt test/log_jobs test/sum2-test1.txt test/sum2-test2.txt test/sum2-test3.txt test/sum-test1.txt test/sum-test2.txt test/sum-test3.txt test/sumpattern2-test1.txt test/sumpattern2-test2.txt test/sumpattern2-test3.txt test/sumpattern-test1.txt test/sumpattern-test2.txt test/sumpattern-test3.txt test/sumpattern.log test/test.txt} test/log_jobs/all.txt.log {testh
1+2=3
3+4+5=12
6+7+8=21
} {6+7+8=21
2
} {Intentional job error}}

test job "basic jobtest nounknowntargets $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init -skipjoberrors 1
	jobtest ../data test testh nounknowntargets
	job_wait
	gridwait
	set result [list \
		[bsort [glob test/*]] \
		[glob test/log_jobs/all.txt.log] \
		[file_read test/all.txt] \
		[file_read test/sum2-test3.txt] \
		[string range [file_read test/log_jobs/joberror.err] 0 20] \
	]
	cd $::testdir
	set result
} {{test/all2.txt test/all.txt test/allp2.txt test/allp.txt test/log_jobs test/sum2-test1.txt test/sum2-test2.txt test/sum2-test3.txt test/sum-test1.txt test/sum-test2.txt test/sum-test3.txt test/sumpattern2-test1.txt test/sumpattern2-test2.txt test/sumpattern2-test3.txt test/sumpattern-test1.txt test/sumpattern-test2.txt test/sumpattern-test3.txt test/sumpattern.log test/test.txt} test/log_jobs/all.txt.log {testh
1+2=3
3+4+5=12
6+7+8=21
} {6+7+8=21
2
} {Intentional job error}}

test job "basic status $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	jobtest ../data test testh nounknowntargets
	job_wait
	gridwait
	set result [list \
		[bsort [glob test/*]] \
		[glob test/log_jobs/all.txt.log] \
		[file_read test/all.txt] \
		[file_read test/sum2-test3.txt] \
		[string range [file_read test/log_jobs/joberror.err] 0 20] \
	]
	cd $::testdir
	set result
	# status
	cd $::testdir/tmp
	job_init -d status
	jobtest ../data test testh nounknowntargets
	job_wait
	glob jobdeps.dot jobdeps.ps
} {jobdeps.dot jobdeps.ps}

test job "--force 0 $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	file mkdir test
	job_init -skipjoberrors
	jobtest --force 0 ../data test testh
	job_wait
	gridwait
	after 1000
	file_write test/all.txt error
	job_init -skipjoberrors
	jobtest --force 0 ../data test testh
	job_wait
	gridwait
	set result [list \
		[bsort [glob test/*]] \
		[file_read test/all.txt] \
		[file_read test/sum2-test3.txt] \
	]
	cd $::testdir
	set result
} {{test/all2.txt test/all.txt test/allp2.txt test/allp.txt test/log_jobs test/sum2-test1.txt test/sum2-test2.txt test/sum2-test3.txt test/sum-test1.txt test/sum-test2.txt test/sum-test3.txt test/sumpattern2-test1.txt test/sumpattern2-test2.txt test/sumpattern2-test3.txt test/sumpattern-test1.txt test/sumpattern-test2.txt test/sumpattern-test3.txt test/sumpattern.log test/test.txt} error {6+7+8=21
2
}}

test job "--force 1 $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	file mkdir test
	job_init -skipjoberrors
	jobtest --force 0 ../data test testh
	job_wait
	gridwait
	after 1000
	file_write test/all.txt error
	test_job_init
	jobtest --force 1 ../data test testh
	job_wait
	gridwait
	set result [list \
		[bsort [glob test/*]] \
		[file_read test/all.txt] \
		[file_read test/sum2-test3.txt] \
	]
	cd $::testdir
	set result
} {{test/all2.txt test/all.txt test/allp2.txt test/allp.txt test/log_jobs test/sum2-test1.txt test/sum2-test2.txt test/sum2-test3.txt test/sum-test1.txt test/sum-test2.txt test/sum-test3.txt test/sumpattern2-test1.txt test/sumpattern2-test2.txt test/sumpattern2-test3.txt test/sumpattern-test1.txt test/sumpattern-test2.txt test/sumpattern-test3.txt test/sumpattern.log test/test.txt} {testh
1+2=3
3+4+5=12
6+7+8=21
} {6+7+8=21
2
}}

test job "time $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	jobtest ../data test testh
	job_wait
	gridwait
	after 1000
	file_write $::testdir/tmp/test/sum-test3.txt "replaced\n"
	test_job_init
	jobtest ../data test testh
	job_wait
	gridwait
	set result [list \
		[bsort [glob test/*]] \
		[file_read test/all.txt] \
		[file_read test/sum2-test3.txt] \
	]
	cd $::testdir
	set result
} {{test/all2.txt test/all2.txt.old test/all.txt test/all.txt.old test/allp2.txt test/allp.txt test/log_jobs test/sum2-test1.txt test/sum2-test2.txt test/sum2-test3.txt test/sum2-test3.txt.old test/sum-test1.txt test/sum-test2.txt test/sum-test3.txt test/sumpattern2-test1.txt test/sumpattern2-test2.txt test/sumpattern2-test3.txt test/sumpattern-test1.txt test/sumpattern-test2.txt test/sumpattern-test3.txt test/sumpattern.log test/test.txt} {testh
1+2=3
3+4+5=12
replaced
} {replaced
2
}}

test job "-skip: not present $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	job testskip -deps {} -targets result.txt -skip {skip1.txt skip2.txt} -code {
		file_write result.txt test
	}
	job_wait
	gridwait
	set result [bsort [glob *]]
	cd $::testdir
	set result
} {log.*.finished result.txt} match

test job "-skip: only one present $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	file_write skip1.txt test1
	job testskip -deps {} -targets result.txt -skip {skip1.txt skip2.txt} -code {
		file_write result.txt test
	}
	job_wait
	gridwait
	set result [bsort [glob *]]
	cd $::testdir
	set result
} {log.*.finished result.txt skip1.txt} match

test job "-skip: all present $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	file_write skip1.txt test1
	file_write skip2.txt test2
	job testskip -deps {} -targets result.txt -skip {skip1.txt skip2.txt} -code {
		file_write result.txt test
	}
	job_wait
	gridwait
	set result [bsort [glob *]]
	cd $::testdir
	set result
} {log.*.finished skip1.txt skip2.txt} match

test job "-skip -skip: none present $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	job testskip -deps {} -targets result.txt -skip skip1.txt -skip skip2.txt -code {
		file_write result.txt test
	}
	job_wait
	gridwait
	set result [bsort [glob *]]
	cd $::testdir
	set result
} {log.*.finished result.txt} match

test job "-skip -skip: one present $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	file_write skip2.txt test2
	job testskip -deps {} -targets result.txt -skip skip1.txt -skip skip2.txt -code {
		file_write result.txt test
	}
	job_wait
	gridwait
	set result [bsort [glob *]]
	cd $::testdir
	set result
} {log.*.finished skip2.txt} match

test job "-skip: chain $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init -skipjoberrors 0
	file_write result.txt result
	job step1 -deps {} -targets step1.txt -skip {result.txt} -code {
		file_write step1.txt test
	}
	job result -deps step1.txt -targets result.txt -code {
		file_write result.txt [file_read $dep]
	}
	job_wait
	gridwait
	set result [bsort [glob *]]
	cd $::testdir
	set result
} {log.*.finished result.txt} match

test job "-skip: sone deps do not exist" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init -skipjoberrors 0
	file_write result.txt result
	job step1 -deps {notpresent.txt} -targets step1.txt -skip {result.txt} -code {
		file_write step1.txt test
	}
	job_wait
	gridwait
	set result [bsort [glob *]]
	cd $::testdir
	set result
} {log.*.finished result.txt} match

test job "-skip: sone deps are newer -> do not move skiptarget to .old" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init -skipjoberrors 0
	file_write dep.txt dep
	file_write skipresult.txt result
	file mtime dep.txt [expr {[file mtime skipresult.txt]+1}]
	job step1 -deps {dep.txt} -targets result.txt -skip {skipresult.txt skipresult2.txt} -code {
		file_write result.txt test
	}
	job_wait
	gridwait
	set result [glob skipresult.txt*]
	cd $::testdir
	set result
} {skipresult.txt}

test job "jobtestnojobs $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	jobtestnojobs $::testdir/tmp
	job_wait
	gridwait
} {}

test job "jobtestlong $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	jobtestlong $::testdir/tmp
	job_wait
	gridwait
} {}

test job "no -targets $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	file_write dep.txt dep
	job testnotarget -deps {dep.txt} -code {
		file_write result.txt test
	}
	job_wait
	gridwait
	set result [bsort [glob *]]
	cd $::testdir
	set result
} {dep.txt log.*.finished result.txt} match

test job "no -targets dep not found $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	job testnotarget -optional 1 -deps {dep.txt} -code {
		file_write result.txt test
	}
	job_wait
	gridwait
	set result [bsort [glob *]]
	cd $::testdir
	set result
} {log.*.finished} match

test job "no -targets dep not found not optional $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init -skipjoberrors 0
	job testnotarget -deps {dep.txt} -code {
		file_write result.txt test
	}
	job_wait
	gridwait
	set result [bsort [glob *]]
	cd $::testdir
	set result
} {error trying to run job testnotarget:
missing dependency "dep.txt"} error

test job "no -checkcompressed 1 (default) dep $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	file_write dep.txt test
	cg razip dep.txt
	job testcheckcompressed -deps {dep.txt} -targets result.txt -code {
		file_write result.txt test
	}
	job_wait
	gridwait
	set result [bsort [glob *]]
	cd $::testdir
	set result
} {dep.txt.rz log.*.finished result.txt} match

test job "no -checkcompressed 0 dep $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init -skipjoberrors 0
	file_write dep.txt test
	cg razip dep.txt
	job testcheckcompressed -checkcompressed 0 -deps {dep.txt} -targets result.txt -code {
		file_write result.txt test
	}
	job_wait
	gridwait
	set result [bsort [glob *]]
	cd $::testdir
	set result
} {error trying to run job testcheckcompressed:
missing dependency "dep.txt"} error

test job "no -checkcompressed 1 (default) targets $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	file_write dep.txt test
	file_write target.txt test
	cg razip target.txt
	job testcheckcompressed -deps {dep.txt} -targets target.txt -code {
		file_write result.txt test
		file_write target.txt test
	}
	job_wait
	gridwait
	set result [bsort [glob *]]
	cd $::testdir
	set result
} {dep.txt log.*.finished target.txt.rz} match

test job "no -checkcompressed 1 (default) dep $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	file_write dep.txt test
	file_write target.txt test
	cg razip target.txt
	job testcheckcompressed -checkcompressed 0 -deps {dep.txt} -targets target.txt -code {
		file_write result.txt test
		file_write target.txt test
	}
	job_wait
	gridwait
	set result [bsort [glob *]]
	cd $::testdir
	set result
} {dep.txt log.*.finished result.txt target.txt target.txt.rz} match

test job "rmtargets1 $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init -dcleanup never
	file_write data1.txt test1
	job data1 -deps {data1.txt} -rmtargets data1.txt -code {
		after 1000
		file delete data1.txt
	}
	job data2 -optional 1 -deps {data1.txt} -targets data2.txt -code {
		file copy data1.txt data2.txt
	}
	job_wait
	gridwait
	set temp [file_read log_jobs/data2.log]
	set result [list [bsort [glob *]] [regexp {missing dependency "data1.txt"} $temp]]
	cd $::testdir
	set result
} {{log.*.finished log_jobs} 1} match

test job "rmtargets2 $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init -dcleanup never
	job data1 -targets data1.txt -code {
		after 1000
		file_write data1.txt test1
	}
	job data2 -deps {data1.txt} -targets data2.txt -code {
		file copy data1.txt data2.txt
	}
	job rmdata1 -deps {data1.txt data2.txt} -rmtargets data1.txt -code {
		after 1000
		file delete data1.txt
	}
	job data3 -optional 1 -deps {data1.txt} -targets data3.txt -code {
		file copy data1.txt data3.txt
	}
	job_wait
	gridwait
	set temp [file_read log_jobs/data3.log]
	set result [list [bsort [glob *]] [regexp {missing dependency "data1.txt"} $temp]]
	cd $::testdir
	set result
} {{data2.txt log.*.finished log_jobs} 1} match

test job "do not run if deps not done $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	job data1 -targets data1.txt -code {
		after 1000
		file_write data1.txt test1
	}
	job data2 -deps {data1.txt} -targets data2.txt -code {
		# intentional error: target not made
	}
	job rmdata1 -deps {data1.txt data2.txt} -rmtargets data1.txt -code {
		after 1000
		file delete data1.txt
	}
	job_wait
	gridwait
	set result [list [bsort [glob *]]]
	cd $::testdir
	set result
} {{data1.txt log.*.error log_jobs}} match

test job "rmtargets with gzip $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	job makedata -targets data.txt -code {
		after 1000
		file_write data.txt test1
	}
	job compress -checkcompressed 0 -deps {data.txt} -targets data.txt.gz -rmtargets data.txt -code {
		exec gzip data.txt
	}
	job result -deps {data.txt} -targets result.txt -code {
		after 1000
		exec {*}[gzcat $dep] $dep > result.txt
	}
	job_wait
	gridwait
	set result [list [bsort [glob *]] [file_read result.txt]]
	cd $::testdir
	set result
} {{data.txt.gz log.*.finished result.txt} test1} match

test job "rmtargets with gzip exists $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	file_write data.txt testpre
	exec gzip data.txt
	job makedata -targets data.txt -skip data.txt -code {
		after 1000
		file_write data.txt test1
	}
	job compress -optional 1 -checkcompressed 0 -deps {data.txt} -targets data.txt.gz -rmtargets data.txt -code {
		exec gzip data.txt
	}
	job result -deps {data.txt} -targets result.txt -code {
		after 1000
		exec {*}[gzcat $dep] $dep > result.txt
	}
	job_wait
	gridwait
	set result [list [bsort [glob *]] [file_read result.txt]]
	cd $::testdir
	set result
} {{data.txt.gz log.*.finished result.txt} testpre} match

test job "rmtargets and -checkcompressed 0 on previous targets $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	job makedata -targets data.txt -code {
		after 1000
		file_write data.txt test1
	}
	job compress -checkcompressed 0 -deps {data.txt} -targets data.txt.rz -rmtargets data.txt -code {
		cg razip data.txt
	}
	job result -optional 1 -checkcompressed 0 -deps {data.txt} -targets result.txt -code {
		after 1000
		exec {*}[gzcat $dep] $dep > result.txt
	}
	job_wait
	gridwait
	set result [bsort [glob *]]
	cd $::testdir
	set result
} {data.txt.rz log.*.finished} match

test job "rmtargets and -checkcompressed 0 on previous targets write one first $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	file_write data.txt testpre
	cg razip data.txt
	job makedata -targets data.txt -code {
		after 1000
		file_write data.txt test1
	}
	job compress -optional 1 -checkcompressed 0 -deps {data.txt} -targets data.txt.gz -rmtargets data.txt -code {
		cg razip data.txt
	}
	job result -optional 1 -checkcompressed 0 -deps {data.txt} -targets result.txt -code {
		after 1000
		exec {*}[gzcat $dep] $dep > result.txt
	}
	job_wait
	gridwait
	set result [bsort [glob *]]
	cd $::testdir
	set result
} {data.txt.rz log.*.finished} match

test job "rmtargets afterwards with gzip exists $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	file_write data.txt testpre
	exec gzip data.txt
	test_job_init
	job makedata -targets data.txt -code {
		after 1000
		file_write data.txt test1
	}
	job result -deps {data.txt} -targets result.txt -code {
		after 1000
		exec {*}[gzcat $dep] $dep > result.txt
	}
	job compress -optional 1 -checkcompressed 0 -deps {data.txt result.txt} -targets data.txt.gz -rmtargets data.txt -code {
		exec gzip data.txt
	}
	job_wait
	gridwait
	set result [list [bsort [glob *]] [file_read result.txt]]
	cd $::testdir
	set result
} {{data.txt.gz log.*.finished result.txt} testpre} match

test job "jobforce $testname" {
	cd $::testdir/tmp
	test_job_init
	file_write data1.txt test1
	job data1 -targets data1.txt -force 1 -code {
		after 1000
		set c [file_read $target]
		file_write data1.txt $c\ntest2
	}
	job data2 -deps {data1.txt} -targets data2.txt -code {
		set c [file_read $dep]
		file_write $target $c\ntest3
		
	}
	job_wait
	gridwait
	set result [file_read data2.txt]
	cd $::testdir
	set result
} {test1
test2
test3}

test job "jobtargetexists 1 $testname" {
	cd $::testdir/tmp
	test_job_init
	writetestfiles data1.txt data2.txt data3.txt
	set result {}
	jobtargetexists data3.txt {data1.txt data2.txt}
} 1

test job "jobtargetexists 2 $testname" {
	cd $::testdir/tmp
	test_job_init
	writetestfiles data1.txt data2.txt data3.txt
	set result {}
	jobtargetexists {data2.txt data3.txt} {data1.txt}
} 1

test job "jobtargetexists 3 $testname" {
	cd $::testdir/tmp
	test_job_init
	writetestfiles data2.txt data3.txt data1.txt
	set result {}
	jobtargetexists data3.txt {data1.txt data2.txt}
} 0

test job "jobtargetexists 4 $testname" {
	cd $::testdir/tmp
	test_job_init
	writetestfiles data2.txt data3.txt data1.txt
	set result {}
	jobtargetexists {data2.txt data3.txt} {data1.txt}
} 0

test job "jobtargetexists 5 $testname" {
	cd $::testdir/tmp
	test_job_init
	writetestfiles data2.txt data1.txt data3.txt
	set result {}
	jobtargetexists data3.txt {data1.txt data2.txt}
} 1

test job "jobtargetexists 6 $testname" {
	cd $::testdir/tmp
	test_job_init
	writetestfiles data2.txt data1.txt data3.txt
	set result {}
	jobtargetexists {data2.txt data3.txt} {data1.txt}
} 0

test job "jobtargetexists -checkdepexists 0 (default) $testname" {
	cd $::testdir/tmp
	test_job_init
	writetestfiles data1.txt data3.txt
	set result {}
	jobtargetexists data3.txt {data1.txt data2.txt}
} 1

test job "jobtargetexists -checkdepexists 1 $testname" {
	cd $::testdir/tmp
	test_job_init
	writetestfiles data1.txt data3.txt
	set result {}
	jobtargetexists -checkdepexists 1 data3.txt {data1.txt data2.txt}
} 0

test job "basic chain rerun $testname" {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	file_write test1.txt test1\n
	after 10
	job job1 -deps {test1.txt} -targets {test2.txt} -code {
		set c [file_read $dep]
		file_write $target ${c}test2\n
	}
	job job2 -deps {test2.txt} -targets {test3.txt} -code {
		error "wrong this time"
	}
	job_wait
	gridwait
	after 10
	set result {}
	lappend result [bsort [glob *]]
	test_job_init
	job job1 -deps {test1.txt} -targets {test2.txt} -code {
		set c [file_read $dep]
		file_write $target ${c}test2\n
	}
	job job2 -deps {test2.txt} -targets {test3.txt} -code {
		set c [file_read $dep]
		file_write $target ${c}test3\n
	}
	job_wait
	gridwait
	lappend result [bsort [glob *]]
	cd $::testdir
	set result
} {{log.*.error log_jobs test1.txt test2.txt} {log.*.error log.*.finished test1.txt test2.txt test3.txt}} match

test job "job_update -r 1 $testname" {
	test_cleantmp
	cd $::testdir/tmp
	test_job_init
	file_write test1.txt test1\n
	after 10
	job job1 -deps {test1.txt} -targets {test2.txt} -code {
		set c [file_read $dep]
		file_write $target ${c}test2\n
	}
	job job2 -deps {test2.txt} -targets {test3.txt} -code {
		error "wrong this time"
	}
	job_wait
	gridwait
	after 10
	set result {}
	lappend result [bsort [glob *]]
	file delete -force log_jobs
	test_job_init
	job job1 -deps {test1.txt} -targets {test2.txt} -code {
		set c [file_read $dep]
		file_write $target ${c}test2\n
	}
	job job2 -deps {test2.txt} -targets {test3.txt} -code {
		set c [file_read $dep]
		file_write $target ${c}test3\n
	}
	job_wait
	gridwait
	cg job_update -r 1 [lindex [glob log.*.finished] 0]
	lappend result [bsort [glob *]]
	cd $::testdir
	set result
} {{log.*.error log_jobs test1.txt test2.txt} {log.*.finished test1.txt test2.txt test3.txt}} match

test job "-cores 2 $testname" {
	test_cleantmp
	test_job_init
	job test -targets tmp/dep.txt -code {
		file_write tmp/dep.txt test
	}
	job test2 -cores 2 -deps {tmp/dep.txt} -targets tmp/target.txt -code {
		file_write tmp/target.txt [file_read $dep]_ok
	}
	job_wait ; gridwait
	cg select -f cores -rc 1 [glob tmp/log.*]
} {cores
1
2
}

test job "job_cleanup_add $testname" {
	test_job_init
	job test -targets tmp/dep.txt -code {
		file_write tmp/dep.txt test
	}
	job_cleanup_add tmp/dep.txt
	job test2 -deps {tmp/dep.txt} -targets tmp/target.txt -code {
		file_write tmp/target.txt [file_read $dep]_ok
	}
	job_wait ; gridwait
	glob tmp/*.txt
} {tmp/target.txt}

test job "job_cleanup_ifempty_add $testname" {
	test_job_init
	job test -targets {tmp/dep.txt tmp/tmp.dir tmp/tmp2.dir} -code {
		file_write tmp/dep.txt test
		file mkdir tmp/tmp.dir
		file_write tmp/tmp.dir/test.txt test
		file mkdir tmp/tmp2.dir
	}
	job_cleanup_ifempty_add tmp/tmp.dir tmp/tmp2.dir
	job test2 -deps {tmp/dep.txt} -targets tmp/target.txt -code {
		file_write tmp/target.txt [file_read $dep]_ok
	}
	job_wait ; gridwait
	glob tmp/*.dir
} {tmp/tmp.dir}

if {$testname eq "-d sge"} {
	test job ", in name $testname" {
		cd $::testdir
		test_cleantmp
		cd $::testdir/tmp
		test_job_init -skipjoberrors 0
		file_write test1.txt test1\n
		job jobmissing -optional 1 -deps {missing.txt} -targets {test3.txt} -code {
			file_write $target missing\n
		}
		job job1 -deps {test1.txt} -targets {test2.txt} -code {
			set c [file_read $dep]
			file_write $target ${c}test2\n
		}
		job_wait
		gridwait
		set result [list [bsort [glob *]] [file_read test2.txt]]
		cd $::testdir
		set result
	} {Cannot submit job to sge: it has a comma in the output file *job1.out, which grid engine sometimes has problems with} match error
}

test job "filename too long $testname" {
	cd $::testdir/tmp
	test_job_init
	file_write data1.txt test1
	set name [string_fill __abcdefghijklmnopqrstuvwxyz 20]__1234567890
	job data1-$name -targets data1.txt -force 1 -code {
		after 1000
		set c [file_read $target]
		file_write data1.txt $c\ntest2
	}
	job data2-$name -deps {data1.txt} -targets data2.txt -code {
		set c [file_read $dep]
		file_write $target $c\ntest3
		
	}
	job_wait
	gridwait
	set result [file_read data2.txt]
	cd $::testdir
	set result
} {test1
test2
test3}

# end of block
}

# only test in direct
foreach {testname initcode} {
	"direct" {uplevel job_init -skipjoberrors 1 $args}
} break 
proc test_job_init {args} $initcode
proc gridwait {} {}

test job {deps both compressed and uncompressed} {
	cd $::testdir
	test_cleantmp
	cd $::testdir/tmp
	file_write dep.txt test
	exec gzip -c dep.txt > dep.txt.gz
	job_init
	job test -deps {dep.txt} -targets result.txt -code {
		foreach dep $deps {lappend result [file tail $dep]}
		file_write $target $result
	}
	file_read result.txt
} {dep.txt}

test job {gzarraynames} {
	array set a {dep.txt 1 dep.txt.gz 1 x 1}
	gzarraynames a dep.*
} {dep.txt dep.txt.gz}

test job {gzarraynames} {
	array set a {dep.txt 1 dep.txt2 1 x 1}
	lsort [gzarraynames a dep.*]
} {dep.txt dep.txt2}

test job_getinfo {job_getinfo basic} {
	file_write tmp/dep.txt test
	job_init
	job_getinfo 1
	job test -deps {tmp/dep.txt} -targets tmp/target.txt -code {
		puts ok
	}
	foreach {deps targets} [job_getinfo 0] break
	list $deps $targets
} {*tmp/dep.txt *tmp/target.txt} match

test job_getinfo {job_getinfo skipped} {
	file_write tmp/dep.txt test
	after 10
	file_write tmp/target.txt target
	job_init
	job_getinfo 1
	job test -deps {tmp/dep.txt} -targets tmp/target.txt -code {
		puts ok
	}
	foreach {deps targets} [job_getinfo 0] break
	list $deps $targets
} {{} {}}

test cgjob_files {basic} {
	array set cgjob_id {abc.tsv 1 abc.tsvx 1 abc.ts 1}
	cgjob_files cgjob_id a*.tsv
} abc.tsv

test cgjob_files {basic} {
	array set cgjob_id {abc.tsv 1 ab/c.tsv 1}
	cgjob_files cgjob_id a*.tsv
} abc.tsv

test job {dont give error on missing dependency if skiptarget exists} {
	file delete dep.txt
	file_write result.txt test
	job_init
	job test -skip result.txt -deps {dep.txt} -targets preresult.txt -code {
		file_write $target test2
	}
	file_read result.txt
} {test}

test job {inf missing dep and target} {
	job_init
	# logverbose 2
	file delete dep.txt result.txt
	job test -deps {dep.txt} -targets result.txt -code {
		file_write $target test2
	}
	file_read result.txt
} {error trying to run job test:
missing dependency "dep.txt"} error

test sge_safename {shorten prefix} {
	set name test
	set prefix jvar_longshot_map-sminimap2-xxxxx_xxxxx_xx_xxxx_xxxx_v4.0.11__f1071ce_xxxxxx_hg38s.bam.2020-07-25_04-44_560
	sge_safename $name $prefix
} {jvar_longshot_map-sminimap2-xxxxx_xxxxx_xx_xxxx_..11__f1071ce_xxxxxx_hg38s.bam.2020-07-25_04-44_560#test}

set ::env(PATH) $keeppath

cd $keepdir

testsummarize

