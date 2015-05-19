#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

# use these for trying out individual tests
set testname "-d direct"
proc test_job_init {} {uplevel job_init -silent -job_skiperrors}
proc gridwait {} {}
if 0 {
	set testname "-d direct"
	proc test_job_init {} {uplevel job_init -silent -job_skiperrors}
	proc gridwait {} {}

	set testname "-d 30"
	proc test_job_init {} {uplevel job_init -silent -d 30}
	proc gridwait {} {}

	set testname "-d sge"
	proc test_job_init {} {uplevel job_init -silent -d sge}
	proc gridwait {} {
		while 1 {
			after 500
			puts -nonewline .
			flush stdout
			if {[exec qstat] eq ""} break
		}
	}
}

catch {file delete -force {*}[glob tmp/*]}

proc jobtest {args} {
	set args [job_args $args]
	foreach {srcdir destdir header} $args break
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
	job sumvar -deps {$srcdir/cgdat*.tsv} \
	-targets {$destdir/sum-$_names.txt} \
	-vars destdir -code {
		for {set i 1} {$i < 5} {incr i} {
			puts "progress $target $i"
			after 250
		}
		set targets {}
		set f [open [glob $dep]]
		while {![eof $f]} {
			set line [split [gets $f] \t]
			if {$line eq ""} continue
			puts line=$line
			set name [lindex $line 0]
			set target $destdir/sum-$name.txt
			set o [open $target.temp w]
			set calc [join [lrange $line 1 end] +]
			puts $o $calc=[expr $calc]
			close $o
			file rename -force $target.temp $target
			lappend targets $target
		}
		close $f
	}
	job sumpattern -deps {$srcdir/cgdat*.tsv} \
	-targets {$destdir/sumpattern.log} \
	-ptargets {$destdir/sumpattern-*.txt} \
	-vars destdir -code {
		# var target contains target
		# var target1 contains first braced part of target (= destdir)
		# var target2 contains second braced part of target (= name)
		# var deps contains a list of all dependencies
		# var dep contains the first element of the first dependency (so you do not have to do [lindex $deps 0] to get it)
		for {set i 1} {$i < 5} {incr i} {
			puts "progress $target $i"
			after 250
		}
		set targets {}
		set f [open [glob $dep]]
		while {![eof $f]} {
			set line [split [gets $f] \t]
			if {$line eq ""} continue
			puts line=$line
			set name [lindex $line 0]
			set target $destdir/sumpattern-$name.txt
			set o [open $target.temp w]
			set calc [join [lrange $line 1 end] +]
			puts $o "p $calc=[expr $calc]"
			close $o
			file rename -force $target.temp $target
			lappend targets $target
		}
		close $f
		file_write $destdir/sumpattern.log [join $targets \n]
	}
	job sumpattern2 -foreach {^$destdir/sumpattern-(.*)\.txt$} \
	  -skip {$destdir/allp2.txt} \
	  -cores 2 \
	  -targets {$destdir/sumpattern2-\1.txt} -code {
		for {set i 1} {$i < 5} {incr i} {
			puts "progress $i"
			after 250
		}
		file copy $dep $target.temp
		exec echo 2 >> $target.temp
		file rename -force $target.temp $target
	}
	job allp2.txt -vars header -deps {^$destdir/sumpattern2-(.*)\.txt$} -targets {$destdir/allp2.txt} -code {
		after 500
		file_write $target.temp $header\n
		exec cat {*}$deps >> $target.temp
		file rename -force $target.temp $target
	}
	job allp.txt -vars header -deps {^$destdir/sumpattern-(.*)\.txt$} -targets {$destdir/allp.txt} -code {
		file_write $target.temp $header\n
		exec cat {*}$deps >> $target.temp
		file rename -force $target.temp $target
	}
	job test -foreach {^$srcdir/cgdat(.*)\.tsv$} -targets {$destdir/test.txt} -code {
		exec wc $dep > $target
	}
	job sum2 -foreach {^$destdir/sum-(.*)\.txt$} \
	  -skip {$destdir/all2.txt} \
	  -targets {$destdir/sum2-\1.txt} -code {
		for {set i 1} {$i < 5} {incr i} {
			puts stderr "progress $i"
			after 250
		}
		file copy $dep $target.temp
		exec echo 2 >> $target.temp
		file rename -force $target.temp $target
	}
	job error_all.txt -deps {$srcdir/notpresent.txt} -targets {$destdir/all.txt} -code {
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
		exec cat {*}$deps >> $target.temp
		file rename -force $target.temp $target
	}
	set keepdir [pwd]
	cd $destdir
	job all2.txt -vars header -deps {^sum2-(.*)\.txt$} -targets {all2.txt} -code {
		file_write $target.temp $header\n
		exec cat {*}$deps >> $target.temp
		file rename -force $target.temp $target
	}
	cd $keepdir
}

proc jobtestnojobs {args} {
	set args [job_args $args]
	foreach destdir $args break
	job_logdir $destdir/log_jobs
	job nodeps1 -deps {abcd} -targets {abcde} -code {
	}
	job nodeps2 -deps {abcde} -targets {abcdef} -code {
	}
}

proc jobtestlong {args} {
	set args [job_args $args]
	foreach destdir $args break
	file mkdir $destdir
	cd $destdir
	job_logdir $destdir/log_jobs
	job long -deps {} -targets {long.txt} -code {
		for {set i 0} {$i < 4} {incr i} {
			after 1000
			puts $i
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

# ------------------------------------
# test in different "processing modes"
# ------------------------------------
foreach {testname initcode} {
	"direct" {uplevel job_init -silent -job_skiperrors}
	"-d 4" {uplevel job_init -silent -d 4}
	"-d 30" {uplevel job_init -silent -d 30}
	"-d sge" {uplevel job_init -silent -d sge}
} {
# start of block

if {$testname eq "-d sge"} {
	if {[catch {exec qstat}]} {
		puts "Cannot test sge option (missing qstat; grid engine not installed?)"
		continue
	}
	proc gridwait {} {
		while 1 {
			after 500
			puts -nonewline .
			flush stdout
			if {[exec qstat] eq ""} break
		}
	}
} else {
	proc gridwait {} {}
}

proc test_job_init {} $initcode

test job "basic chain $testname" {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
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
	set result [list [lsort -dict [glob *]] [file_read test3.txt]]
	cd $::testdir
	set result
} {{log_jobs test1.txt test2.txt test3.txt} {test1
test2
test3
}}

test job "foreach $testname" {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
	cd $::testdir/tmp
	test_job_init
	file_write test1.txt test1\n
	file_write test2.txt test2\n
	job job1 -foreach {^test(.*)\.txt$} -targets {rtest\1.txt} -code {
		set c [file_read $dep]
		file_write $target r${c}
	}
	job_wait
	gridwait
	set result [list [lsort -dict [glob *]] [file_read rtest1.txt]]
	cd $::testdir
	set result
} {{log_jobs rtest1.txt rtest2.txt test1.txt test2.txt} {rtest1
}}

test job "chained foreach $testname" {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
	cd $::testdir/tmp
	test_job_init
	job init -deps {} -targets {test1.txt test2.txt} -code {
		file_write test1.txt test1\n
		file_write test2.txt test2\n
	}
	job job1 -foreach {^test(.*)\.txt$} -targets {rtest\1.txt} -code {
		set c [file_read $dep]
		file_write $target r${c}
	}
	job job2 -deps {rtest1.txt} -targets {final1.txt} -code {
		set c [file_read $dep]
		file_write $target f${c}
	}
	job_wait
	gridwait
	set result [list [lsort -dict [glob *]] [file_read final1.txt]]
	cd $::testdir
	set result
} {{final1.txt log_jobs rtest1.txt rtest2.txt test1.txt test2.txt} {frtest1
}}

test job "basic chain --force 0 $testname" {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
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
	set result [list [lsort -dict [glob *]] [file_read test3.txt]]
	cd $::testdir
	set result
} {{log_jobs test1.txt test2.txt test3.txt} {error3
}}

test job "basic chain --force 1 $testname" {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
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
	set result [list [lsort -dict [glob *]] [file_read test2.txt] [file_read test3.txt]]
	cd $::testdir
	set result
} {{log_jobs test1.txt test2.txt test3.txt} {test1
test2
} {test1
test2
test3
}}

test job "time chain $testname" {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
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
	set result [list [lsort -dict [glob *]] [file_read test2.txt] [file_read test3.txt] [file_read test3.txt.old]]
	cd $::testdir
	set result
} {{log_jobs test1.txt test2.txt test3.txt test3.txt.old} {error2
} {error2
test3
} {error3
}}

test job "basic $testname" {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
	cd $::testdir/tmp
	test_job_init
	jobtest ../data test testh
	job_wait
	gridwait
	set result [list \
		[lsort -dict [glob test/*]] \
		[glob test/log_jobs/all.txt.log] \
		[file_read test/all.txt] \
		[file_read test/sum2-test3.txt] \
		[string range [file_read test/log_jobs/joberror.err] 0 20] \
	]
	cd $::testdir
	set result
} {{test/all.txt test/all2.txt test/allp.txt test/allp2.txt test/log_jobs test/sum-test1.txt test/sum-test2.txt test/sum-test3.txt test/sum2-test1.txt test/sum2-test2.txt test/sum2-test3.txt test/sumpattern-test1.txt test/sumpattern-test2.txt test/sumpattern-test3.txt test/sumpattern.log test/sumpattern2-test1.txt test/sumpattern2-test2.txt test/sumpattern2-test3.txt test/test.txt} test/log_jobs/all.txt.log {testh
1+2=3
3+4+5=12
6+7+8=21
} {6+7+8=21
2
} {Intentional job error}}

test job "basic status $testname" {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
	cd $::testdir/tmp
	test_job_init
	jobtest ../data test testh
	job_wait
	gridwait
	set result [list \
		[lsort -dict [glob test/*]] \
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
	jobtest ../data test testh
	job_wait
	glob jobdeps.dot jobdeps.ps
} {jobdeps.dot jobdeps.ps}

test job "--force 0 $testname" {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
	cd $::testdir/tmp
	file mkdir test
	job_init -silent -job_skiperrors
	jobtest --force 0 ../data test testh
	job_wait
	gridwait
	after 1000
	file_write test/all.txt error
	test_job_init
	jobtest --force 0 ../data test testh
	job_wait
	gridwait
	set result [list \
		[lsort -dict [glob test/*]] \
		[file_read test/all.txt] \
		[file_read test/sum2-test3.txt] \
	]
	cd $::testdir
	set result
} {{test/all.txt test/all2.txt test/allp.txt test/allp2.txt test/log_jobs test/sum-test1.txt test/sum-test2.txt test/sum-test3.txt test/sum2-test1.txt test/sum2-test2.txt test/sum2-test3.txt test/sumpattern-test1.txt test/sumpattern-test2.txt test/sumpattern-test3.txt test/sumpattern.log test/sumpattern2-test1.txt test/sumpattern2-test2.txt test/sumpattern2-test3.txt test/test.txt} error {6+7+8=21
2
}}

test job "--force 1 $testname" {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
	cd $::testdir/tmp
	file mkdir test
	job_init -silent -job_skiperrors
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
		[lsort -dict [glob test/*]] \
		[file_read test/all.txt] \
		[file_read test/sum2-test3.txt] \
	]
	cd $::testdir
	set result
} {{test/all.txt test/all2.txt test/allp.txt test/allp2.txt test/log_jobs test/sum-test1.txt test/sum-test2.txt test/sum-test3.txt test/sum2-test1.txt test/sum2-test2.txt test/sum2-test3.txt test/sumpattern-test1.txt test/sumpattern-test2.txt test/sumpattern-test3.txt test/sumpattern.log test/sumpattern2-test1.txt test/sumpattern2-test2.txt test/sumpattern2-test3.txt test/test.txt} {testh
1+2=3
3+4+5=12
6+7+8=21
} {6+7+8=21
2
}}

test job "time $testname" {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
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
		[lsort -dict [glob test/*]] \
		[file_read test/all.txt] \
		[file_read test/sum2-test3.txt] \
	]
	cd $::testdir
	set result
} {{test/all.txt test/all.txt.old test/all2.txt test/all2.txt.old test/allp.txt test/allp2.txt test/log_jobs test/sum-test1.txt test/sum-test2.txt test/sum-test3.txt test/sum2-test1.txt test/sum2-test2.txt test/sum2-test3.txt test/sum2-test3.txt.old test/sumpattern-test1.txt test/sumpattern-test2.txt test/sumpattern-test3.txt test/sumpattern.log test/sumpattern2-test1.txt test/sumpattern2-test2.txt test/sumpattern2-test3.txt test/test.txt} {testh
1+2=3
3+4+5=12
replaced
} {replaced
2
}}

test job "-skip: not present $testname" {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
	cd $::testdir/tmp
	test_job_init
	job testskip -deps {} -targets result.txt -skip {skip1.txt skip2.txt} -code {
		file_write result.txt test
	}
	job_wait
	gridwait
	set result [lsort -dict [glob *]]
	cd $::testdir
	set result
} {log_jobs result.txt}

test job "-skip: only one present $testname" {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
	cd $::testdir/tmp
	test_job_init
	file_write skip1.txt test1
	job testskip -deps {} -targets result.txt -skip {skip1.txt skip2.txt} -code {
		file_write result.txt test
	}
	job_wait
	gridwait
	set result [lsort -dict [glob *]]
	cd $::testdir
	set result
} {log_jobs result.txt skip1.txt}

test job "-skip: all present $testname" {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
	cd $::testdir/tmp
	test_job_init
	file_write skip1.txt test1
	file_write skip2.txt test2
	job testskip -deps {} -targets result.txt -skip {skip1.txt skip2.txt} -code {
		file_write result.txt test
	}
	job_wait
	gridwait
	set result [lsort -dict [glob *]]
	cd $::testdir
	set result
} {log_jobs skip1.txt skip2.txt}

test job "-skip -skip: none present $testname" {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
	cd $::testdir/tmp
	test_job_init
	job testskip -deps {} -targets result.txt -skip skip1.txt -skip skip2.txt -code {
		file_write result.txt test
	}
	job_wait
	gridwait
	set result [lsort -dict [glob *]]
	cd $::testdir
	set result
} {log_jobs result.txt}

test job "-skip -skip: one present $testname" {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
	cd $::testdir/tmp
	test_job_init
	file_write skip2.txt test2
	job testskip -deps {} -targets result.txt -skip skip1.txt -skip skip2.txt -code {
		file_write result.txt test
	}
	job_wait
	gridwait
	set result [lsort -dict [glob *]]
	cd $::testdir
	set result
} {log_jobs skip2.txt}

test job {jobtestnojobs} {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
	cd $::testdir/tmp
	test_job_init
	jobtestnojobs $::testdir/tmp
	job_wait
	gridwait
} {}

test job {jobtestlong} {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
	cd $::testdir/tmp
	test_job_init
	jobtestlong $::testdir/tmp
	job_wait
	gridwait
} {}

test job {no -targets} {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
	cd $::testdir/tmp
	test_job_init
	file_write dep.txt dep
	job testnotarget -deps {dep.txt} -code {
		file_write result.txt test
	}
	job_wait
	gridwait
	set result [lsort -dict [glob *]]
	cd $::testdir
	set result
} {dep.txt log_jobs result.txt}

test job {no -targets, dep not found} {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
	cd $::testdir/tmp
	test_job_init
	job testnotarget -deps {dep.txt} -code {
		file_write result.txt test
	}
	job_wait
	gridwait
	set result [lsort -dict [glob *]]
	cd $::testdir
	set result
} {log_jobs}

test job {no -checkcompressed 1 (default), dep} {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
	cd $::testdir/tmp
	test_job_init
	file_write dep.txt test
	cg razip dep.txt
	job testcheckcompressed -deps {dep.txt} -targets result.txt -code {
		file_write result.txt test
	}
	job_wait
	gridwait
	set result [lsort -dict [glob *]]
	cd $::testdir
	set result
} {dep.txt.rz log_jobs result.txt}

test job {no -checkcompressed 1 (default), dep} {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
	cd $::testdir/tmp
	test_job_init
	file_write dep.txt test
	cg razip dep.txt
	job testcheckcompressed -checkcompressed 0 -deps {dep.txt} -targets result.txt -code {
		file_write result.txt test
	}
	job_wait
	gridwait
	set result [lsort -dict [glob *]]
	cd $::testdir
	set result
} {dep.txt.rz log_jobs}

test job {no -checkcompressed 1 (default), targets} {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
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
	set result [lsort -dict [glob *]]
	cd $::testdir
	set result
} {dep.txt log_jobs target.txt.rz}

test job {no -checkcompressed 1 (default), dep} {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
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
	set result [lsort -dict [glob *]]
	cd $::testdir
	set result
} {dep.txt log_jobs result.txt target.txt target.txt.rz}

test job "rmtargets1 $testname" {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
	cd $::testdir/tmp
	test_job_init
	file_write data1.txt test1
	job data1 -deps {data1.txt} -rmtargets data1.txt -code {
		after 1000
		file delete data1.txt
	}
	job data2 -deps {data1.txt} -targets data2.txt -code {
		file copy data1.txt data2.txt
	}
	job_wait
	gridwait
	set temp [file_read log_jobs/data2.log]
	set result [list [lsort -dict [glob *]] [regexp {missing dependency data1.txt} $temp]]
	cd $::testdir
	set result
} {log_jobs 1}

test job "rmtargets2 $testname" {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
	cd $::testdir/tmp
	test_job_init
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
	job data3 -deps {data1.txt} -targets data3.txt -code {
		file copy data1.txt data3.txt
	}
	job_wait
	gridwait
	set temp [file_read log_jobs/data3.log]
	set result [list [lsort -dict [glob *]] [regexp {missing dependency data1.txt} $temp]]
	cd $::testdir
	set result
} {{data2.txt log_jobs} 1}

test job "do not run if deps not done $testname" {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
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
	set result [list [lsort -dict [glob *]]]
	cd $::testdir
	set result
} {{data1.txt log_jobs}}

test job "rmtargets with gzip $testname" {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
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
	set result [list [lsort -dict [glob *]] [file_read result.txt]]
	cd $::testdir
	set result
} {{data.txt.gz log_jobs result.txt} test1}

test job "rmtargets with gzip exists $testname" {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
	cd $::testdir/tmp
	test_job_init
	file_write data.txt testpre
	exec gzip data.txt
	job makedata -targets data.txt -skip data.txt -code {
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
	set result [list [lsort -dict [glob *]] [file_read result.txt]]
	cd $::testdir
	set result
} {{data.txt.gz log_jobs result.txt} testpre}

test job "rmtargets and -checkcompressed 0 on previous targets $testname" {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
	cd $::testdir/tmp
	test_job_init
	job makedata -targets data.txt -code {
		after 1000
		file_write data.txt test1
	}
	job compress -checkcompressed 0 -deps {data.txt} -targets data.txt.gz -rmtargets data.txt -code {
		cg_razip data.txt
	}
	job result -checkcompressed 0 -deps {data.txt} -targets result.txt -code {
		after 1000
		exec {*}[gzcat $dep] $dep > result.txt
	}
	job_wait
	gridwait
	set result [lsort -dict [glob *]]
	cd $::testdir
	set result
} {data.txt.rz log_jobs}

test job "rmtargets and -checkcompressed 0 on previous targets, write one first $testname" {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
	cd $::testdir/tmp
	test_job_init
	file_write data.txt testpre
	cg_razip data.txt
	job makedata -targets data.txt -code {
		after 1000
		file_write data.txt test1
	}
	job compress -checkcompressed 0 -deps {data.txt} -targets data.txt.gz -rmtargets data.txt -code {
		cg_razip data.txt
	}
	job result -checkcompressed 0 -deps {data.txt} -targets result.txt -code {
		after 1000
		exec {*}[gzcat $dep] $dep > result.txt
	}
	job_wait
	gridwait
	set result [lsort -dict [glob *]]
	cd $::testdir
	set result
} {data.txt.rz log_jobs}

test job "rmtargets afterwards with gzip exists $testname" {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
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
	job compress -checkcompressed 0 -deps {data.txt result.txt} -targets data.txt.gz -rmtargets data.txt -code {
		exec gzip data.txt
	}
	job_wait
	gridwait
	set result [list [lsort -dict [glob *]] [file_read result.txt]]
	cd $::testdir
	set result
} {{data.txt.gz log_jobs result.txt} testpre}

# end of block
}

# only test in direct
foreach {testname initcode} {
	"direct" {uplevel job_init -silent -job_skiperrors}
} break 

test job {deps both compressed and uncompressed} {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
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
} {dep.txt}

test job {gzarraynames} {
	array set a {dep.txt 1 dep.txt2 1 x 1}
	lsort [gzarraynames a dep.*]
} {dep.txt dep.txt2}

set ::env(PATH) $keeppath

testsummarize

