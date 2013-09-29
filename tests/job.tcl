#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

# use these for trying out individual tests
set testname direct
proc test_job_init {} {uplevel job_init}
proc gridwait {} {}
if 0 {
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
	set srcdir [file normalize $srcdir]
	set destdir [file normalize $destdir]
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
			file rename $target.temp $target
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
		file rename $target.temp $target
	}
	job allp2.txt -vars header -deps {^$destdir/sumpattern2-(.*)\.txt$} -targets {$destdir/allp2.txt} -code {
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
		file rename $target.temp $target
	}
	job error_all.txt -deps {$srcdir/notpresent.txt} -targets {$destdir/all.txt} -code {
		error "This should not be executed, as the dependencies are not fullfilled, the other target is used"
	}
	job joberror -deps {^$srcdir/cgdata\.tsv$} -targets {$destdir/joberror.txt} -code {
		# this (intentionally) causes an error, and should by identified by an error in the logs
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
	"direct" {uplevel job_init -silent}
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
	job_init -silent
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
	job_init -silent
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
} {{test/all.txt test/all.txt.old test/all2.txt test/all2.txt.old test/allp.txt test/allp.txt.old test/allp2.txt test/allp2.txt.old test/log_jobs test/sum-test1.txt test/sum-test2.txt test/sum-test3.txt test/sum2-test1.txt test/sum2-test1.txt.old test/sum2-test2.txt test/sum2-test2.txt.old test/sum2-test3.txt test/sum2-test3.txt.old test/sumpattern-test1.txt test/sumpattern-test2.txt test/sumpattern-test3.txt test/sumpattern.log test/sumpattern2-test1.txt test/sumpattern2-test1.txt.old test/sumpattern2-test2.txt test/sumpattern2-test2.txt.old test/sumpattern2-test3.txt test/sumpattern2-test3.txt.old test/test.txt} {testh
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

test job "rmtargets $testname" {
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

test job "rmtargets $testname" {
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

# end of block
}

set ::env(PATH) $keeppath

testsummarize

