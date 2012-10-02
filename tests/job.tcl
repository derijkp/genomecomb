#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

catch {file delete -force {*}[glob tmp/*]}

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

proc jobtest {args} {
	set args [job_args $args]
	foreach {srcdir destdir header} $args break
	set srcdir [file normalize $srcdir]
	set destdir [file normalize $destdir]
	job_logdir $destdir/log_jobs
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
		set targets
	}
	job sumpattern -deps {$srcdir/cgdat*.tsv} \
	-targets {$destdir/sumpattern.log} \
	-ptargets {$destdir/sumpattern-(.*).txt} \
	-vars destdir -code {
		# var target contains target
		# var target1 contains first braced part of target (= destdir)
		# var target1 contains second braced part of target (= name)
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
			file rename $target.temp $target
			lappend targets $target
		}
		close $f
		file_write $destdir/sumpattern.log [join $targets \n]
	}
	job test -foreach {^$srcdir/cgdat(.*)\.tsv$} -targets {$destdir/test.txt} -code {
		exec wc $dep > $target
	}
	job sum2 -foreach {^$destdir/sum-(.*)\.txt$} \
	  -skip {$destdir/all2.txt} \
	  -targets {$destdir/sum2-\1.txt} -code {
		for {set i 1} {$i < 5} {incr i} {
			puts "progress $i"
			after 250
		}
		file copy $dep $target.temp
		exec echo 2 >> $target.temp
		file rename $target.temp $target
	}
	job error_all.txt -deps {$srcdir/notpresent.txt} -targets {$destdir/all.txt} -code {
		error "This should not be executed, as the dependencies are not fullfilled, the other target is used"
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


test job {basic} {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
	cd $::testdir/tmp
	jobtest ../data test testh
	set result [list \
		[lsort -dict [glob test/*]] \
		[file_read test/all.txt] \
		[file_read test/sum2-test3.txt] \
	]
	cd $::testdir
	set result
} {{test/all.txt test/all2.txt test/log_jobs test/sum-test1.txt test/sum-test2.txt test/sum-test3.txt test/sum2-test1.txt test/sum2-test2.txt test/sum2-test3.txt test/sumpattern-test1.txt test/sumpattern-test2.txt test/sumpattern-test3.txt test/sumpattern.log test/test.txt} {testh
1+2=3
3+4+5=12
6+7+8=21
} {6+7+8=21
2
}}

test job {--force 0} {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
	cd $::testdir/tmp
	file mkdir test
	file_write test/all.txt error
	jobtest --force 0 ../data test testh
	set result [list \
		[lsort -dict [glob test/*]] \
		[file_read test/all.txt] \
		[file_read test/sum2-test3.txt] \
	]
	cd $::testdir
	set result
} {{test/all.txt test/all2.txt test/log_jobs test/sum-test1.txt test/sum-test2.txt test/sum-test3.txt test/sum2-test1.txt test/sum2-test2.txt test/sum2-test3.txt test/sumpattern-test1.txt test/sumpattern-test2.txt test/sumpattern-test3.txt test/sumpattern.log test/test.txt} error {6+7+8=21
2
}}

test job {--force 1} {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
	cd $::testdir/tmp
	file mkdir test
	file_write test/all.txt error
	jobtest --force 1 ../data test testh
	set result [list \
		[lsort -dict [glob test/*]] \
		[file_read test/all.txt.old1] \
		[file_read test/all.txt] \
		[file_read test/sum2-test3.txt] \
	]
	cd $::testdir
	set result
} {{test/all.txt test/all.txt.old1 test/all2.txt test/log_jobs test/sum-test1.txt test/sum-test2.txt test/sum-test3.txt test/sum2-test1.txt test/sum2-test2.txt test/sum2-test3.txt test/sumpattern-test1.txt test/sumpattern-test2.txt test/sumpattern-test3.txt test/sumpattern.log test/test.txt} error {testh
1+2=3
3+4+5=12
6+7+8=21
} {6+7+8=21
2
}}

set ::env(PATH) $keeppath

testsummarize

