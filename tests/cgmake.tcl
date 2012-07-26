#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

catch {file delete -force {*}[glob tmp/*]}

test cgmake {cgmake_expandvars} {
	unset -nocomplain ::a; unset -nocomplain ::b; unset -nocomplain ::c; unset -nocomplain ::d;
	set string {$a($b($c))}
	set ::a(B(C)) A(B(C))
	set ::b(C) B(C)
	set ::c C
	set ::d D
	cgmake_expandvars $string
} {A(B(C))}

test cgmake {cgmake_expandvars} {
	unset -nocomplain ::a; unset -nocomplain ::b; unset -nocomplain ::c; unset -nocomplain ::d;
	set string {$a($b($_c))}
	set ::a(B(C1)) A(B(C1))
	set ::a(B(C2)) A(B(C2))
	set ::b(C1) B(C1)
	set ::b(C2) B(C2)
	set ::c {C1 C2}
	set ::d D
	cgmake_expandvars $string
} {A(B(C1)) A(B(C2))}

test cgmake {cgmake_expandvars} {
	unset -nocomplain ::a; unset -nocomplain ::b; unset -nocomplain ::c; unset -nocomplain ::d;
	set string {$a($_b($_c))}
	set ::a(B1(C1)) A(B1(C1))
	set ::a(B2(C1)) A(B2(C1))
	set ::a(B(C2)) A(B(C2))
	set ::b(C1) {B1(C1) B2(C1)}
	set ::b(C2) B(C2)
	set ::c {C1 C2}
	set ::d D
	cgmake_expandvars $string
} {A(B1(C1)) A(B2(C1)) A(B(C2))}

test cgmake {basic} {
	cd $::testdir
	catch {file delete -force {*}[glob tmp/*]}
	cd $::testdir/tmp
	source ../cgmaketest
	cg_maketest ../data test testh
	set result [list [lsort -dict [glob test/*]] [file_read test/all.txt] [file_read test/sum2-test3.txt]]
	cd $::testdir
	set result
} {{test/all.txt test/sum-test1.txt test/sum-test2.txt test/sum-test3.txt test/sum2-test1.txt test/sum2-test2.txt test/sum2-test3.txt} {testh
1+2=3
3+4+5=12
6+7+8=21
} {6+7+8=21
2
}}

set ::env(PATH) $keeppath

testsummarize

