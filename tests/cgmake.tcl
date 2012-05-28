#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

catch {file delete -force {*}[glob tmp/*]}

test cgmake {basic} {
	catch {file delete -force {*}[glob tmp/*]}
	cd $::testdir/tmp
	cg exec ../cgmaketest ../data test
	set result [list [lsort -dict [glob *]] [file_read all.txt] [file_read sum2-test3.txt]]
	cd $::testdir
	set result
} {{all.txt sum-test1.txt sum-test2.txt sum-test3.txt sum2-test1.txt sum2-test2.txt sum2-test3.txt} {test
1+2=3
3+4+5=12
6+7+8=21
} {6+7+8=21
2
}}

set ::env(PATH) $keeppath

testsummarize

