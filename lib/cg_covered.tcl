# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

package require Extral

proc cg_covered args {
	cg_options covered args {
		-n - -namecol {
			set namecol $value
		}
	} {regfile} 0 1
	if {[info exists regfile]} {
		set f [gzopen $regfile]
	} else {
		set f stdin
	}
	set header [tsv_open $f]
	catch {tsv_basicfields $header 3 0} poss
	if {[info exists namecol]} {
		set pos [lsearch $header $namecol]
		lset poss 0 $pos
	}
	if {[lsearch $poss -1] != -1} {
		exiterror "header error: some fields (or alternatives) not found"
	}
	chanexec $f stdout "covered $poss"
}

if {[info exists argv0] && [file tail [info script]] eq [file tail $argv0]} {
	package require pkgtools
	set appdir [file dir [pkgtools::startdir]]
	lappend auto_path $appdir/lib
	append env(PATH) :[file dir [file dir $appdir]]/bin:$appdir/bin
	package require Extral
	set ::base [file tail [info script]]
	cg_covered {*}$argv
}
