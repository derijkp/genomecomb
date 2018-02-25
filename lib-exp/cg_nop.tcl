#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

proc cg_nop args {
}

if {[info exists argv0] && [file tail [info script]] eq [file tail $argv0]} {
	package require pkgtools
	set appdir [file dir [pkgtools::startdir]]
	lappend auto_path $appdir/lib
	append env(PATH) [pathsep][file dir [file dir $appdir]]/bin[pathsep]$appdir/bin
	package require Extral
	set ::base [file tail [info script]]
	cg_nop {*}$argv
}
