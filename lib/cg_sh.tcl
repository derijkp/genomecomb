#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

package require Extral

proc cg_wish {args} {
	package require Tk
	set tk 1
	package require Tclx
	signal -restart error SIGINT
	if {[info commands "console"] == "console"} {
		console show
	} else {
		package require ClassyTk
		Classy::cmd
	}
}

proc cg_sh {args} {
	if {[lsearch $args tk] != -1} {
		cg_wish {*}$args
		return
	}
	package require Tclx
	signal -restart error SIGINT
	if {[info commands "console"] == "console"} {
		console show
	} else {
		uplevel #0 {commandloop -prompt1 {puts -nonewline "% "} -prompt2 {puts -nonewline ""}}
	}
}

proc cg_exec {file args} {
	set ::argv $args
	uplevel #0 source $file
}

if {[info exists argv0] && [file tail [info script]] eq [file tail $argv0]} {
	set scriptname [info script]
	package require pkgtools
	set appdir [file dir [pkgtools::startdir]]
	lappend auto_path $appdir/lib
	append env(PATH) :[file dir [file dir $appdir]]/bin:$appdir/bin
	package require Extral
	set ::base $scriptname
	cg_sh {*}$argv
}

