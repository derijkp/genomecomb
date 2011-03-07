#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

package require Extral

proc cg_sh {args} {
	if {[lsearch $args tk] != -1} {
		package require TK
	}
	package require Tclx
	signal -restart error SIGINT
	if {[info commands "console"] == "console"} {
		console show
	} else {
		uplevel #0 {commandloop -prompt1 {puts -nonewline "% "} -prompt2 {puts -nonewline ""}}
	}
}

if {[info exists argv]} {
	set scriptname [info script]
	package require pkgtools
	set appdir [file dir [pkgtools::startdir]]
	lappend auto_path $appdir/lib
	append env(PATH) :[file dir [file dir $appdir]]/bin:$appdir/bin
	package require Extral
	set ::base $scriptname
	cg_sh {*}$argv
}

