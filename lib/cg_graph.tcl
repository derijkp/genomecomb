#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg_graph {args} {
	global graphd
	package require Tk
	package require ClassyTk
	set object .g
	graphwidget .g
	pack .g -fill both -expand yes
	if {[lindex $args 0] eq "-line"} {
		set opts [list line [lindex $args 1]]
		set args [lrange $args 2 end]
	}
	foreach file $args {
		if {![file exists $file]} {puts stderr "$file does not exist"}
		$object defsettings $file {*}$opts
		.g open $file
	}
}
