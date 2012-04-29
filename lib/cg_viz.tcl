#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg_viz {args} {
	package require Tk
	uplevel #0 [list set argv $args]
	uplevel #0 source $::appdir/cg_viz/cg_viz.tcl
}
