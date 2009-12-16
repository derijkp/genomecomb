#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

lappend auto_path /home/peter/bin/tcl
lappend auto_path /home/peter/dev/completegenomics/lib

set object .g

graphwidget .g
pack .g -fill both -expand yes

set file /complgen/sv/chr20.smoothed.100
set file /complgen/sv/chr1.smoothed.100
set file /complgen/sv/sv79-20s.sv
set file /complgen/sv/sv70-20s.sv
if {[llength [get argv ""]]} {
	set file [lindex $argv 0]
}
.g open $file

