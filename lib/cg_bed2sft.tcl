#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

package require Extral

proc cg_bed2sft {args} {
	if {([llength $args] < 0) || ([llength $args] > 2)} {
		errorformat bed2sft
		exit 1
	}
	if {[llength $args] > 0} {
		set filename [lindex $args 0]
		set f [gzopen $filename]
	} else {
		set f stdin
	}
	if {[llength $args] > 1} {
		set outfile [lindex $args 1]
		set o [open $outfile w]
	} else {
		set o stdout
	}
	while {![eof $f]} {
		set line [gets $f]
		if {[inlist {browser track} [lindex $line 0]]} {
			puts $o #$line
		} else {
			break
		}
	}
	set bedfields {chromosome begin end name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts}
	set size [expr {[llength $line]-1}]
	puts $o [join [lrange $bedfields 0 $size] \t]
	puts $o $line
	fconfigure $f -translation binary
	fconfigure $o -translation binary
	fcopy $f $o
	if {$o ne "stdout"} {catch {close $o}}
	if {$f ne "stdin"} {catch {close $f}}
}

if {[info exists argv0] && [file tail [info script]] eq [file tail $argv0]} {
	package require pkgtools
	set appdir [file dir [pkgtools::startdir]]
	lappend auto_path $appdir/lib
	append env(PATH) :[file dir [file dir $appdir]]/bin:$appdir/bin
	package require Extral
	set ::base $scriptname
	cg_clc2sft {*}$argv
}
