#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

package require Extral

proc cg_clc2sft {args} {
	set coveragecutoff 0
	set pos 0
	set minfreq 0.25
	foreach {key value} $args {
		switch -- $key {
			-coverage {
				set coveragecutoff $value
			}
			-minfreq {
				set minfreq $value
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]	
	if {([llength $args] < 0) || ([llength $args] > 2)} {
		errorformat clc2sft
		exit 1
	}
	catch {gzclose $f} ; catch {close $o}
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
	puts $o [join [list chromosome begin end type ref alt status freq alleleSeq1 alleleSeq2 sequenced coverage other otherfreq countA countC countG countT countN countD] \t]
	set line [gets $f]
	set num 0
	array set basica {A 1 T 1 G 1 C 1 N 1 - 1}
	while {![eof $f]} {
		if {$line ne ""}	{error "format error at line $num: expected empty line"}
		incr num
		set line [gets $f]
		if {$line eq "  Zero coverage regions:"} {
			incr num
			set line [gets $f]
			if {$line ne ""}	{error "format error at line $num: expected empty line"}
			while {![eof $f]} {
				incr num
				set line [gets $f]
				if {$line eq ""} break
			}	
		} else {
			set chrom [string trimright [lindex $line 0] :]
			incr num
			set line [gets $f]
			if {$line ne ""}	{error "format error at line $num: expected empty line"}
			while {![eof $f]} {
				incr num
				set line [gets $f]
				if {$line eq ""} break
				foreach {pos status ref temp alt} $line break
				unset -nocomplain count
				set total 0
				foreach {key num} [lrange $line 5 end] {
					set key [string range $key 0 end-1]
					set count($key) $num
					incr total $num
				}
				set list {}
				foreach base [array names count] {
					if {$base eq $ref} continue
					lappend list [list $base [expr {double($count($base))/$total}]]
				}
				set list [lsort -real -decreasing -index 1 $list]
				set alt [lindex $list 0 0]
				set freq [format %.3f [lindex $list 0 1]]
				set begin [expr {$pos-1}]
				set end $pos
				if {$ref eq "-"} {
					set type ins
					set ref ""
					set end $begin
				} elseif {$alt eq "-"} {
					set type del
					set end [expr {$begin+[string length $ref]}]
					set alt ""
				} else {
					set type snp
				}
				set otherallele {}
				set otherfreq {}
				list_foreach {a fr} [lrange $list 1 end] {
					if {[info exists basica($a)]} continue
					lappend otherallele $a
					lappend otherfreq [format %.4f $fr]
				}
				if {$freq < $minfreq} {
					set alleleSeq1 $ref
					set alleleSeq2 $ref
				} else {
					if {$freq > 0.8} {
						set alleleSeq1 $alt
					} else {
						set alleleSeq1 $ref
					}
					set alleleSeq2 $alt
				}
				if {$total < $coveragecutoff} {
					set sequenced u
				} elseif {$freq < $minfreq} {
					set sequenced r
				} else {
					set sequenced v
				}
				puts $o [join [list $chrom $begin $end $type $ref $alt $status $freq $alleleSeq1 $alleleSeq2 $sequenced $total [join $otherallele ,] [join $otherfreq ,] [get count(A) 0] [get count(C) 0] [get count(G) 0] [get count(T) 0] [get count(N) 0] [get count(-) 0]] \t]
			}
		}
	}	
	if {$o ne "stdout"} {catch {close $o}}
	if {$f ne "stdin"} {catch {gzclose $f}}
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
