#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

package require Extral

if 0 {
	set dbdir /complgen/refseq/hg19/
	set regionfile data/reg_genome_seq.tsv
}


proc cg_genome_seq {args} {
	set extraseq 124
	set freql 0
	set freqN 0.2
	set delsize 5
	set repeats s
	set gc -1
	set id {}
	set pos 0
	foreach {key value} $args {
		switch -- $key {
			-f - --freq {
				set freql $value
			}
			-n - --freqn {
				set freqN $value
			}
			-d - --delsize {
				set delsize $value
			}
			-r - --repeatmasker {
				set repeats $value
			}
			-i - --id {
				set id $value
			}
			-g - --gc {
				set gc $value
			}
			-c - --concat {
				set concat $value
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	if {([llength $args] != 2)} {
		errorformat genome_seq
		exit 1
	}
	foreach {regionfile dbdir} $args break
	#
	catch {close $f}; catch {close $fg}
	set fg [genome_open [lindex [glob $dbdir/genome_*.ifas] 0]]
	set f [gzopen $regionfile]
	set header [tsv_open $f]
	set poss [tsv_basicfields $header 3]
	if {$id ne ""} {
		set idpos [lsearch $header $id]
		if {$idpos == -1} {
			error "id Column $id not found in header"
		}
	} else {
		set idpos -1
	}
	if {[info exists concat]} {
		puts "\>$regionfile concatenated"
		set firstline 1
	}
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set sub [list_sub $line $poss]
		putslog $sub
		foreach {chr estart eend} $sub break
		regsub ^chr $chr {} chr
		set seq [genome_get $fg $chr [expr {$estart}] [expr {$eend}]]
		set seq [genome_mask $dbdir $seq $chr [expr {$estart}] [expr {$eend}] $freql $freqN $delsize $repeats]
		
		if {![info exists concat]} {
			set name [join [list_sub $sub {0 1 2}] -]
			if {$idpos != -1} {
				set name "[lindex $line $idpos] $name"
				if {$gc == 0} {
					append name " GC:[format %.1f [seq_gc $seq]]"
				} elseif {$gc != -1} {
					set maxgc [lmath_max [seq_gc $seq $gc]]
					append name " GC:[format %.1f [seq_gc $seq]] maxGC($gc):[format %.1f $maxgc]"
				}
			}
			puts \>$name
		} elseif {!$firstline} {
			puts $concat
		}
		puts $seq
		set firstline 0
	}
	close $f; close $fg
}

if {[info exists argv0] && [file tail [info script]] eq [file tail $argv0]} {
	package require pkgtools
	set appdir [file dir [pkgtools::startdir]]
	lappend auto_path $appdir/lib
	append env(PATH) :[file dir [file dir $appdir]]/bin:$appdir/bin
	package require Extral
	set ::base $scriptname
	cg_genome_seq {*}$argv
}
