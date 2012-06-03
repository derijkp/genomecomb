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
	set makemap 0
	set concatlen -1
	set econcatlen 0
	set aconcat ""
	set aconcatlen 0
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
				set concatlen [string length $concat]
			}
			-e - -ce - --concatend {
				set econcat $value
				set econcatlen [string length $econcat]
			}
			-ca - --concatadj {
				set aconcat $value
				set aconcatlen [string length $econcat]
			}
			-m - --mapfile {
				set mapfile $value
				set makemap 1
			}
			--namefield {
				set namefield $value
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
	if {[info exists namefield]} {
		set namepos [lsearch $header $namefield]
		if {$namepos == -1} {
			puts stderr "namefield $namefiled not found"
			exit 1
		}
	} else {
		set namepos [lsearch $header name]
	}
	if {$concatlen >= 0} {
		puts "\>$regionfile concatenated"
		set name $regionfile
		set firstline 1
	} else {
		set name concat
	}
	if {$makemap} {
		set fm [open $mapfile w]
		puts $fm [join {chromosome begin end destchromosome destbegin destend name} \t]
	}
	set fstart 0
	if {$concatlen >= 0 && $econcatlen} {
		puts -nonewline $econcat
		incr fstart $econcatlen
	}
	set pchr {}
	set pend {}
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set sub [list_sub $line $poss]
		putslog $sub
		foreach {chr estart eend} $sub break
		regsub ^chr $chr {} chr
		set seq [genome_get $fg $chr [expr {$estart}] [expr {$eend}]]
		set seq [genome_mask $dbdir $seq $chr [expr {$estart}] [expr {$eend}] $freql $freqN $delsize $repeats]
		
		if {$concatlen == -1} {
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
			puts \n\>$name
		} elseif {!$firstline} {
			foreach {chr begin end} $sub break
			if {[comparechr $chr $pchr] == 0 && $begin == $pend} {
				puts -nonewline $aconcat
				incr fstart $aconcatlen
			} else {
				puts -nonewline $concat
				incr fstart $concatlen
			}
			set pchr $chr
			set pend $end
		}
		puts -nonewline $seq
		if {$makemap} {
			puts $fm $name\t$fstart\t[expr {$fstart+[string length $seq]}]\t[join $sub \t]\t[lindex $line $namepos]
		}
		incr fstart [string length $seq]
		set firstline 0
	}
	if {$concatlen >= 0 && $econcatlen} {
		puts -nonewline $econcat
		incr fstart $econcatlen
	}
	puts ""
	if {$makemap} {
		close $fm
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
