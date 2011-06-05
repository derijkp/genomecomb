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
	cd /complgen/tests/sequenom
	set dbdir /complgen/refseq/hg18
	set compar_file test.csv
	set resultfile sequenom.tsv
	cg makesequenom test.csv sequenom.tsv /complgen/refseq/hg18

	cd /complgen/tests/sequenom
	set dbdir /complgen/refseq/hg19
	set compar_file data/testvars.tsv
	cg makesequenom data/testvars.tsv sequenom.tsv /complgen/refseq/hg19
}


proc cg_makesequenom {args} {
	set extraseq 124
	set freql 0
	set freqN 0.2
	set pos 0
	foreach {key value} $args {
		switch -- $key {
			-f - --freq {
				set freql $value
			}
			-n - --freqn {
				set freqN $value
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	if {([llength $args] != 3)} {
		errorformat makesequenom
		exit 1
	}
	foreach {compar_file resultfile dbdir} $args break
	#
	catch {close $f}; catch {close $fg}
	set fg [genome_open [lindex [glob $dbdir/genome_*.ifas] 0]]
	set f [gzopen $compar_file]
	set header [tsv_open $f]
	set poss [tsv_basicfields $header]
	set o [open $resultfile w]
	puts $o SNP_ID\tSEQUENCE
	set dbsnpfiles [gzfiles $dbdir/var_*snp*.tsv.gz]
	set dbsnpposs {}
	foreach dbsnp $dbsnpfiles {
		set dbsnpheader [cg select -h $dbsnp]
		set temp [tsv_basicfields $dbsnpheader 4]
		lappend temp [lsearch $dbsnpheader freq]
		lappend dbsnpposs $temp
	}
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set sub [list_sub $line $poss]
		putslog $sub
		foreach {chr start end type ref alt} $sub break
		regsub ^chr $chr {} chr
		set name [join [list_sub $sub {0 1 2}] -]
		set estart [expr {$start-$extraseq}]
		if {$estart < 0} {set estart 0}
		set eend [expr {$end+$extraseq-1}]
		set seq [genome_get $fg $chr [expr {$estart}] [expr {$eend}]]
		set rstart [expr {$start-$estart}]
		set rend [expr {$rstart+$end-$start-1}]
		set test [string range $seq $rstart $rend]
		if {[string toupper $ref] ne [string toupper $test]} {
			error "ref in vars ($ref) different from ref in genome ($test) for:\n$sub"
		}
		# repeats are already masked
		# mask snps
		set list {}
		foreach snpposs $dbsnpposs dbsnp $dbsnpfiles {
			set temp [split [exec tabix $dbsnp chr$chr:$estart-$eend] \n]
			lappend list {*}[list_subindex $temp $snpposs]
		}
		set list [lsort -dict -decreasing $list]
		list_foreach {c s e type freq} $list {
			set start [expr {$s-$estart}]
			if {$freq eq ""} {set freq 0}
			if {$freq <= $freql} continue
			if {$type eq "ins"} {
				set end $start
			} else {
				set end [expr {$e-$estart-1}]
				if {$end < $start} {set end $start}
			}
			if {$type eq "del" && [expr {$end-$start}] > 5} continue
			set base [string range $seq $start $end]
			if {$freq > $freqN} {
				regsub -all . $base N base
			} else {
				set base [string tolower $base]
			}
			set seq [string_replace $seq $start $end $base]
		}
		if {$ref eq ""} {set ref -}
		if {$alt eq ""} {set alt -}
		set list [list $ref {*}[split $alt ,]]
		set list [list_change $list {{} -}]
		set seq [string_replace $seq $rstart $rend \[[join $list /]\]]
		puts $o $name\t$seq
		flush $o
	}
	close $o
}

if {[info exists argv0] && [file tail [info script]] eq [file tail $argv0]} {
	package require pkgtools
	set appdir [file dir [pkgtools::startdir]]
	lappend auto_path $appdir/lib
	append env(PATH) :[file dir [file dir $appdir]]/bin:$appdir/bin
	package require Extral
	set ::base $scriptname
	cg_makesequenom {*}$argv
}
