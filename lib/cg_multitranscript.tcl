#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg_multitranscript {args} {
	cg_options multitranscript args {
	} compar_file 2
	set isoformfiles $args

	foreach file $isoformfiles {
		catch {close $a(f,$file)}
	}
	unset -nocomplain a
	set header1 {}
	set common {}
	foreach file $isoformfiles {
		set a(f,$file) [gzopen $file]
		set header [tsv_open $a(f,$file) comment]
		set poss [tsv_basicfields $header 3]
		foreach pos $poss field {chromosome begin end} {
			if {$pos == -1} continue
			lset header $pos $field
		}
		set a(h,$file) $header
		set a(id,$file) [list_cor $header {chromosome begin end strand exonStarts exonEnds}]
		if {[inlist $a(id,$file) -1]} {
			set missing [list_sub {chromosome begin end strand exonStarts exonEnds} [list_find $a(id,$file) -1]]
			error "file $file is missing essential fields: $missing"
		}
		set a(data,$file) [list_find -glob $header *-*]
		set a(empty,$file) [list_fill [llength $a(data,$file)] 0.0]
		set a(status,$file) [gets $a(f,$file) a(curline,$file)]
		set a(curline,$file) [split $a(curline,$file) \t]
		set a(curid,$file) [list_sub $a(curline,$file) $a(id,$file)]
		set other [list_sub $header -exclude [list {*}$a(id,$file) {*}$a(data,$file)]]
		if {$header1 eq ""} {
			set header1 [list_sub $header -exclude $a(data,$file)]
			set common $other
		} else {
			set common [list_common $common $other]
#			if {[list_sub $header -exclude $a(data,$file)] ne $header1} {
#				error "$file has a different basic header (not data/count fields) from [lindex $isoformfiles 0]"
#			}
		}
	}
	foreach file $isoformfiles {
		set a(common,$file) [list_cor $a(h,$file) $common]
	}
	set header {chromosome begin end strand exonStarts exonEnds}
	lappend header {*}$common
	foreach file $isoformfiles {
		lappend header {*}[list_sub $a(h,$file) $a(data,$file)]
	}
	set o [open $compar_file.temp w]
	if {$comment ne ""} {puts -nonewline $o $comment}
	puts $o [join $header \t]
	while 1 {
		set curids {}
		foreach file $isoformfiles {
			lappend curids $a(curid,$file)
		}
		set curid [lindex [bsort $curids] 0]
		if {$curid eq {{} {} {} {} {} {}}} break
		set pos [lsearch $curids $curid]
#join [bsort $curids] \n\n
#flush $o
#join $curids \n\n
#puts [list_subindex $curids 6]
#if {[lindex $curid 2] eq "32128246"} {error stop}
		set mfile [lindex $isoformfiles $pos]
		set line $curid
		lappend line {*}[list_sub $a(curline,$mfile) $a(common,$mfile)]
		foreach file $isoformfiles {
			if {$a(curid,$file) eq $curid} {
				lappend line {*}[list_sub $a(curline,$file) $a(data,$file)]
				set a(status,$file) [gets $a(f,$file) a(curline,$file)]
				set a(curline,$file) [split $a(curline,$file) \t]
				set temp [list_sub $a(curline,$file) $a(id,$file)]
				if {[lindex [bsort [list $curid $temp]] 1] ne $temp} {
					error "file $file not sorted correctly; should be sorted on: chromosome begin end strand exonStarts exonEnds"
				}
				set a(curid,$file) $temp
			} else {
				lappend line {*}$a(empty,$file)
			}
		}
		puts $o [join $line \t]
	}

	close $o
	foreach file $isoformfiles {
		catch {close $a(f,$file)}
	}
	file rename -force $compar_file.temp $compar_file

}
