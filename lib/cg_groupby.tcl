#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg_groupby {args} {
	set sumfields {}
	set statsfields {}
	set sorted 1
	set usefields {}
	cg_options groupby args {
		-sumfields {
			set sumfields $value
		}
		-stats {
			set statsfields $value
		}
		-sorted {
			set sorted $value
		}
		-f {
			set usefields $value
		}
	} {fields filename outfile} 1 3
	if {[info exists filename]} {
		set f [gzopen $filename]
	} else {
		set f stdin
	}
	if {[info exists outfile]} {
		set o [open $outfile w]
	} else {
		set o stdout
	}
	set header [tsv_open $f]
	set poss [list_cor $header $fields]
	set sumposs [list_cor $header $sumfields]
	set statsposs [list_cor $header $statsfields]
	set listfields [list_sub $header -exclude [list_concat $poss $sumposs $statsposs]]
	if {[llength $usefields]} {
		set listfields [list_common $listfields $usefields]
	}
	set listposs [list_cor $header $listfields]
	set header [list_concat $fields $sumfields $listfields]
	foreach field $statsfields {
		lappend header ${field}_min ${field}_total ${field}_count ${field}_max
	}
	puts $o [join $header \t]
	if {$sorted} {
		set sumposs [list_change $sumposs {-1 x}]
		chanexec $f $o [list groupby $poss $listposs $sumposs $statsposs]
	} else {
		groupby_unsorted $f $o $poss $listposs $sumposs $statsposs
	}
}

proc groupby_unsorted {f o poss listposs sumposs statsposs} {
	set next 1000000; set num 0
	unset -nocomplain a
	set start [list_concat [list_fill [llength $sumposs] 0] [list_fill [llength $listposs] {}] [list_concat [list_fill [llength $statsposs] {x 0 0 x}]]]
	set prevsumposs [list_fill [llength $sumposs] 0 1]
	set prevlistposs [list_fill [llength $listposs] [llength $sumposs] 1]
	set prevstatsposs [list_fill [expr {4*[llength $statsposs]}] [expr {[llength $sumposs]+[llength $listposs]}] 1]
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		incr num
		if {$num >= $next} {putsprogress $num; incr next 1000000}
		set cur [list_sub $line $poss]
		if {[info exists a($cur)]} {
			set prev $a($cur)
		} else {
			set prev $start
		}
		set result {}
		foreach p [list_sub $prev $prevsumposs] n [list_sub $line $sumposs] {
			if {[isint $n]} {
				lappend result [expr {$p+$n}]
			} else {
				lappend result [expr {$p+1}]
			}
		}
		foreach p [list_sub $prev $prevlistposs] n [list_sub $line $listposs] {
			if {$p eq ""} {set p $n} else {append p ,$n}
			lappend result $p
		}
		foreach {min tot count max} [list_sub $prev $prevstatsposs] n [list_sub $line $statsposs] {
			if {[isint $n]} {
				if {![isdouble $min]} {set min $n}
				if {$n < $min} {lappend result $n} else {lappend result $min}
				lappend result [expr {$tot+$n}]
				lappend result [expr {$count+1}]
				if {![isdouble $max]} {set max $n}
				if {$n > $max} {lappend result $n} else {lappend result $max}
			} else {
				lappend result $min $tot $count $max
			}
		}
		set a($cur) $result
	}
	foreach cur [ssort -natural [array names a]] {
		puts $o [join $cur \t]\t[join $a($cur) \t]
	}
	close $f
	close $o
}

if {[info exists argv0] && [file tail [info script]] eq [file tail $argv0]} {
	package require pkgtools
	set appdir [file dir [pkgtools::startdir]]
	lappend auto_path $appdir/lib
	append env(PATH) :[file dir [file dir $appdir]]/bin:$appdir/bin
	package require Extral
	set ::base $scriptname
	cg_groupby {*}$argv
}
