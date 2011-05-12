#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

if 0 {
cd ~/dev/completegenomics/tests
set file table2.tsv
set fields s1

cg groupby s1 < ~/dev/completegenomics/tests/table2.tsv
cg select -q '$symbol != ""' /complgen/projects/test/annottest_compar.tsv | cg groupby symbol > /complgen/projects/test/groupby_symbol_test_compar.tsv
}

proc cg_groupby {args} {
	set pos 0
	set sumfields {}
	set sorted 1
	foreach {key value} $args {
		switch -- $key {
			-sumfields {
				set sumfields $value
			}
			-sorted {
				set sorted $value
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	if {([llength $args] < 1)} {
		errorformat groupby
		exit 1
	}
	set fields [lindex $args 0]
	set args [lrange $args 1 end]
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
	set header [tsv_open $f]
	set poss [list_cor $header $fields]
	set sumposs [list_cor $header $sumfields]
	set listfields [list_sub $header -exclude [list_concat $poss $sumposs]]
	set listposs [list_cor $header $listfields]
	puts $o [join [list_concat $fields $sumfields $listfields] \t]
	if {$sorted} {
		chanexec $f $o [list groupby $poss $listposs $sumposs]
	} else {
		groupby_unsorted $f $o $poss $listposs $sumposs
	}
}

proc groupby_unsorted {f o poss listposs sumposs} {
	set next 1000000; set num 0
	unset -nocomplain a
	set start [list_concat [list_fill [llength $sumposs] 0] [list_fill [llength $listposs] {}]]
	set prevsumposs [list_fill [llength $sumposs] 0 1]
	set prevlistposs [list_fill [llength $listposs] [llength $sumposs] 1]
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
			lappend result [expr {$p+$n}]
		}
		foreach p [list_sub $prev $prevlistposs] n [list_sub $line $listposs] {
			if {$p eq ""} {set p $n} else {append p ,$n}
			lappend result $p
		}
		set a($cur) $result
	}
	foreach cur [lsort -dict [array names a]] {
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
