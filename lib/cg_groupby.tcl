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

proc cg_groupby_help {} {
set help [file_read $::appdir/lib/cg_groupby.help]
puts [string_change $help [list @BASE@ [get ::base {[info source]}]]]
}

proc cg_groupby {args} {
	set pos 0
	set sumfields {}
	foreach {key value} $args {
		switch -- $key {
			-sumfields {
				set sumfields $value
			}
			-h - --help {
				cg_groupby_help
				exit 0
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
		puts "Wrong number of arguments"
		cg_groupby_help
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
	chanexec $f $o [list groupby $poss $listposs $sumposs]
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
