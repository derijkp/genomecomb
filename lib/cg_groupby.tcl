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
	if {([llength $args] < 1)} {
		puts "Wrong number of arguments"
		cg_groupby_help
		exit 1
	}
	set fields [lindex $args 0]
	set args [lrange $args 1 end]
	if {[llength $args] > 0} {
		set filename [lindex $args 0]
		set f [rzopen $filename]
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
	puts $o [join [list_concat $fields [list_sub $header -exclude $poss]] \t]
	set line [split [gets $f] \t]
	set prev [list_sub $line $poss]
	set grouped $line
	set next 100000; set num 0
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		incr num
		if {$num >= $next} {putslog $num; incr next 100000}
		set cur [list_sub $line $poss]
		if {$cur != $prev} {
			puts $o [join [list_concat $prev [list_sub $grouped -exclude $poss]] \t]
			set grouped $line
			set prev $cur
		} else {
			set temp {}
			foreach g $grouped v $line {
				lappend temp $g,$v
			}
			set grouped $temp
		}
	}
	puts $o [join [list_concat $prev [list_sub $grouped -exclude $poss]] \t]
	if {$o ne "stdout"} {catch {close $o}}
	if {$f ne "stdin"} {catch {close $f}}
}

if {[info exists argv]} {
	package require pkgtools
	set appdir [file dir [pkgtools::startdir]]
	lappend auto_path $appdir/lib
	append env(PATH) :[file dir [file dir $appdir]]/bin:$appdir/bin
	package require Extral
	set ::base $scriptname
	cg_groupby {*}$argv
}
