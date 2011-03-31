# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

package require Extral

proc cg_covered_help {} {
	set help [file_read $::appdir/lib/cg_covered.help]
	puts [string_change $help [list @BASE@ [get ::base {[info source]}]]]
}

proc cg_covered args {
	set pos 0
	foreach {key value} $args {
		switch -- $key {
			-h - --help {
				cg_covered_help
				exit 0
			}
			default {
				break
			}
		}
		incr pos 2
	}
	if {$pos} {set args [lrange $args $pos end]}
	if {[llength $args] > 1} {
		puts "Wrong number of arguments"
		cg_covered_help
		exit 1
	}
	if {[llength $args] > 0} {
		set regfile [lindex $args 0]
		set f [rzopen $regfile]
	} else {
		set f stdin
	}
	set poss [open_region $f]
	chanexec $f stdout "covered $poss"
}

if {[info exists argv0] && [file tail [info script]] eq [file tail $argv0]} {
	package require pkgtools
	set appdir [file dir [pkgtools::startdir]]
	lappend auto_path $appdir/lib
	append env(PATH) :[file dir [file dir $appdir]]/bin:$appdir/bin
	package require Extral
	set ::base [file tail [info script]]
	cg_select {*}$argv
}
