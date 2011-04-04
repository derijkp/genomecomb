# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

package require Extral

proc cg_covered args {
	if {[llength $args] > 1} {
		errorformat covered
		exit 1
	}
	if {[llength $args] > 0} {
		set regfile [lindex $args 0]
		set f [gzopen $regfile]
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
