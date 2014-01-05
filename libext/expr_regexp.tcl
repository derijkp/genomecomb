proc tcl::mathfunc::regexp {args} {
	if {[llength $args] < 2} {error "function regexp must have at least 2 arguments"}
	foreach {value pattern} [lrange $args end-1 end] break
	set args [lrange $args 0 end-2]
	::regexp {*}$args $pattern $value
}

proc tcl::mathfunc::ncregexp {args} {
	if {[llength $args] < 2} {error "function regexp must have at least 2 arguments"}
	foreach {value pattern} [lrange $args end-1 end] break
	set args [lrange $args 0 end-2]
	::regexp -nocase {*}$args $pattern $value
}

proc tcl::mathfunc::matches {args} {
	if {[llength $args] < 2} {error "function matches must have at least 2 arguments"}
	foreach {value pattern} [lrange $args end-1 end] break
	set args [lrange $args 0 end-2]
	::string match {*}$args $pattern $value
}

proc tcl::mathfunc::ncmatches {args} {
	if {[llength $args] < 2} {error "function matches must have at least 2 arguments"}
	foreach {value pattern} [lrange $args end-1 end] break
	set args [lrange $args 0 end-2]
	::string match -nocase {*}$args $pattern $value
}
