proc tcl::mathfunc::regexp {args} {
	if {[::llength $args] < 2} {error "function regexp must have at least 2 arguments"}
	foreach {value pattern} [::lrange $args end-1 end] break
	set args [::lrange $args 0 end-2]
	::regexp {*}$args -- $pattern $value
}

proc tcl::mathfunc::ncregexp {args} {
	if {[::llength $args] < 2} {error "function regexp must have at least 2 arguments"}
	foreach {value pattern} [::lrange $args end-1 end] break
	set args [::lrange $args 0 end-2]
	::regexp -nocase {*}$args -- $pattern $value
}

proc tcl::mathfunc::regextract {args} {
	if {[::llength $args] < 2} {error "function regextract must have at least 2 arguments"}
	foreach {value} $args break
	foreach pattern [::lrange $args 1 end] {
		set matches [::regexp -inline -- $pattern $value]
		if {[::llength $matches]} break
	}
	join [lrange $matches 1 end] ,
}

proc tcl::mathfunc::regsub {args} {
	if {[::llength $args] < 3} {error "function regsub must have at least 3 arguments"}
	foreach {value pattern replace} [::lrange $args end-2 end] break
	set args [::lrange $args 0 end-3]
	::regsub {*}$args -- $pattern $value $replace temp
	return $temp
}

proc tcl::mathfunc::matches {args} {
	if {[::llength $args] < 2} {error "function matches must have at least 2 arguments"}
	foreach {value pattern} [::lrange $args end-1 end] break
	set args [::lrange $args 0 end-2]
	::string match {*}$args $pattern $value
}

# same as matches
interp alias {} tcl::mathfunc::match {} tcl::mathfunc::matches

proc tcl::mathfunc::ncmatches {args} {
	if {[::llength $args] < 2} {error "function matches must have at least 2 arguments"}
	foreach {value pattern} [::lrange $args end-1 end] break
	set args [::lrange $args 0 end-2]
	::string match -nocase {*}$args $pattern $value
}
