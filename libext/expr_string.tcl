proc tcl::mathfunc::concat {args} {
	join $args ""
}

proc tcl::mathfunc::length {string} {
	string length $string
}

proc tcl::mathfunc::split {string {sep ,}} {
	::split $string $sep
}

proc tcl::mathfunc::toupper {string} {
	::string toupper $string
}
