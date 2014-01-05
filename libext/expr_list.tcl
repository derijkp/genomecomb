proc tcl::mathfunc::lindex {list pos} {
	::lindex [split $list ";, "] $pos
}

proc tcl::mathfunc::slist_if_ {args} {
	set result {}
	foreach {if value} $args {
		if {$if} {lappend result $value}
	}
	return [join $result ,]
}
