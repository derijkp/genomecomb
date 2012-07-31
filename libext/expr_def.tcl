proc tcl::mathfunc::def {value def} {
	if {[isdouble $value]} {return $value} else {return $def}
}

proc tcl::mathfunc::vdef {value def} {
	set result {}
	foreach el [split $value ,] {
		if {[isdouble $el]} {
			lappend result $el
		} else {
			lappend result $def
		}
	}
	return [join $result ,]
}
