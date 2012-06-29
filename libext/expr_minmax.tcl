proc tcl::mathfunc::min args {
	set min Inf
	foreach v $args {
		if {[isdouble $v]} {
			if {$v < $min} {set min $v}
		} else {
			foreach v [split $v ";, "] {
				if {[isdouble $min] && $v < $min} {set min $v}
			}
		}
	}
	return $min
}

proc tcl::mathfunc::max args {
	set max -Inf
	foreach v $args {
		if {[isdouble $v]} {
			if {$v > $max} {set max $v}
		} else {
			foreach v [split $v ";, "] {
				if {[isdouble $v] && $v > $max} {set max $v}
			}
		}
	}
	return $max
}
