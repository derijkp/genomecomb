proc tcl::mathfunc::lmin args {
	set min Inf
	foreach v $args {
		if {[isdouble $v]} {
			if {$v < $min} {set min $v}
		} else {
			foreach v [split $v ";, "] {
				if {[isdouble $v] && $v < $min} {set min $v}
			}
		}
	}
	return $min
}

proc tcl::mathfunc::lmax args {
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

proc tcl::mathfunc::lmind {args} {
	set def [list_pop args]
	set min Inf
	foreach v $args {
		if {[isdouble $v]} {
			if {$v < $min} {set min $v}
		} else {
			foreach v [split $v ";, "] {
				if {![isdouble $v]} {set v $def}
				if {$v < $min} {set min $v}
			}
		}
	}
	return $min
}

proc tcl::mathfunc::lmaxd {args} {
	set def [list_pop args]
	set max -Inf
	foreach v $args {
		if {[isdouble $v]} {
			if {$v > $max} {set max $v}
		} else {
			foreach v [split $v ";, "] {
				if {![isdouble $v]} {set v $def}
				if {$v > $max} {set max $v}
			}
		}
	}
	return $max
}
