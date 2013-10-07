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
	if {$min == Inf} {return NaN}
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
	if {$min == -Inf} {return NaN}
	return $max
}

proc tcl::mathfunc::lmaxpos args {
	set max -Inf
	set maxpos NaN
	foreach v $args {
		set curpos 0
		foreach v [split $v ";, "] {
			if {[isdouble $v] && $v > $max} {
				set max $v
				set maxpos $curpos
			}
			incr curpos
		}
	}
	return $maxpos
}

proc tcl::mathfunc::lminpos args {
	set min Inf
	set minpos NaN
	foreach v $args {
		set curpos 0
		foreach v [split $v ";, "] {
			if {[isdouble $v] && $v < $min} {
				set min $v
				set minpos $curpos
			}
			incr curpos
		}
	}
	return $minpos
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
