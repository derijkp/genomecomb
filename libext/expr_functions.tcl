proc tcl::mathfunc::avg args {
	set tot 0
	set num 0
	foreach v $args {
		if {[isdouble $v]} {
			set tot [expr {$tot + $v}]
			incr num
		}
	}
	if {$num == 0} {
		return NaN
	} else {
		expr {double($tot)/$num}
	}
}

proc tcl::mathfunc::sum args {
	set tot 0
	set num 0
	foreach v $args {
		if {[isdouble $v]} {
			set tot [expr {$tot + $v}]
			incr num
		}
	}
	if {$num == 0} {
		return NaN
	} else {
		expr $tot
	}
}
