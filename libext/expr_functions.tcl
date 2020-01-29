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

proc tcl::mathfunc::median args {
	::median $args
}

proc tcl::mathfunc::q1 args {
	set args [bsort $args]
	set len [::llength $args]
	set len [expr {$len/2}]
	if {[expr {$len % 2}]} {
		# uneven
		set pos [expr {($len-1)/2}]
		return [lindex $args $pos]
	} else {
		set pos [expr {$len/2-1}]
		return [lmath_average [::lrange $args $pos [expr {$pos+1}]]]
	}
}

proc tcl::mathfunc::q3 args {
	set args [bsort $args]
	set len [::llength $args]
	set len [expr {$len/2}]
	if {[expr {$len % 2}]} {
		# uneven
		set pos [expr {($len-1)/2}]
		return [lindex $args end-$pos]
	} else {
		set pos [expr {$len/2}]
		return [lmath_average [::lrange $args end-$pos end-[expr {$pos-1}]]]
	}
}
