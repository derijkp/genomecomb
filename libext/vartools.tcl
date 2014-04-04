proc concat_numericvect {arguments} {
	set data {}
	foreach v $arguments {
		if {[isdouble $v]} {
			lappend data $v
		} else {
			foreach v [split $v ";, "] {
				if {[isdouble $v]} {lappend data $v}
			}
		}
	}
	return $data
}

proc concat_vect {arguments} {
	set data {}
	foreach v $arguments {
		lappend data {*}[split $v ";, "]
	}
	return $data
}

proc median {list} {
	set list [lsort -real $list]
	set len [llength $list]
	if {[expr {$len%2}]} {
		lindex $list [expr {$len/2}]
	} else {
		set pos [expr {$len/2}]
		lmath_average [::lrange $list [expr {$pos-1}] $pos]
	}
}

proc mode {list} {
	foreach el $list {
		if {![info exists a($el)]} {
			set a($el) 1
		} else {
			incr a($el)
		}
	}
	set max 0
	set result {}
	foreach v [array names a] {
		if {$a($v) > $max} {
			set max $a($v)
			set result [list $v]
		} elseif {$a($v) == $max} {
			lappend result $v
		}
	}
	return $result
}
