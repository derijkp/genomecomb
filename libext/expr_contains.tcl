proc tcl::mathfunc::contains {list value} {
	set result 1
	set list [split $list ";, "]
	foreach v $list {
		if {$v eq $value} {
			return 1
		}
	}
	return 0
}

proc tcl::mathfunc::shares {list valuelist} {
	set result 1
	set list [split $list ";, "]
	foreach v [split $valuelist ";, "] {set a($v) 1}
	foreach v $list {
		if {[info exists a($v)]} {
			return 1
		}
	}
	return 0
}

proc tcl::mathfunc::oneof args {
	set value [list_shift args]
	::inlist $args $value
}
