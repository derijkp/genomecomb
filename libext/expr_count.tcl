proc tcl::mathop::~ {value pattern} {
	regexp $pattern $value
}

proc tcl::mathfunc::hasone {list operand value} {
	set result 1
	set list [split $list ",'\;"]
	set cmd tcl::mathfunc::$operand
	foreach v $list {
		if {[$cmd $v $value]} {
			return 1
		}
	}
	return 0
}

proc tcl::mathfunc::hasall {list operand value} {
	set result 1
	set list [split $list ",'\;"]
	set cmd tcl::mathfunc::$operand
	foreach v $list {
		if {![$cmd $v $value]} {
			return 0
		}
	}
	return 1
}
