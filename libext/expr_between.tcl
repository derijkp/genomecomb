proc tcl::mathfunc::between {value list} {
	foreach {min max} $list break
	if {$value >= $min && $value <= $max} {return 1} else {return 0}
}
