proc tcl::mathfunc::between {value args} {
	if {[llength $args] == 1} {
		foreach {min max} [::lindex $args 0] break
	} else {
		foreach {min max} $args break
	}
	if {$value >= $min && $value <= $max} {return 1} else {return 0}
}
