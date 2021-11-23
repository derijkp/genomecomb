proc tcl::mathfunc::isnum value {
	isdouble $value
}

proc tcl::mathfunc::isint value {
	::isint $value
}

proc tcl::mathfunc::format args {
	::format {*}$args
}

proc tcl::mathfunc::trimformat args {
	string trimright [string trimright [::format {*}$args] 0] .
}

proc tcl::mathfunc::percent args {
	foreach {value prec} $args break
	if {![::isint $prec]} {set prec 1}
	::format %.${prec}f [expr {100*$value}]
}
