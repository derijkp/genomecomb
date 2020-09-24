proc tcl::mathfunc::lt {v1 v2} {
	if {[nat_compare $v1 $v2] < 0} {return 1} else {return 0}
}

proc tcl::mathfunc::le {v1 v2} {
	if {[nat_compare $v1 $v2] <= 0} {return 1} else {return 0}
}

proc tcl::mathfunc::gt {v1 v2} {
	if {[nat_compare $v1 $v2] > 0} {return 1} else {return 0}
}

proc tcl::mathfunc::ge {v1 v2} {
	if {[nat_compare $v1 $v2] >= 0} {return 1} else {return 0}
}
