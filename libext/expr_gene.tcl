proc tcl::mathfunc::codingcat {args} {
	if {[regexp CDS $args]} {
		return C
	} elseif {[regexp UTR $args]} {
		return U
	} elseif {[regexp RNA $args]} {
		return R
	} elseif {[regexp splice $args]} {
		return s
	} else {
		return -
	}
}

proc tcl::mathfunc::coding {args} {
	tcl::mathfunc::codingcat $args
}

proc tcl::mathfunc::zyg args {
	::zyg {*}$args
}

