proc tcl::mathfunc::ltchecked {v1 v2} {
	if {![string is double $v1]} {error "$v1 is not a number"}
	if {![string is double $v2]} {error "$v2 is not a number"}
	expr {$v1 < $v2}
}

proc tcl::mathfunc::ltechecked {v1 v2} {
	if {![string is double $v1]} {error "$v1 is not a number"}
	if {![string is double $v2]} {error "$v2 is not a number"}
	expr {$v1 <= $v2}
}

proc tcl::mathfunc::gtchecked {v1 v2} {
	if {![string is double $v1]} {error "$v1 is not a number"}
	if {![string is double $v2]} {error "$v2 is not a number"}
	expr {$v1 > $v2}
}

proc tcl::mathfunc::gtechecked {v1 v2} {
	if {![string is double $v1]} {error "$v1 is not a number"}
	if {![string is double $v2]} {error "$v2 is not a number"}
	expr {$v1 >= $v2}
}
