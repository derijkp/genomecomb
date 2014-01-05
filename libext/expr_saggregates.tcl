proc tcl::mathfunc::slist_cond_ {args} {
	set result {}
	foreach {if value} $args {
		if {$if} {lappend result $value}
	}
	return [join $result ,]
}

proc tcl::mathfunc::sdistinct_cond_ {args} {
	set result {}
	foreach {if value} $args {
		if {$if} {lappend result $value}
	}
	return [join [list_remdup $result] ,]
}

proc tcl::mathfunc::smin_cond_ {args} {
	set result {}
	foreach {if value} $args {
		if {$if} {lappend result $value}
	}
	::tcl::mathfunc::lmin $result
}

proc tcl::mathfunc::smax_cond_ {args} {
	set result {}
	foreach {if value} $args {
		if {$if} {lappend result $value}
	}
	::tcl::mathfunc::lmax $result
}

#; **smin(?condition?,value)**: returns the minimum of results of value for each sample for which (if given) **condition** is true
#; **smax(?condition?,value)**: returns the maximum of results of value for each sample for which (if given) **condition** is true
#; **ssum(?condition?,value)**: returns the sum of results of value for each sample for which (if given) **condition** is true
#; **savg(?condition?,value)**: returns the average of results of value for each sample for which (if given) **condition** is true
#; **sstddev(?condition?,value)**: returns the standard deviation of results of value for each sample for which (if given) **condition** is true
#; **smedian(?condition?,value)**: returns the median of results of value for each sample for which (if given) **condition** is true
#; **smode(?condition?,value)**: returns the mode of results of value for each sample for which (if given) **condition** is true
#; **percent(condition1,condition2)**: returns 100.0*(number of samples for which condition1 and condition2 are true)/(number of samples for which condition1 is true)
