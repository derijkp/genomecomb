proc tcl::mathfunc::slist_cond_ {args} {
	set result {}
	foreach {if value} $args {
		if {$if} {lappend result $value}
	}
	join $result ,
}

proc tcl::mathfunc::sdistinct_cond_ {args} {
	set result {}
	foreach {if value} $args {
		if {$if} {set a($value) 1}
	}
	join [array names a] ,
}

proc tcl::mathfunc::ucount {args} {
	llength [list_remdup $args]
}

proc tcl::mathfunc::sucount_cond_ {args} {
	set result {}
	foreach {if value} $args {
		if {$if} {set a($value) 1}
	}
	llength [array names a]
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

proc tcl::mathfunc::ssum_cond_ {args} {
	set result {}
	foreach {if value} $args {
		if {$if} {lappend result $value}
	}
	::tcl::mathfunc::lsum $result
}

proc tcl::mathfunc::savg_cond_ {args} {
	set result {}
	foreach {if value} $args {
		if {$if} {lappend result $value}
	}
	::tcl::mathfunc::lavg $result
}

proc tcl::mathfunc::sstdev_cond_ {args} {
	set result {}
	foreach {if value} $args {
		if {$if} {lappend result $value}
	}
	::tcl::mathfunc::lstddev $result
}

proc tcl::mathfunc::smedian_cond_ {args} {
	set result {}
	foreach {if value} $args {
		if {$if} {lappend result $value}
	}
	::tcl::mathfunc::lmedian $result
}

proc tcl::mathfunc::smode_cond_ {args} {
	set result {}
	foreach {if value} $args {
		if {$if} {lappend result $value}
	}
	::tcl::mathfunc::lmode $result
}

proc tcl::mathfunc::spercent_ {args} {
	set total 0
	set selected 0
	foreach {tot sel} $args {
		if {$tot} {
			incr total
			if {$sel} {incr selected}
		}
	}
	if {$total == 0} {return NaN}
	expr {100.0*$selected/$total}
}
