proc matchlist {v1 v2} {
	set v1 [split $v1 ";, "]
	set len1 [llength $v1]
	set v2 [split $v2 ";, "]
	set len2 [llength $v2]
	if {$len1 != $len2} {
		if {$len1 == 1} {
			set v1 [list_fill $len2 $v1]
		} elseif {$len2 == 1} {
			set v2 [list_fill $len2 $v2]
		} else {
			error "lists of different length: $v1 and $v2"
		}
	}
	list $v1 $v2	
}

proc tcl::mathfunc::lone {args} {
	foreach value $args {
		foreach el [split $value ";, "] {
			if {[isint $el]} {
				if {$el} {return 1}
			}
		}
	}
	return 0
}

proc tcl::mathfunc::lall {args} {
	foreach value $args {
		foreach el [split $value ";, "] {
			if {[isint $el]} {
				if {!$el} {return 0}
			}
		}
	}
	return 1
}

proc tcl::mathfunc::lcount {args} {
	set result 0
	foreach value $args {
		foreach el [split $value ";, "] {
			if {[isint $el]} {
				if {$el} {incr result}
			}
		}
	}
	return $result
}

proc tcl::mathfunc::llen {args} {
	set result 0
	foreach value $args {
		incr result [llength [split $value ";, "]]
	}
	return $result
}

proc tcl::mathfunc::vector {args} {
	join $args ,
}

proc tcl::mathfunc::lavg args {
	set data {}
	foreach v $args {
		if {[isdouble $v]} {
			lappend data $v
		} else {
			foreach v [split $v ";, "] {
				if {[isdouble $v]} {lappend data $v}
			}
		}
	}
	lmath_average $data
}

proc tcl::mathfunc::lsum args {
	set data {}
	foreach v $args {
		if {[isdouble $v]} {
			lappend data $v
		} else {
			foreach v [split $v ";, "] {
				if {[isdouble $v]} {lappend data $v}
			}
		}
	}
	lmath_sum $data
}

proc tcl::mathfunc::vif {args} {
	if {[llength $args] < 3 || [expr {[llength $args]%2}] != 1} {
		error "wrong # args for function vif, must be: vif(condition1,true1,?condition2?,?true2?,...,false)"
	}
	set args [matchlistsize {*}$args]
	set len [llength [::lindex $args 0]]
	set result {}
	for {set pos 0} {$pos < $len} {incr pos} {
		set line [list_subindex $args $pos]
		unset -nocomplain value
		foreach {cond true} $line {
			if {[true $cond]} {
				set value $true
				break
			}
		}
		if {[info exists value]} {
			lappend result $value
		} else {
			lappend result [::lindex $line end]
		}
	}
	return [join $result ,]
}

proc tcl::mathfunc::vavg {args} {
	set result [split [::lindex $args 0] ";, "]
	set result [list_fill [expr {[llength $result]*2}] 0]
	foreach value $args {
		set newresult {}
		foreach el [split $value ";, "] {p t} $result {
			if {[isdouble $el]} {
				set p [expr {$p + $el}]
				incr t
			}
			lappend newresult $p $t
		}
		set result $newresult
	}
	set newresult {}
	foreach {p t} $result {
		if {![catch {
			expr {$p/double($t)}
		} res]} {
			lappend newresult $res
		} else {
			lappend newresult NaN
		}
	}
	return [join $newresult ,]
}

proc tcl::mathfunc::vabs {value} {
	set result {}
	foreach el [split $value ";, "] {
		if {[catch {
			lappend result [expr {abs($el)}]
		}]} {
			lappend result NaN
		}
	}
	return [join $result ,]
}

proc tcl::mathfunc::vmax args {
	set len [llength [split [::lindex $args 0] ";, "]]
	set max [list_fill $len -Inf]
	foreach v $args {
		set pos 0
		foreach v [split $v ";, "] m $max {
			if {[isdouble $v] && $v > $m} {lset max $pos $v}
			incr pos
		}
	}
	return [join $max ,]
}

proc tcl::mathfunc::vmin args {
	set len [llength [split [::lindex $args 0] ";, "]]
	set min [list_fill $len Inf]
	foreach v $args {
		set pos 0
		foreach v [split $v ";, "] m $min {
			if {[isdouble $v] && $v < $m} {lset min $pos $v}
			incr pos
		}
	}
	return [join $min ,]
}
