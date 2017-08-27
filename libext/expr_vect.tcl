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

proc tcl::mathfunc::lindex {list pos} {
	::lindex [split $list ";, "] $pos
}

proc tcl::mathfunc::lsearch {list item args} {
	::lsearch {*}$args [split $list ";, "] $item
}

proc tcl::mathfunc::lrange {list start end} {
	::join [::lrange [split $list ";, "] $start $end] ,
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

proc tcl::mathfunc::distinct {args} {
	join [list_remdup $args] ,
}

proc tcl::mathfunc::vdistinct args {
	join [distinct {*}[concat_vect $args]] ,
}

proc tcl::mathfunc::lavg args {
	lmath_average [concat_numericvect $args]
}

proc tcl::mathfunc::lstddev args {
	lmath_stdev [concat_numericvect $args]
}

proc tcl::mathfunc::lstdev args {
	lmath_stdev [concat_numericvect $args]
}

proc tcl::mathfunc::lsum args {
	lmath_sum [concat_numericvect $args]
}

proc tcl::mathfunc::lmedian args {
	median [concat_numericvect $args]
}

proc tcl::mathfunc::lmode args {
	join [mode [concat_vect $args]] ,
}

proc tcl::mathfunc::vfunc {args} {
	set function [lindex $args 0]
	set temp {}
	set len 1
	foreach value [lrange $args 1 end] {
		set s [split $value ";, "]
		set l [llength $s]
		if {$l != 1} {
			if {$len != 1} {error "vconcat error: $value has a different number of elements"}
			set len $l
		}
		lappend temp $s
	}
	set result {}
	if {$len == 1} {
		foreach v $temp {
			lappend result [tcl::mathfunc::$function $v]
		}
	} else {
		set i 0
		foreach line $temp {
			if {[llength $line] == 1} {lset temp $i [list_fill $len [lindex $line 0]]}
			incr i
		}
		for {set i 0} {$i < $len} {incr i} {
			lappend result [tcl::mathfunc::$function {*}[list_subindex $temp $i]]
		}
	}
	return [join $result ,]
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

proc tcl::mathfunc::vformat {args} {
	set format [list_shift args]
	set result {}
	set len 1
	set todo {}
	foreach el $args {
		set el [split $el ";, "]
		set testlen [llength $el]
		if {$testlen != 1} {
			if {$len == 1} {set len $testlen} elseif {$testlen != $len} {error "args of vformat are vectors of different size (not 1)"}
		}
		lappend todo $el
	}
	if {$len != 1} {
		set pos 0
		foreach el $todo {
			set testlen [llength $el]
			if {$testlen == 1} {
				lset todo $pos [list_fill $len $el]
			}
			incr pos
		}
	}
	for {set pos 0} {$pos < $len} {incr pos} {
		set args [list_subindex $todo $pos]
		lappend result [::format $format {*}$args]
	}
	::join  $result ,
}

