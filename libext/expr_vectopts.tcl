proc matchlistsize {args} {
	set data {}
	set lengths {}
	foreach v $args {
		set line [split $v ";, "]
		lappend lengths [llength $line]
		lappend data $line
	}
	set size [list_remove [list_remdup $lengths] 1]
	if {[llength $size] > 1} {
		error "some lists of different length: $args"
	}
	set result {}
	if {[llength $size] == 0} {
		return $data
	} else {
		foreach line $data len $lengths {
			if {$len == 1} {
				lappend result [list_fill $size $line]
			} else {
				lappend result $line
			}
		}
		return $result
	}
}

proc tcl::mathfunc::vminus {v1 v2} {
	foreach {v1 v2} [::matchlistsize $v1 $v2] break
	set result {}
	foreach e1 $v1 e2 $v2 {
		if {[catch {
			lappend result [expr {$e1 - $e2}]
		}]} {
			lappend result NaN
		}
	}
	return [join $result ,]
}

proc tcl::mathfunc::vplus {v1 v2} {
	foreach {v1 v2} [::matchlistsize $v1 $v2] break
	set result {}
	foreach e1 $v1 e2 $v2 {
		if {[catch {
			lappend result [expr {$e1 + $e2}]
		}]} {
			lappend result NaN
		}
	}
	return [join $result ,]
}

proc tcl::mathfunc::vtimes {v1 v2} {
	foreach {v1 v2} [::matchlistsize $v1 $v2] break
	set result {}
	foreach e1 $v1 e2 $v2 {
		if {[catch {
			lappend result [expr {$e1 * $e2}]
		}]} {
			lappend result NaN
		}
	}
	return [join $result ,]
}

proc tcl::mathfunc::vdivide {v1 v2} {
	foreach {v1 v2} [::matchlistsize $v1 $v2] break
	set result {}
	foreach e1 $v1 e2 $v2 {
		if {[catch {
			lappend result [expr {$e1 / $e2}]
		}]} {
			lappend result NaN
		}
	}
	return [join $result ,]
}

proc tcl::mathfunc::vmod {v1 v2} {
	foreach {v1 v2} [::matchlistsize $v1 $v2] break
	set result {}
	foreach e1 $v1 e2 $v2 {
		if {[catch {
			lappend result [expr {$e1 % $e2}]
		}]} {
			lappend result NaN
		}
	}
	return [join $result ,]
}

proc tcl::mathfunc::vpower {v1 v2} {
	foreach {v1 v2} [::matchlistsize $v1 $v2] break
	set result {}
	foreach e1 $v1 e2 $v2 {
		if {[catch {
			lappend result [expr {$e1 ** $e2}]
		}]} {
			lappend result NaN
		}
	}
	return [join $result ,]
}

proc tcl::mathfunc::vgt {v1 v2} {
	foreach {v1 v2} [::matchlistsize $v1 $v2] break
	set result {}
	foreach e1 $v1 e2 $v2 {
		if {[catch {
			lappend result [expr {$e1 > $e2}]
		}]} {
			lappend result NaN
		}
	}
	return [join $result ,]
}

proc tcl::mathfunc::vlt {v1 v2} {
	foreach {v1 v2} [::matchlistsize $v1 $v2] break
	set result {}
	foreach e1 $v1 e2 $v2 {
		if {[catch {
			lappend result [expr {$e1 < $e2}]
		}]} {
			lappend result NaN
		}
	}
	return [join $result ,]
}

proc tcl::mathfunc::vgte {v1 v2} {
	foreach {v1 v2} [::matchlistsize $v1 $v2] break
	set result {}
	foreach e1 $v1 e2 $v2 {
		if {[catch {
			lappend result [expr {$e1 >= $e2}]
		}]} {
			lappend result NaN
		}
	}
	return [join $result ,]
}

proc tcl::mathfunc::vlte {v1 v2} {
	foreach {v1 v2} [::matchlistsize $v1 $v2] break
	set result {}
	foreach e1 $v1 e2 $v2 {
		if {[catch {
			lappend result [expr {$e1 <= $e2}]
		}]} {
			lappend result NaN
		}
	}
	return [join $result ,]
}

proc tcl::mathfunc::veq {v1 v2} {
	foreach {v1 v2} [::matchlistsize $v1 $v2] break
	set result {}
	foreach e1 $v1 e2 $v2 {
		if {[catch {
			lappend result [expr {$e1 == $e2}]
		}]} {
			lappend result NaN
		}
	}
	return [join $result ,]
}

proc tcl::mathfunc::vne {v1 v2} {
	foreach {v1 v2} [::matchlistsize $v1 $v2] break
	set result {}
	foreach e1 $v1 e2 $v2 {
		if {[catch {
			lappend result [expr {$e1 != $e2}]
		}]} {
			lappend result NaN
		}
	}
	return [join $result ,]
}

proc tcl::mathfunc::vand {v1 v2} {
	foreach {v1 v2} [::matchlistsize $v1 $v2] break
	set result {}
	foreach e1 $v1 e2 $v2 {
		lappend result [expr {$e1 && $e2}]
	}
	return [join $result ,]
}

proc tcl::mathfunc::vor {v1 v2} {
	foreach {v1 v2} [::matchlistsize $v1 $v2] break
	set result {}
	foreach e1 $v1 e2 $v2 {
		lappend result [expr {$e1 || $e2}]
	}
	return [join $result ,]
}

proc tcl::mathfunc::vin {v1 v2} {
	foreach {v1 v2} [::matchlistsize $v1 $v2] break
	set result {}
	foreach e1 $v1 {
		lappend result [inlist $v2 $e1]
	}
	return [join $result ,]
}

proc tcl::mathfunc::vni {v1 v2} {
	foreach {v1 v2} [::matchlistsize $v1 $v2] break
	set result {}
	foreach e1 $v1 {
		if {![inlist $v2 $e1]} {
			lappend result 1
		} else {
			lappend result 0
		}
	}
	return [join $result ,]
}

