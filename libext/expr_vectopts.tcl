proc matchlistsize {v1 v2} {
	set v1 [split $v1 ";, "]
	set len1 [llength $v1]
	set v2 [split $v2 ";, "]
	set len2 [llength $v2]
	if {$len1 != $len2} {
		if {$len1 == 1} {
			set v1 [list_fill $len2 $v1]
		} elseif {$len2 == 1} {
			set v2 [list_fill $len1 $v2]
		} else {
			error "lists of different length: $v1 and $v2"
		}
	}
	list $v1 $v2	
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
