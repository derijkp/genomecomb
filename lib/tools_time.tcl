proc timescan {time} {
	if {$time eq ""} {return $time}
	time_scan $time
}

proc time_comp {time1 time2} {
	if {$time1 eq ""} {
		if {$time2 eq ""} {return 0} else {return 1}
	}
	if {$time2 eq ""} {return -1}
	set diff [expr {[lindex $time1 0] - [lindex $time2 0]}]
	if {$diff != 0} {return $diff}
	expr {[lindex $time1 1] - [lindex $time2 1]}
}
