proc timescan {timeVar {warning {}}} {
	upvar $timeVar time
	if {$time eq ""} {return $time}
	if {$warning eq ""} {
		return [time_scan $time]
	} else {
		if {[catch {
			set code [time_scan $time]
		} msg]} {
			puts stderr "$warning: $msg"
			set code {} ; set time {}
		}
		return $code
	}
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
