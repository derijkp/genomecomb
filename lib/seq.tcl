proc seq_gc {seq {window {}}} {
	set len [string length $seq]
	if {$window eq "" || $len <= $window} {
		set num [regexp -all {[GCgc]} $seq]
		expr {100.0*$num/$len}
	} else {
		set result {}
		set sub [string range $seq 0 [expr {$window-1}]]
		set num [regexp -all {[GCgc]} $sub]
		lappend result [expr {100.0*$num/$window}]
		set prev 0
		for {set cur $window} {$cur < $len} {incr cur} {
			set gone [string index $seq $prev]
			if {[regexp {[GCgc]} $gone]} {
				incr num -1
			}
			set new [string index $seq $cur]
			if {[regexp {[GCgc]} $new]} {
				incr num
			}
			lappend result [expr {100.0*$num/$window}]
			incr prev
		}
		return $result	
	}
}

