proc opencgifile {file headerVar {numlinesVar {}}} {
	if {$numlinesVar ne ""} {
		upvar $numlinesVar numlines
	}
	set numlines 0
	global cache
	upvar $headerVar header
	if {[file extension $file] eq ".gz"} {
		set f [open "|zcat $file"]
	} else {
		set f [open $file]
	}
	while {![eof $f]} {
		set line [gets $f]
		incr numlines
		if {[string length $line] && [string index $line 0] ne "#"} break
	}
	if {[string index $line 0] eq ">"} {
		set header [string range $line 1 end]
	} else {
		set header $line
	}
	while {![eof $f]} {
		set line [split [gets $f] \t]
		incr numlines
		if {[llength $line]} break
	}
	set cache($f) $line
	incr numlines -1
	return $f
}

proc opencgifile {file headerVar {numlinesVar {}}} {
	if {$numlinesVar ne ""} {
		upvar $numlinesVar numlines
	}
	set numlines 0
	global cache
	upvar $headerVar header
	if {[file extension $file] eq ".gz"} {
		set f [open "|zcat $file"]
	} else {
		set f [open $file]
	}
	while {![eof $f]} {
		set line [gets $f]
		incr numlines
		if {[string length $line] && [string index $line 0] ne "#"} break
	}
	if {[string index $line 0] eq ">"} {
		set header [string range $line 1 end]
	} else {
		set header $line
	}
	while {![eof $f]} {
		set line [split [gets $f] \t]
		incr numlines
		if {[llength $line]} break
	}
	set cache($f) $line
	incr numlines -1
	return $f
}

proc cggets {f} {
	global cache
	set line {}
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {[llength $line]} break
	}
	set result $cache($f)
	set cache($f) $line
	return $result
}

proc assert {check message} {
	if {![uplevel expr $check]} {
		error $message
	}
}

proc readcgimap {f} {
	global cache
	set line $cache($f)
	while {![eof $f]} {
		foreach {flags chromosome offsetInChr gap1 gap2 gap3 weight mateRec} $line break
		set end [expr {$offsetInChr + 35 + $gap1 + $gap2 + $gap3}]
		binary scan $weight c weight
		incr weight -33
		set lastdnbrecord [expr {$flags & 0x01}]
		if {[expr {$flags & 0x02}]} {set side r} else {set side l}
		if {[expr {$flags & 0x04}]} {set strand -} else {set strand +}
		lappend list [list $flags $chromosome $offsetInChr $end $strand $side $gap1 $gap2 $gap3 $weight $mateRec]
		set line [split [gets $f] \t]
		if {$lastdnbrecord} {
			set cache($f) $line
			return $list
		}
	}
}

proc region2bins {start end {level 0}} {
	if {$level == 0} {
		set bins 0
		set stop 0
	} else {
		set bins {}
		array set levels {1 1 2 9 3 73 4 585 5 4681}
		set stop $levels($level)
	}
	set pre 4681
	set start [expr {$start >> 14}]
	set end [expr {$end >> 14}]
	while {$pre > $stop} {
		for {set i $start} {$i <= $end} {incr i} {
			lappend bins [expr {$pre+$i}]
		}
		set pre [expr {$pre >> 3}]
		set start [expr {$start >> 3}]
		set end [expr {$end >> 3}]
	}
	return [lsort -integer $bins]
}

proc histogram {list aVar} {
	upvar $aVar a
	unset -nocomplain a
	foreach el $list {
		if {![info exists a($el)]} {
			set a($el) 1
		} else {
			incr a($el)
		}
	}
}

proc max {args} {
	lmath_max $args
}

proc min {args} {
	lmath_min $args
}

proc opensqlite3 {dbfile query} {
	set f [open "| sqlite3 -separator \"\t\" $dbfile \"$query\""]
}

proc rzopen {file {pos -1}} {
	if {[inlist {.rz} [file extension $file]]} {
		if {$pos == -1} {
			set f [open "| razip -d -c $file"]
		} else {
			set f [open "| razip -d -c -b $pos $file"]
		}
	} else {
		set f [open $file]
		if {$pos != -1} {
			seek $f $pos
		}
	}
	return $f
}
