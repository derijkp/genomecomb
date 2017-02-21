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

