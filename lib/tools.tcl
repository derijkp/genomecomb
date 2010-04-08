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

proc rzroot filename {
	if {[inlist {.rz} [file extension $filename]]} {
		return [file root $filename]
	} else {
		return $filename
	}
}

proc overlap {start1 end1 start2 end2} {
	if {$start2 >= $end1} {return [expr {$end1-$start2}]}
	if {$end2 < $start1} {return [expr {$end2-$start1}]}
	if {$start2 > $start1} {set start1 $start2}
	if {$end2 < $end1} {set end1 $end2}
	expr {$end1-$start1}
}

proc putslog {args} {
}

proc putslog {args} {
	foreach message $args {
		puts stderr $message
	}
}

proc chrindexseek {file f chr} {
	set indexfile [rzroot $file].chrindex
	if {![file exists $indexfile]} {
		set tf [rzopen $file]
		set header [gets $tf]
		set chrpos [lsearch $header chromosome]
		set prevchr {}
		set list {}
		set o [open $indexfile w]
		while {![eof $tf]} {
			set pos [tell $tf]
			set line [gets $tf]
			set chr [lindex $line $chrpos]
			if {$chr ne $prevchr} {
				puts $o $chr\t$pos
				set prevchr $chr
			}
		}
		close $tf
		close $o
	}
	set trfchrpos [split [string trim [file_read $indexfile]] \n\t]
	if {[dict exists $trfchrpos chr$chr]} {
		set fpos [dict get $trfchrpos chr$chr]
		seek $f $fpos start
	} else {
		set fpos [dict get $trfchrpos $chr]
		seek $f $fpos start
	}
}
