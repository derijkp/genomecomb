package require Extral

proc covered regfile {
	set f [opencgifile $regfile header]
	if {[lrange [split $header \t] 0 2] ne "chromosome begin end"} {
		error "header error in $regfile"
	}
	set num 0
	unset -nocomplain a
	while {![eof $f]} {
		incr num
		if {![expr $num%100000]} {puts stderr $num}
		set line [cggets $f]
		foreach {chr start end} $line break
		if {[info exists a($chr)]} {
			set a($chr) [expr {$a($chr) + $end - $start}]
		} else {
			set a($chr) [expr {$end - $start}]
		}
	}
	close $f
	puts chr\tbases
	set list [list_remove [lsort -dict [array names a]] X Y M]
	lappend list X Y M
	set total 0
	foreach chr $list {
		if {[info exists a($chr)]} {
			puts $chr\t$a($chr)
			incr total $a($chr)
		}
	}
	puts ""
	puts total\t$total
}
