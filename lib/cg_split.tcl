proc cg_split {args} {
	set pos 0
	while 1 {
		set key [lindex $args $pos]
		switch -- $key {
			-f {
				incr pos
				set field [lindex $args $pos]
			}
			-- break
			default {
				if {[string index $key 0] eq "-"} {error "unknown option \"$key\""}
				break
			}
		}
		incr pos
	}
	set args [lrange $args $pos end]
	if {[llength $args] < 2 || [llength $args] > 3} {
		errorformat split
		exit 1
	}
	set postfix {}
	foreach {file prefix postfix} $args break
	set f [gzopen $file]
	set header [tsv_open $f comment]
	if {[info exists field]} {
		set pos [lsearch $header $field]
		set chrclip 0
		if {$pos == -1} {error "field $field was not found in the file"}
	} else {
		set pos [tsv_basicfields $header 1]
		set chrclip 1
		if {$pos == -1} {error "no chromosome field was found in the file"}
	}
	while {[gets $f line] >= 0} {
		set value [lindex [split $line \t] $pos]
		if {$chrclip} {set value [chr_clip $value]}
		if {![info exists a($value)]} {
			set o [open $prefix$value$postfix w]
			puts -nonewline $o $comment
			puts $o [join $header \t]
			set a($value) $o
		}
		puts $a($value) $line
	}
	close $f
	foreach value [array names a] {
		close $a($value)
	}
}
