proc cg_split {args} {
	set postfix {}
	cg_options split args {
		-f {
			set field $value
		}
	} {file prefix postfix} 2 3
	if {$file ne "-"} {
		set f [gzopen $file]
	} else {
		set f stdin
	}
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
	gzclose $f
	foreach value [array names a] {
		close $a($value)
	}
}
