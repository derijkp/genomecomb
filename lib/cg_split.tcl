proc cg_split {args} {
	set postfix {}
	set sorted 0
	cg_options split args {
		-f - -field {
			set field $value
		}
		-s - -sorted {
			set sorted $value
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
	set skip 0
	while {[gets $f line] >= 0} {
		set value [lindex [split $line \t] $pos]
		if {$chrclip} {set value [chr_clip $value]}
		if {![info exists a($value)]} {
			if {$sorted} {
				if {[info exists prevvalue] && !$skip} {
					gzclose $a($prevvalue)
					file rename $prefix$prevvalue$postfix.temp[gzext $postfix] $prefix$prevvalue$postfix
				}
				if {[file exists $prefix$value$postfix]} {
					puts stderr "File $prefix$value$postfix exists, skipping"
					set skip 1
				} else {
					set skip 0
				}
				set prevvalue $value
			}
			if {!$skip} {
				set o [wgzopen $prefix$value$postfix.temp[gzext $postfix]]
				puts -nonewline $o $comment
				puts $o [join $header \t]
				set a($value) $o
			}
		}
		if {$skip} continue
		puts $a($value) $line
	}
	if {$sorted} {
		if {[info exists prevvalue] && !$skip} {
			gzclose $a($prevvalue)
			file rename $prefix$prevvalue$postfix.temp[gzext $postfix] $prefix$prevvalue$postfix
		}
	} else {
		foreach value [array names a] {
			gzclose $a($value)
			file rename $prefix$value$postfix.temp[gzext $postfix] $prefix$value$postfix
		}
	}
	gzclose $f
}
