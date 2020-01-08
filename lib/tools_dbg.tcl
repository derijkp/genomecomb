proc eputsvars {args} {
	foreach var $args {
		if {[catch {uplevel [list set $var]} value]} {
			puts stderr [list unset $var]
		} else {
			puts stderr [list set $var $value]
		}
	}
}

proc fputsvars {file args} {
	set o [open $file a]
	foreach var $args {
		if {[catch {uplevel [list set $var]} value]} {
			puts stderr [list unset $var]
			puts $o [list unset $var]
		} else {
			puts stderr [list set $var $value]
			puts $o [list set $var $value]
		}
	}
	puts stderr "-- fputsvars to $file"
	close $o
}
