proc eputsvars {args} {
	foreach var $args {
		if {[catch {uplevel [list set $var]} value]} {
			puts stderr [list unset $var]
		} else {
			puts stderr [list set $var $value]
		}
	}
}
