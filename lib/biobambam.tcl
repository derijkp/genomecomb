proc biobambam {args} {
	if (![info exists ::biobambamset]) {
		if {[info exists ::env(LD_LIBRARY_PATH)]} {
			set ::env(LD_LIBRARY_PATH) $::externdir/biobambam2:${::env(LD_LIBRARY_PATH)}
		} else {
			set ::env(LD_LIBRARY_PATH) $::externdir/biobambam2
		}
		set ::biobambamset 1
	}
	exec -ignorestderr {*}$args
}
