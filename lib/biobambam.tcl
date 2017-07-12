proc biobambam {args} {
	if {[info exists ::env(LD_LIBRARY_PATH)]} {
		set ::env(LD_LIBRARY_PATH) $::externdir/biobambam2:${::env(LD_LIBRARY_PATH)}
	} else {
		set ::env(LD_LIBRARY_PATH) $::externdir/biobambam2
	}
	exec {*}$args
}
