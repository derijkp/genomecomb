proc biobambam {args} {
	set ::env(LD_LIBRARY_PATH) $::externdir/biobambam2:$::env(LD_LIBRARY_PATH)
	exec {*}$args
}
