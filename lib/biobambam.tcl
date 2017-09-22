proc biobambam {args} {
	if {[catch {exec {*}$args} msg] && ![regexp {\[V\] MemUsage} $msg]} {
		error $msg
	}
}
