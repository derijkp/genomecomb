proc cg_zcat {args} {
	if {[lindex $args 0] in {-p -pos}} {
		set pos [lindex $args 1]
		set args [lrange $args 2 end]
		foreach file $args {
			catchchildkilled_exec {*}[gzcatra $file $pos] >@ stdout 2>@ stderr
			set pos 0
		}
	} elseif {[lindex $args 0] eq "--"} {
		set args [lrange $args 1 end]
	}
	foreach file $args {
		catchchildkilled_exec {*}[gzcat $file] $file >@ stdout 2>@ stderr
	}
}
