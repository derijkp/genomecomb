proc cg_zcat {args} {
	foreach file $args {
		exec {*}[gzcat $file] $file >@ stdout 2>@ stderr
	}
}
