proc cg_histo {args} {
	set header 1
	cg_options histo args {
		-header {set header $value}
	} 1 1
	if {$header} {
		foreach {field} $args break
		set header [tsv_open stdin]
		set pos [lsearch $header $field]
		if {$pos == -1} {error "Field $field not found"}
	} else {
		set pos [lindex $args 0]
	}
	puts value\tcount
	set o [open "| tsv_histo $pos 0 >@ stdout" w]
	fcopy stdin $o
	close $o
}

