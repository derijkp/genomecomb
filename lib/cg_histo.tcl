proc cg_histo {args} {
	set header 1
	cg_options histo args {
		-header {set header $value}
	} pos 1 1
	if {$header} {
		set header [tsv_open stdin]
		set pos [lsearch $header $pos]
		if {$pos == -1} {error "Field $field not found"}
	}
	puts value\tcount
	set o [open "| tsv_histo $pos 0 >@ stdout" w]
	fcopy stdin $o
	close $o
}

