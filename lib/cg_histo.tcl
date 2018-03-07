proc cg_histo {args} {
	set header 1
	cg_options histo args {
		-header {set header $value}
	} field 1 1 {
		makes a fast histogram of values in the given field. With -header 0, the position (number) must be given
	}
	if {$header} {
		set header [tsv_open stdin]
		set pos [lsearch $header $field]
		if {$pos == -1} {error "Field $field not found"}
	} else {
		set pos $field
	}
	puts value\tcount
	set o [open "| tsv_histo $pos 0 >@ stdout" w]
	fcopy stdin $o
	close $o
}

