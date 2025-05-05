proc cg_mm2tsv {args} {
	upvar job_logdir job_logdir
	set refseq {}
	set src -
	set dest -
	set col1file {}
	set col2file {}
	set col1fields {}
	set col2fields {}
	cg_options matrixmarket2tsv args {
		-col1file {set col1file $value}
		-col2file {set col2file $value}
		-col1fields {set col1fields $value}
		-col2fields {set col2fields $value}
	} {src dest} 0 2 {
		convert matrixmarket format to tsv
	}
	if {$src eq "-"} {
		set f stdin
	} else {
		set f [gzopen $src]
	}
	if {$dest eq "-"} {
		set o stdout
	} else {
		set o [wgzopen $dest]
	}
	if {$col1file eq "" && $col2file eq ""} {
		puts $o m\tn\tvalue
		while 1 {
			if {[gets $f line] == -1} break
			if {[string index $line 0] ne "%"} break
		}
		if {$line eq ""} return
		while 1 {
			puts $o $line
			if {[gets $f line] == -1} break
		}
	} else {
		unset -nocomplain col1a
		unset -nocomplain col2a
		if {$col1file ne ""} {
			set ff [gzopen $col1file]
			set header1 [tsv_open $ff]
			if {$col1fields eq ""} {set col1fields $header1}
			set corr [list_cor $header1 $col1fields]
			if {-1 in $corr} {error "unknown fields in -col1fields"}
			set num 1
			while {[gets $ff line] != -1} {
				set col1a($num) [list_sub $line $corr]
				incr num
			}
			gzclose $ff
		}
		if {$col2file ne ""} {
			set ff [gzopen $col2file]
			set header2 [tsv_open $ff]
			if {$col2fields eq ""} {set col2fields $header1}
			set corr [list_cor $header2 $col2fields]
			if {-1 in $corr} {error "unknown fields in -col2fields"}
			set num 1
			while {[gets $ff line] != -1} {
				set col2a($num) [list_sub $line $corr]
				incr num
			}
			gzclose $ff
		}
		puts $o [join $col1fields \t]\t[join $col2fields \t]\tvalue
		while 1 {
			if {[gets $f line] == -1} break
			if {[string index $line 0] ne "%"} break
		}
		if {$line eq ""} return
		foreach {m n value} [split $line \t] break
		puts stderr "col1nr: $m col2nr $n datanr: $value"
		while 1 {
			if {[gets $f line] == -1} break
			foreach {m n value} [split $line \t] break
			puts $o [join [get col1a($m) $m] \t]\t[join [get col2a($n) $n] \t]\t$value
		}
	}
	if {$src ne "-"} {
		gzclose $f
	}
	if {$dest ne "-"} {
		gzclose $o
	}
}

proc cg_matrixmarket2tsv {args} {
	cg_mm2tsv {*}$args
}
