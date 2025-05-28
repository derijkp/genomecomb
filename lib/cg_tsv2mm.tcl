proc cg_tsv2mm {args} {
	upvar job_logdir job_logdir

	set refseq {}
	set src -
	set dest -
	set col1file {}
	set col2file {}
	set col1fields {}
	set col2fields {}
	set valuefields {}
	set limitbarcodes {}
	set valuetypes {}
	cg_options matrixmarket2tsv args {
		-col1file {set col1file $value}
		-col2file {set col2file $value}
		-col1fields {set col1fields $value}
		-col2fields {set col2fields $value}
		-valuefields {set valuefields $value}
		-limitbarcodes {set limitbarcodes $value}
		-valuetypes {set valuetypes $value}
	} {src dest} 0 2 {
		convert tsv to matrixmarket
	}

	catch {gzclose $f} ; catch {gzclose $o} ; catch {gzclose $o1} ; catch {gzclose $o2} ; 

	unset -nocomplain limitbca
	if {$limitbarcodes ne ""} {
		set f [gzopen $limitbarcodes]
		while {[gets $f barcode] != -1} {
			set limitbca($barcode) 1
		}
		gzclose $f
	}

	if {$src eq "-"} {
		set f stdin
	} else {
		set f [gzopen $src]
	}
	set header [tsv_open $f]
	if {$col1fields eq ""} {
		set col1fields [list_sub $header [tsv_basicfields $header]]
	}
	if {$col2fields eq ""} {
		set col2fields [lindex $header [lsearch $header cell]]
		if {$col2fields eq ""} {
			set col2fields [lindex $header [lsearch $header cellbarcode]]
		}
	}
	if {$valuefields eq ""} {
		set valuefields [list_sub $header [list_find -regexp $header count]]
		set valuetypes [list_fill [llength $valuefields] real]
	}
	set poss [list_cor $header $valuefields]
	set pos [lindex $poss 0]
	if {-1 in $poss} {error "some valuefields ([list_sub $valuefields [list_find $poss -1]]) not found in src $src"}
	set poss1 [list_cor $header $col1fields]
	if {-1 in $poss1} {error "some fields in col1fields not found"}
	set poss2 [list_cor $header $col2fields]
	if {-1 in $poss2} {error "some fields in col2fields not found"}
	unset -nocomplain oa
	if {$dest eq "-"} {
		if {[llength $valuefields] > 1} {
			error "Cannot have more than one valuefields for output to stdout"
		}
		set valuefield [lindex $valuefields 0]
		set oa($valuefield) stdout
	} else {
		if {[llength $valuefields] == 1} {
			set valuefield [lindex $valuefields 0]
			set oa($valuefield) [wgzopen $dest]
		} else {
			foreach valuefield $valuefields valuetype $valuetypes {
				set filenamea($valuefield) [tempfile]
				set oa($valuefield) [wgzopen $filenamea($valuefield)]
			}
		}
	}
	if {$col1file ne ""} {
		set o1 [wgzopen $col1file]
		puts $o1 [join $col1fields \t]
		set write1 1
	} else {
		set write1 0
	}
	if {$col2file ne ""} {
		set o2 [wgzopen $col2file]
		puts $o2 [join $col2fields \t]
		set write2 1
	} else {
		set write2 0
	}

	unset -nocomplain col1a
	unset -nocomplain col2a
	set col1num 0
	set col2num 0
	set total 0

	while 1 {
		if {[gets $f line] == -1} break
		set line [split $line \t]
		set v2 [list_sub $line $poss2]
		if {$limitbarcodes ne ""} {
			if {![info exists limitbca($v2)]} continue
		}
		set v1 [list_sub $line $poss1]
		set values [list_sub $line $poss]
		if {![info exists col1a($v1)]} {
			incr col1num
			set col1a($v1) $col1num
			if {$write1} {
				puts $o1 [join $v1 \t]
			}
		}
		if {![info exists col2a($v2)]} {
			incr col2num
			set col2a($v2) $col2num
			if {$write2} {
				puts $o2 [join $v2 \t]
			}
		}
		incr total
		foreach valuefield $valuefields value $values {
			puts $oa($valuefield) $col1a($v1)\t$col2a($v2)\t$value
		}
	}
	if {$write1} {gzclose $o1}
	if {$write2} {gzclose $o2}
	foreach valuefield $valuefields {
		gzclose $oa($valuefield)
	}
	gzclose $f
	foreach valuefield $valuefields {
		set filename [file root [gzroot $dest]]-$valuefield[file extension [gzroot $dest]][gzext $dest]
		set tempfile $filenamea($valuefield)
		set o [wgzopen $filename]
		puts $o "%%MatrixMarket matrix coordinate $valuetype general"
		puts $o "$col1num $col2num $total"
		set f [gzopen $tempfile]
		fcopy $f $o
		gzclose $f
		gzclose $o
		
	}
}

proc cg_tsv2matrixmarket {args} {
	cg_tsv2mm {*}$args
}
