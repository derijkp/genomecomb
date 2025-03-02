proc cg_tsv2fasta {args} {
	set fields {}
	set infile -
	set outfile -
	set idfield {}
	set sequencefield {}
	cg_options cg_fasta args {
		-idfield {set idfield $value}
		-sequencefield {set sequencefield $value}
	} {infile outfile} 0 2
	if {$infile ne "-"} {
		set f [gzopen $infile]
	} else {
		set f stdin
	}
	if {$outfile ne "-"} {
		if {[file exists $outfile]} {error "file exists: $outfile"}
		set tempfile $outfile.temp[file extension $outfile]
		set o [wgzopen $tempfile w]
	} else {
		set o stdout
	}
	set header [tsv_open $f]
	if {$idfield eq ""} {
		set idfield [list_common {id qname name} $header]
	}	
	if {$sequencefield eq ""} {
		set sequencefield [list_common {sequence seq} $header]
	}	
	set fields [list $idfield $sequencefield]
	set cor [list_cor $header $fields]
	if {-1 in $cor} {
		error "tsv2fasta missing fields: [list_sub $fields [list_find $cor -1]]"
	}
	while {1} {
		set read [gets $f line]
		if {$read == -1} break
		foreach {id seq} [list_sub [split $line \t] $cor] break
		puts $o \>$id\n$seq
	}
	if {$f ne "stdin"} {
		gzclose $f
	}
	if {$o ne "stdout"} {
		gzclose $o
		file rename $tempfile $outfile
	} else {
		# do this, otherwise does not out output last \n when run e.g. from exec
		flush stdout
	}
}