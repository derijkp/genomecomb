proc cg_fasta2tsv {args} {
	set fields {}
	set infile -
	set outfile -
	cg_options fasta2tsv args {
		-f - -fields {
			set fields $value
		}
	} {infile outfile} 0
	if {[llength $args]} {
		set infiles [list $infile $outfile]
		set outfile [list_pop args]
		lappend infiles {*}$args
	} else {
		set infiles [list $infile]
	}
	if {$outfile ne "-"} {
		if {[file exists $outfile]} {error "file exists: $outfile"}
		set tempfile $outfile.temp[file extension $outfile]
		set o [wgzopen $tempfile]
	} else {
		set o stdout
	}
	puts $o [join {id sequence} \t]
	foreach infile $infiles {
		if {$infile ne "-"} {
			set f [gzopen $infile]
		} else {
			set f stdin
		}
		set filename [file tail $infile]
		set read [gets $f line]
		puts -nonewline $o [string range $line 1 end]\t
		while {1} {
			set read [gets $f line]
			if {$read == -1} break
			if {[string index $line 0] eq ">"} {
				puts -nonewline $o \n[string range $line 1 end]\t
			} else {
				puts -nonewline $o $line
			}
		}
		puts -nonewline $o \n
		if {$f ne "stdin"} {gzclose $f}
	}
	if {$o ne "stdout"} {
		gzclose $o
		file rename $tempfile $outfile
	} else {
		# do this, otherwise does not out output last \n when run e.g. from exec
		puts ""
	}
}