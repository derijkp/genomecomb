proc cg_gb2fasta {args} {
	set fields {}
	set infile -
	set outfile -
	set clean_names 1
	cg_options cg_fasta args {
		-clean_names {
			set clean_names $value
		}
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
	while 1 {
		set c {}
		while 1 {
			set read [gets $f line]
			if {$read == -1} break
			if {[regexp ^ORIGIN $line]} break
			append c $line\n
		}
		if {$read == -1} break
		if {![regexp {LOCUS *([^\n]+)} $c temp name] || [lindex $name 0] in {{} . Exported}} {
			if {![regexp {VERSION *([^\n]+)} $c temp name] || [lindex $name 0] in {{} . Exported}} {
				if {![regexp {ACCESSION *([^\n]+)} $c temp name] || [lindex $name 0] in {{} . Exported}} {
					if {![regexp {KEYWORDS *([^\n]+)} $c temp name] || [lindex $name 0] in {{} . Exported}} {
						if {![regexp {DEFINITION *([^\n]+)} $c temp name] || [lindex $name 0] in {{} . Exported}} {
							error "Could not find name (tried LOCUS, VERSION, ACCESSION, KEYWORDS, DEFINITION)"
						}
					}
				}
			}
		}
		if {$clean_names} {
			regsub -all {[^A-Za-z0-9 _]} $name _ name
		}
		puts $o >$name
		while 1 {
			set read [gets $f line]
			if {$read == -1} break
			if {[string range $line 0 1] eq "//"} break
			regsub -all {[ 0-9]} $line {} seq
			puts -nonewline $o $seq
		}
		puts $o ""
		if {$read == -1} break
	}
	if {$f ne "stdin"} {
		gzclose $f
	}
	if {$o ne "stdout"} {
		gzclose $o
		file rename $tempfile $outfile
	} else {
		# do this, otherwise does not out output last \n when run e.g. from exec
		puts ""
	}
}

proc cg_gb2fas {args} {
	cg_gb2fasta {*}$args
}

proc cg_gbk2fasta {args} {
	cg_gb2fasta {*}$args
}

proc cg_gbk2fas {args} {
	cg_gb2fasta {*}$args
}
