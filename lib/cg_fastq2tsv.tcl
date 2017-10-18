proc cg_fastq2tsv {args} {
	set fields {}
	set infile -
	set outfile -
	cg_options fastq2tsv args {
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
	set header [join {id sequence temp quality} \t]
	if {$outfile ne "-"} {
		if {[file exists $outfile]} {error "file exists: $outfile"}
		set compresspipe [compresspipe $outfile]
		if {$compresspipe ne ""} {
			set o [open [list {*}$compresspipe > $outfile.temp] w]
		} else {
			set o [open $outfile.temp w]
		}
	} else {
		set o stdout
	}
	if {[llength $fields]} {
		if {$outfile ne "-"} {
			set o [open $outfile w]
		} else {
			set o stdout
		}
		puts $o [join $fields \t]
		foreach infile $infiles {
			if {$infile ne "-"} {
				set f [gzopen $infile]
			} else {
				set f stdin
			}
			set filename [file tail $infile]
			while {![eof $f]} {
				set nameline [gets $f]
				set seqline [gets $f]
				set temp [gets $f]
				set qualityline [gets $f]
				if {![string length $qualityline]} break
				set result {}
				foreach field $fields {
					if {$field eq "name"} {
						lappend result $name
					} elseif {$field eq "id"} {
						lappend result [lindex $nameline 0]
					} elseif {$field eq "sequence"} {
						lappend result $seqline
					} elseif {$field eq "avgquality"} {
						binary scan $qualityline c* qualities
						set qualities [lmath_calc $qualities - 33]
						lappend result [format %3f [lmath_average $qualities]]
					} elseif {$field eq "quality"} {
						lappend result $qualityline
					} elseif {$field eq "readlength"} {
						lappend result [string length $seqline]
					} elseif {[regexp "$field=(\[^\\n \]+)(\[\n \]|$)" $nameline temp value]} {
						lappend result $value
					} elseif {$field eq "file"} {
						lappend result $filename
					} else {
						set value ""
					}
				}
				puts $o [join $result \t]
			}
			if {$f ne "stdin"} {gzclose $f}
		}
		if {$o ne "stdout"} {close $o}
	} elseif {$infile eq "-"} {
		puts $header
		exec paste - - - - <@ stdin >@ stdout
	} else {
		puts $o $header
		foreach infile $infiles {
			exec {*}[gzcat $infile] $infile | paste - - - - >@ $o
		}
		if {$outfile ne "-"} {
			close $o
			file rename -force $outfile.temp $outfile
		}
	}
}