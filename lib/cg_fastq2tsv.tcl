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
		puts $o [join $fields \t]
		foreach infile $infiles {
			if {$infile ne "-"} {
				set f [gzopen $infile]
			} else {
				set f stdin
			}
			set filename [file tail $infile]
			while {![eof $f]} {
				set nameline [split [gets $f] \t]
				set name [lindex $nameline 0]
				set comments [lrange $nameline 1 end]
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
					} elseif {$field eq "comments"} {
						lappend result $comments
					} else {
						set pos [list_find -glob $comments $field:*]
						if {$pos == -1} {
							lappend result ""
						} else {
							set value [lindex $comments $pos]
							regsub "$field:\[^:\]*:" $value {} value
							lappend result $value
							set comments [list_sub $comments -exclude $pos]
						}
					}
				}
				puts $o [join $result \t]
			}
			if {$f ne "stdin"} {gzclose $f}
		}
	} else {
		puts $o [join {name sequence quality comments} \t]
		foreach infile $infiles {
			if {$infile ne "-"} {
				set f [gzopen $infile]
			} else {
				set f stdin
			}
			set filename [file tail $infile]
			while {![eof $f]} {
				set nameline [split [gets $f] \t]
				set name [lindex $nameline 0]
				set comments [lrange $nameline 1 end]
				set seqline [gets $f]
				set temp [gets $f]
				set qualityline [gets $f]
				if {![string length $qualityline]} break
				puts $o $name\t$seqline\t$qualityline\t$comments
			}
			if {$f ne "stdin"} {gzclose $f}
		}
	}
	if {$outfile ne "-"} {
		gzclose $o
		file rename -force -- $outfile.temp $outfile
	}
}