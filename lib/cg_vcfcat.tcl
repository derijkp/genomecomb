proc cg_vcfcat {args} {
	set index 0
	set o stdout
	set threads 1
	cg_options catvcf args {
		-o {
			set outfile $value
		}
		-i {
			set index $value
		}
		-threads {
			set threads $value
		}
	} {} 1 ... {
		concatenate vcf files that must have the same basic header:
		The header of the first file is used without checking for compatibility!
	}
	if {$index && ![info exists outfile]} {
		error "cg vcfcat cannot index (-i 1) if no outputfile is given (-o)"
	}
	if {[info exist outfile]} {
		if {[file extension $outfile] eq ".gz"} {
			set o [open "| bgzip -c -@ $threads > $outfile.temp" w]
		} else {
			set compress [compresspipe $outfile]
			if {$compress ne ""} {
				set o [open "$compress > $outfile.temp" w]
			} else {
				set o [open $outfile.temp w]
			}
		}
	}
	if {[llength $args] == 1} {
		set f [gzopen [lindex $args 0]]
		fcopy $f $o
		close $o
		gzclose $f
	} else {
		set header 1
		foreach file $args {
			if {![file size $file]} continue
			set f [gzopen $file]
			if {$header} {
				fcopy $f $o
				set header 0
			} else {
				while {[gets $f line] != -1} {
					if {[string index $line 0] eq {#}} continue
					puts $o $line
					break
				}
				fcopy $f $o
			}
			gzclose $f
		}
		close $o
	}
	if {[info exists outfile]} {
		file rename -force -- $outfile.temp $outfile
		if {$index} {
			cg_gatk_index $outfile
		}
	}
}
