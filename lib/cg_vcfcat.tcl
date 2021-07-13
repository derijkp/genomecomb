proc cg_vcfcat {args} {
	set index 0
	set o stdout
	set threads 1
	set sort 0
	cg_options catvcf args {
		-o {
			set outfile $value
		}
		-i {
			set index $value
		}
		-s - -sort {
			set sort $value
		}
		-threads {
			set threads $value
		}
	} {} 1 ... {
		concatenate vcf files that must have the same basic header:
		The header of the first (non-empty) file is used without checking for compatibility!
	}
	if {$index && ![info exists outfile]} {
		error "cg vcfcat cannot index (-i 1) if no outputfile is given (-o)"
	}
	set pipe {}
	if {$sort} {
		# get headersize
		foreach file $args {
			if {![file size $file]} continue
			set f [gzopen $file]
			set header {}
			while {[gets $f line] != -1} {
				if {[string index $line 0] ne "\#"} break
				lappend header $line
			}
			gzclose $f
			break
		}
		set headersize [llength $header]
		append pipe "| gnusort8 --header-lines $headersize --parallel $threads -T \"[scratchdir]\" -t \\t -s -N "
	}
	if {[info exist outfile]} {
		set compress [compresspipe $outfile {} $threads]
		if {$compress ne ""} {
			append pipe $compress
		}
		if {$pipe ne ""} {
			set o [open "$pipe > $outfile.temp" w]
		} else {
			set o [open $outfile.temp w]
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
