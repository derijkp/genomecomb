proc cg_vcfcat {args} {
	set index 0
	set o stdout
	set threads 1
	set sort 0
	set sample {}
	cg_options vcfcat args {
		-o {
			set outfile $value
		}
		-i {
			set index $value
		}
		-s - -sort {
			set sort $value
		}
		-sample {
			set sample $value
		}
		-threads {
			set threads $value
		}
	} {} 1 ... {
		concatenate vcf files that must have the same basic header:
		An error will be given if the headers are different, unless
		one/some headers are strictly larger (only one containing extra lines, not
		extra/different lines in both both headers).
		In this case, the largest will be used as header for the output.
	}
	if {$index && ![info exists outfile]} {
		error "cg vcfcat cannot index (-i 1) if no outputfile is given (-o)"
	}
	set pipe {}
	# check if headers match
	unset -nocomplain firstheader
	foreach file $args {
		if {![file size $file]} continue
		set f [gzopen $file]
		set header {}
		while {[gets $f line] != -1} {
			if {[string index $line 0] ne "\#"} break
			lappend header $line
		}
		gzclose $f
		if {![info exists firstheader]} {
			set firstheader $header
			set prevfile $file
		} elseif {$header ne $firstheader} {
			if {[llength $firstheader] > [llength $header]} {
				set temp [list_lremove [lrange $header 0 end-1] $firstheader]
				set temp [list_sub $temp -exclude [list_find -regexp $temp {^##(fileDate|CommandLine|GATKCommandLine|arguments)=}]]
				if {[llength $temp]} {
					error "error concatenating vcf files: $file has a different header from $prevfile"
				}
			} else {
				set temp [list_lremove [lrange $firstheader 0 end-1] $header]
				set temp [list_sub $temp -exclude [list_find -regexp $temp {^##(fileDate|CommandLine|GATKCommandLine|arguments)=}]]
				if {[llength $temp]} {
					error "error concatenating vcf files: $file has a different header from $prevfile"
				}
				set firstheader $header
				set prevfile $file
			}
		}
	}
	set header $firstheader
	if {$sort} {
		# get headersize
		set headersize [llength $firstheader]
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
		set cols [list_pop firstheader]
		if {$firstheader ne ""} {
			puts $o [join $firstheader \n]
		}
		if {$sample eq ""} {
			puts $o $cols
		} else {
			set cols [split $cols \t]
			if {[llength $cols] > 10} {
				error "error making $target: cannot use -sample option with multisample file ($file)"
			}
			lset cols 9 $sample
			puts $o [join $cols \t]
		}
		foreach file $args {
			if {![file size $file]} continue
			set f [gzopen $file]
			while 1 {
				set r [gets $f line]
				if {$r == -1 || [string index $line 0] ne {#}} break
			}
			if {$r != -1} {
				puts $o $line
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
