proc cg_vcfcat {args} {
	set index 0
	set o stdout
	set threads 1
	set sort 0
	set sample {}
	set mergeheaders 1
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
		-mergeheaders {
			set mergeheaders [true $value]
		}
	} {} 1 ... {
		concatenate vcf files that must have the same basic header:
		If -mergeheaders is 1 (default), the headers of the different
		vcf files are merged (if they have differences)
	}
	if {$index && ![info exists outfile]} {
		error "cg vcfcat cannot index (-i 1) if no outputfile is given (-o)"
	}
	set pipe {}
	# check if headers match
	if {$mergeheaders} {
		unset -nocomplain headera
		set keys {}
		set colslist {}
		foreach file $args {
			if {![file size $file]} continue
			set f [gzopen $file]
			while {[gets $f line] != -1} {
				if {[string index $line 0] ne "\#"} break
				if {![regexp {##([^=]+)=} $line temp key]} {
					set key ""
					if {[regexp ^#CHROM\t $line]} {
						lappend colslist $line
						continue
					}
				}
				if {[info exists headera($key)] && $key in {fileDate CommandLine GATKCommandLine arguments}} continue
				list_addnew headera($key) $line
				list_addnew keys $key
			}
			gzclose $f
		}
		if {[llength [get headera(fileformat) ""]] > 1} {
			error "error concatenating vcf files: different versions of vcfs (according to headers)"
		}
		set len [llength [split [lindex $colslist 0] \t]]
		foreach line $colslist {
			set tlen [llength [split $line \t]]
			if {$tlen != $len} {
				error "error concatenating vcf files: different number of samples (according to header)"
			}
		}
		set header {}
		foreach key $keys {
			foreach line $headera($key) {
				lappend header $line
			}
		}
		lappend header [lindex $colslist 0]
	} else {
		set header {}
		foreach file [lindex $args 0] {
			if {[file size $file]} break
		}
		set f [gzopen $file]
		while {[gets $f line] != -1} {
			if {[string index $line 0] ne "\#"} break
			lappend header $line
		}
		gzclose $f
	}
	if {$sort} {
		# get headersize
		set headersize [llength $header]
		append pipe "| gnusort8 --header-lines $headersize --parallel $threads -T \"[scratchdir]\" --buffer-size=500M --compress-program=zstd-mt-1 -t \\t -s -N "
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
		set cols [list_pop header]
		if {$header ne ""} {
			puts $o [join $header \n]
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
