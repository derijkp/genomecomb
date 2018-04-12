proc cg_keyvalue {args} {
	set samplefields {}
	set samplefieldslen 0
	set idfields {}
	set infile {}
	set outfile {}
	cg_options keyvalue args {
		-idfields {
			set idfields $value
		}
		-samplefields {
			set samplefields $value
			set samplefieldslen [llength $samplefields]
		}
	} {file outfile} 0 2
	if {$file eq ""} {
		set f stdin
	} else {
		set f [gzopen $file]
	}
	if {$outfile eq ""} {
		set o stdout
	} else {
		set tempoutfile [filetemp $outfile]
		set o [open $tempoutfile w]
	}
	set header [tsv_open $f comment]

	if {![llength $idfields]} {
		if {[lsearch $header id] != -1} {
			set idfields id
		} else {
			set idfields sample
		}
	}
	set sample [lsearch $idfields sample]
	if {$sample != -1 && [lsearch $header sample] != -1} {
		set sample -1
	}
	set idfields [list_common $idfields $header]
	if {![llength $idfields] && $sample != -1} {
		set idfields [list_sub $header -exclude [list_find -regexp $header -]]
	}
	set keepposs [list_cor $header $idfields]
	set valuefields [list_sub $header -exclude $keepposs]
	#
	# make output
	set nh {}
	if {$sample != -1} {
		if {$samplefieldslen} {
			lappend nh {*}$samplefields
		} else {
			lappend nh sample
		}
	}
	lappend nh {*}$idfields key value
	puts $o [join $nh \t]
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set keepvalues [list_sub $line $keepposs]
		if {$sample != -1} {
			foreach field $valuefields value [list_sub $line -exclude $keepposs] {
				regexp {^([^-]+)-?(.*)$} $field temp field sample
				if {$samplefieldslen} {
					set sample [split $sample -]
					if {[llength $sample] != $samplefieldslen} {error "sample [join $sample -] has a different number of fields from $samplefieldslen"}
					set sample [join $sample \t]
				}
				puts $o $sample\t[join $keepvalues \t]\t$field\t$value
			}
		} else {
			foreach field $valuefields value [list_sub $line -exclude $keepposs] {
				puts $o [join $keepvalues \t]\t$field\t$value
			}
		}
	}
	#
	# finish
	if {$o ne "stdout"} {
		close $o
		file rename -force $tempoutfile $outfile
	}
	if {$f ne "stdout"} {gzclose $f}
}

if 0 {
	set tsvfile tests/data/reg4.tsv
	set tsvfile tests/data/testvars.tsv
	set o stdout
}
