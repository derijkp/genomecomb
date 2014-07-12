proc cg_long {args} {
	set norm 0
	set pos 0
	foreach {key value} $args {
		switch -- $key {
			-norm {set norm $value}
			default break
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	if {[llength $args] > 2} {
		errorformat long
		exit 1
	}
	set file {}
	set outfile {}
	foreach {file outfile} $args break
	if {$file eq ""} {
		set f stdin
	} else {
		set f [gzopen $file]
	}
	if {$outfile eq ""} {
		if {$norm} {
			error "You cannot use the -norm option without giving an outfile"
		}
		set o stdout
	} else {
		set o [open $outfile.temp w]
		if {$norm} {
			set no [open $outfile.sampledata.tsv.temp w]
		}
	}
	set header [tsv_open $f comment]

	unset -nocomplain fieldsa
	unset -nocomplain posa
	unset -nocomplain afieldsa
	unset -nocomplain aposa
	set aposa(post) {}
	set aposa(pre) {}
	set afieldsa(post) {}
	set afieldsa(pre) {}
	set samples {}
	set samplecols {}
	set fpos pre
	set colpos 0
	foreach col $header {
		set pos [string first - $col]
		if {$pos != -1} {
			set field [string range $col 0 [expr {$pos-1}]]
			list_addnew samplecols $field
			incr pos
			set sample [string range $col $pos end]
			lappend samples $sample
			lappend fieldsa($sample) $field
			lappend posa($sample) $colpos
			set fpos post
		} else {
			lappend afieldsa($fpos) $col
			lappend aposa($fpos) $colpos
		}
		incr colpos
	}
	set samples [list_remdup $samples]
	set samplecols [list_union [list_common {sequenced zyg quality totalScore1 totalScore2 alleleSeq1 alleleSeq2} $samplecols] $samplecols]
	set todo {}
	set maincor [list_concat $aposa(pre) $aposa(post)]
	foreach sample $samples {
		set cor [list_cor $fieldsa($sample) $samplecols]
		set cor [list_sub $posa($sample) $cor]
		set cor [list_change $cor {{} -1}]
		if {!$norm} {
			lappend todo [list_concat $aposa(pre) $cor $aposa(post)]
		} else {
			lappend todo $cor
		}
	}
	set common [list_common $samplecols $afieldsa(pre)]
	if {[llength $common]} {
		foreach field $common {
			set newfield ${field}_global
			set num 1
			while {[inlist $samplecols $newfield]} {
				set newfield ${field}_global$num
				incr num
			}
			set afieldsa(pre) [list_change $afieldsa(pre) [list $field $newfield]]
		}
	}
	set common [list_common $samplecols $afieldsa(post)]
	if {[llength $common]} {
		foreach field $common {
			set newfield ${field}_global
			set num 1
			while {[inlist $samplecols $newfield]} {
				set newfield ${field}_global$num
				incr num
			}
			set afieldsa(post) [list_change $afieldsa(post) [list $field $newfield]]
		}
	}
	if {!$norm} {
		puts $o [join [list_concat sample $afieldsa(pre) $samplecols $afieldsa(post)] \t]
	} else {
		puts $o [join [list_concat id $afieldsa(pre) $afieldsa(post)] \t]
		puts $no [join [list_concat id sample $samplecols] \t]
	}
	if {!$norm} {
		while {![eof $f]} {
			set line [split [gets $f] \t]
			if {![llength $line]} continue
			foreach sample $samples cor $todo {
				puts $o $sample\t[join [list_sub $line $cor] \t]
			}
		}
	} else {
		set id 1
		while {![eof $f]} {
			set line [split [gets $f] \t]
			if {![llength $line]} continue
			puts $o $id\t[join [list_sub $line $maincor] \t]
			foreach sample $samples cor $todo {
				puts $no $id\t$sample\t[join [list_sub $line $cor] \t]
			}
			incr id
		}
	}
	if {$o ne "stdout"} {
		close $o
		file rename -force $outfile.temp $outfile
	}
	if {$norm} {
		file rename -force $outfile.sampledata.temp $outfile.sampledata.tsv
	}
	if {$f ne "stdout"} {close $f}
}

if 0 {
	set tsvfile tests/data/reg4.tsv
	set tsvfile tests/data/testvars.tsv
	set o stdout
}
