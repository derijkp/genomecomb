proc cg_long {args} {
	if {[llength $args] > 2} {
		puts stderr "format is cg long ?tsvfile? ?outfile?"
		# errorformat select
		exit 1
	}
	foreach {file outfile} $args break
	if {$file eq ""} {
		set f stdin
	} else {
		set f [gzopen $file]
	}
	if {$outfile eq ""} {
		set o stdout
	} else {
		set o [open $outfile.temp w]
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
	foreach sample $samples {
		set cor [list_cor $fieldsa($sample) $samplecols]
		set cor [list_sub $posa($sample) $cor]
		set cor [list_change $cor {{} -1}]
		lappend todo [list_concat $aposa(pre) $cor $aposa(post)]
	}
	puts $o [join [list_concat sample $afieldsa(pre) $samplecols $afieldsa(post)] \t]
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		foreach sample $samples cor $todo {
			puts $o $sample\t[join [list_sub $line $cor] \t]
		}
	}
	if {$o ne "stdout"} {
		close $o
		file rename -force $outfile.temp $outfile
	}
	if {$f ne "stdout"} {close $f}
}

if 0 {
	set tsvfile tests/data/reg4.tsv
	set tsvfile tests/data/testvars.tsv
	set o stdout
}
