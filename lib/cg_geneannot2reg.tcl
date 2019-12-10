proc cg_geneannot2reg {args} {
	set genecol name2
	set transcriptcol name
	if {[llength $args] != 2 && [llength $args] != 3} {
		errorformat geneannot2reg
	}
	set file {}
	set outfile {}
	foreach {generegfile geneannotfile outfile} $args break
	set genefiles [lrange $args 1 end]
	set tempfile [tempfile]

	catch {close $f}; catch {close $o};
	if {$outfile eq ""} {
		set o stdout
	} else {
		set o [open $outfile.temp w]
	}

	set f [gzopen $geneannotfile]
	set header [tsv_open $f]
	set genepos [lsearch $header gene]
	if {$genepos == -1} {set genepos [lsearch $header {HGNC symbol}]}
	if {$genepos == -1} {set genepos 0}
	# only take first other field (stage support for multiple fields later)
	set fields [lrange [list_sub $header -exclude $genepos] 0 0]
	set fieldsposs [list_cor $header $fields]
	unset -nocomplain a
	while {[gets $f line] != -1} {
		set line [split $line \t]
		set key [lindex $line $genepos]
		set values [list_sub $line $fieldsposs]
		foreach value $values field $fieldsposs {
			if {$value eq ""} continue
			lappend a($key,$field) $value
		}
	}
	gzclose $f
	set f [gzopen $generegfile]
	set header [tsv_open $f]
	set poss [tsv_basicfields $header 3]
	lappend poss [lsearch $header geneid]
	puts $o [join {chromosome begin end gene annotation} \t]
	while {[gets $f line] != -1} {
		set line [split $line \t]
		set line [list_sub $line $poss]
		set genes [split [lindex $line end] ,]
		set annot {}
		foreach value $values field $fieldsposs {
			set list {}
			foreach gene $genes {
				lappend list {*}[get a($gene,$field) ""]
			}
			lappend annot [join $list ,]
		}
		if {[list_remove $annot {}] eq ""} continue
		puts $o [join $line \t]\t[join $annot \t]
	}
	gzclose $f
	if {$o ne "stdout"} {
		close $o
		file rename -force -- $outfile.temp $outfile
	}
}
