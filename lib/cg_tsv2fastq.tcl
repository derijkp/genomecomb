proc cg_tsv2fastq {args} {
	set fields {}
	set infile -
	set outfile -
	cg_options fastq2tsv args {
	} {infile outfile} 0
	if {$infile ne "-"} {
		set f [gzopen $infile]
	} else {
		set f stdin
	}
	set header [tsv_open $f]
	if {$outfile ne "-"} {
		if {[file exists $outfile]} {error "file exists: $outfile"}
		set o [wgzopen $outfile]
	} else {
		set o stdout
	}
	set idpos [lsearch $header id]
	if {$idpos == -1} {
		set idpos [lsearch $header qname]
	}
	if {$idpos == -1} {
		set idpos [lsearch $header name]
	}
	set seqpos [lsearch $header seq]
	if {$seqpos == -1} {
		set seqpos [lsearch $header sequence]
	}
	set qualpos [lsearch $header qual]
	if {$qualpos == -1} {
		set qualpos [lsearch $header quality]
	}
	set poss [list $idpos $seqpos $qualpos]
	while {[gets $f line] != -1} {
		foreach {id seq qual} [list_sub [split $line \t] $poss] break
		if {$qual eq ""} {
			set qual [string_fill ! [string length $seq]]
		}
		puts $o $id\n$seq\n+\n$qual
	}
	gzclose $f
	gzclose $o
}
