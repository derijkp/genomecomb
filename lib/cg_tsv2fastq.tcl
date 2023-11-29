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
	set commentspos [lsearch $header comments]
	set poss [list $idpos $seqpos $qualpos $commentspos]
	while {[gets $f line] != -1} {
		foreach {id seq qual comments} [list_sub [split $line \t] $poss] break
		if {$qual eq ""} {
			set qual [string_fill ! [string length $seq]]
		}
		if {$comments ne ""} {
			set id "$id\t$comments"
		}
		puts $o $id\n$seq\n+\n$qual
	}
	gzclose $f
	gzclose $o
}
