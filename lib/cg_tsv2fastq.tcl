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
	incr idpos
	set seqpos [lsearch $header seq]
	if {$seqpos == -1} {
		set seqpos [lsearch $header sequence]
	}
	incr seqpos
	set qualpos [lsearch $header qual]
	if {$qualpos == -1} {
		set qualpos [lsearch $header quality]
	}
	incr qualpos
	set awkcode [subst {print \$$idpos "\\n" \$$seqpos "\\n+\\n" \$$qualpos}]
	chanexec $f $o [list awk -F {\t} \{$awkcode\}]
}