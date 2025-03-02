proc cg_fasta2cramref {args} {
	cg_options fasta2cramref args {
	} {fastafile outdir} 2 2 {
		fill a reference directory with md5 coded sequences from a fasta file for
		use in cram (REF_PATH)
	}
	set tail [file tail $fastafile]
	# make new one
	if {![file exists $outdir]} {file mkdir $outdir}
	set f [gzopen $fastafile]
	set line [gets $f]
	set id [lindex [string range $line 1 end] 0]
	set o [open $outdir/temp w]
	if {![file exists $outdir/mapping.tsv]} {
		set oinfo [open $outdir/mapping.tsv w]
		puts $oinfo [join {reffile chromosome md5 size} \t]
	} else {
		set oinfo [open $outdir/mapping.tsv a]
	}
	puts [string range $line 1 end]
	while 1 {
		set line [gets $f]
		if {[string index $line 0] eq ">" || [eof $f]} {
			close $o
			set md5 [lindex [exec md5sum $outdir/temp] 0]
			file rename -force -- $outdir/temp $outdir/$md5
			puts $oinfo [join [list $tail $id $md5 [file size $outdir/$md5]] \t]
			if {[eof $f]} break
			puts [string range $line 1 end]
			set id [lindex [string range $line 1 end] 0]
			set o [open $outdir/temp w]
		} elseif {$line eq ""} {
			continue
		} else {
			puts -nonewline $o [string toupper $line]
		}
	}
	gzclose $f
	close $oinfo
}
