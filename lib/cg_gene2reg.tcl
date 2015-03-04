proc cg_gene2reg {args} {
	set genecol name2
	set transcriptcol name
	if {[llength $args] > 2} {
		puts stderr "format is cg gene2reg ?genefile? ?outfile?"
		errorformat gene2reg
		exit 1
	}
	set file {}
	set outfile {}
	foreach {file outfile} $args break

	catch {close $f}; catch {close $o};
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
	set header [open_genefile $f dposs]
	set rest [list_sub $header -exclude $dposs]
	set restposs [list_cor $header $rest]
	

	puts $o [join [list_concat {chromosome begin end type element rna_start rna_end protein_start protein_end gene transcript} $rest] \t]
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set geneobj [annotatevar_gene_makegeneobj {} $line $dposs]
		set ftlist [dict get $geneobj ftlist]
		set chr [lindex $line [lindex $dposs 0]]
		set gene [dict get $geneobj genename]
		set transcript [dict get $geneobj transcriptname]
		set rest [list_sub $line $restposs]
		set ftlist [lrange [dict get $geneobj ftlist] 1 end-1]
		foreach dline $ftlist {
			lset dline 1 [expr {[lindex $dline 1]+1}]
			puts $o [join [list $chr {*}$dline $gene $transcript {*}$rest] \t]
		}
	}
	if {$o ne "stdout"} {
		close $o
		file rename -force $outfile.temp $outfile
	}
	if {$f ne "stdout"} {close $f}
}

if 0 {
	set file tmp/genetmp.tsv
	set outfile tmp/out.tsv
	set args [list $file $outfile]
}

