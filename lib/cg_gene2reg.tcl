proc cg_gene2reg {args} {
	set file {}
	set outfile {}
	set upstream 0
	cg_options gene2reg args {
		-upstream {set upstream $value}
	}  {file outfile} 0 2
	if {[llength $args] > 2} {
		errorformat gene2reg
	}
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
		set geneobj [annotatevar_gene_makegeneobj {} $line $dposs 2000]
		set ftlist [dict get $geneobj ftlist]
		set chr [lindex $line [lindex $dposs 0]]
		set gene [dict get $geneobj genename]
		set transcript [dict get $geneobj transcriptname]
		set rest [list_sub $line $restposs]
		if {$upstream} {
			set dline [lindex $ftlist 0]
			lset dline 1 [expr {[lindex $dline 1]+1}]
			lset dline 0 [expr {[lindex $dline 1]-$upstream}]
			lset dline 3 -1
			lset dline 4 -1
			lset dline 5 -1
			puts $o [join [list $chr {*}$dline $gene $transcript {*}$rest] \t]
		}
		foreach dline [lrange $ftlist 1 end-1] {
			lset dline 1 [expr {[lindex $dline 1]+1}]
			puts $o [join [list $chr {*}$dline $gene $transcript {*}$rest] \t]
		}
		if {$upstream} {
			set dline [lindex $ftlist end]
			lset dline 1 [expr {[lindex $dline 0]+$upstream}]
			lset dline 3 -1
			lset dline 4 -1
			lset dline 5 -1
			puts $o [join [list $chr {*}$dline $gene $transcript {*}$rest] \t]
		}
	}
	if {$o ne "stdout"} {
		close $o
		file rename -force $outfile.temp $outfile
	}
	if {$f ne "stdout"} {gzclose $f}
}

if 0 {
	set file tmp/genetmp.tsv
	set outfile tmp/out.tsv
	set args [list $file $outfile]
}

