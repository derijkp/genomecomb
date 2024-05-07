proc cg_genepredtsv2fasta {args} {
	set file {}
	set outfile {}
	set refseq {}
	cg_options genepredtsv2fasta args {
		-refseq {set refseq $value}
	}  {file outfile} 0 2
	if {[llength $args] > 2} {
		errorformat genepredtsv2fasta
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
	set refseq [refseq $refseq]
	set genomef [genome_open $refseq]
	set header [open_genefile $f dposs]
	set idposs [list_remove [list {*}[lrange $dposs end-1 end] {*}[lrange $dposs 0 end-2]] -1]
	set idfields [list_sub $header $idposs]
	set rest [list_sub $header -exclude $dposs]
	set restposs [list_cor $header $rest]
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		foreach {chrom start end strand exonStarts exonEnds cdsStart cdsEnd transcriptname genename} [list_sub $line $dposs] break
		set geneobj [annotatevar_gene_makegeneobj $genomef $line $dposs 2000 1]
		set seq [annotatevar_gene_rnaseq geneobj]
		set id [join [list_sub $line $idposs] |]
		puts $o ">$transcriptname gene:$genename chromosome:$chrom begin:$start end:$end strand:$strand exonStarts:$exonStarts exonEnds:$exonEnds cdsStart:$cdsStart cdsEnd:$cdsEnd"
		puts $o "$seq"
	}
	if {$o ne "stdout"} {
		close $o
		file rename -force -- $outfile.temp $outfile
	}
	if {$f ne "stdout"} {gzclose $f}
}

if 0 {
	set file tmp/genetmp.tsv
	set outfile tmp/out.tsv
	set args [list $file $outfile]
}


