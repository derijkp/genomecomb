#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

if 0 {
	# testdata
	mkdir $smalltestdir/ori/rnaseqc
	cd $smalltestdir/ori/rnaseqc
	exec wget https://storage.googleapis.com/agraubert/broadinstitute/rnaseqc/test_inputs.tar.gz
	exec tar xzf test_inputs.tar.gz
	file delete test_inputs.tar.gz
}

test count_rnaseqc {count_rnaseqc basic} {
	test_cleantmp
	cg_bam_sort data/star-p.sam tmp/ali_star.bam
	set bam tmp/ali_star.bam
	set refseq $::refseqdir/hg38
	cg count_rnaseqc -refseq $::refseqdir/hg38 tmp/ali_star.bam
	cg select -q {$counts-rnaseqc-ali_star > 0} tmp/counts-rnaseqc-ali_star.tsv
} {#filetype tsv/countfile
#fileversion    0.99
#fields field	number	type	description
#fields geneid	1	String	id field
#fields genename	1	String	name field
#fields counts	1	Integer	gene counts per sample
geneid	genename	counts-rnaseqc-ali_star
ENSG00000160191.18	PDE9A	5
ENSG00000225731.1	PDE9A-AS1	2
ENSG00000159958.7	TNFRSF13C	25}

testsummarize

