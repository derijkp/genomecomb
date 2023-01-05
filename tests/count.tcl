#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test counters {count_rnaseqc basic} {
	test_cleantmp
	cg_bam_sort data/star-p.sam tmp/ali_star.bam
	set bam tmp/ali_star.bam
	set refseq $::refseqdir/hg19
	cg count_rnaseqc -refseq $::refseqdir/hg19 tmp/ali_star.bam
	cg select -q {$counts-rnaseqc-ali_star > 0} tmp/gene_counts-rnaseqc-ali_star.tsv
} {#filetype tsv/countfile
#fileversion    0.99
#fields field	number	type	description
#fields geneid	1	String	id field
#fields genename	1	String	name field
#fields counts	1	Integer	gene counts per sample
geneid	genename	counts-rnaseqc-ali_star
ENSG00000183486.8	MX2	38
ENSG00000100412.11	ACO2	1
ENSG00000100413.12	POLR3H	8}

test counters {count_qorts basic} {
	test_cleantmp
	cg_bam_sort data/star-p.sam tmp/ali_star.bam
	set bam tmp/ali_star.bam
	set refseq $::refseqdir/hg19
	cg count_qorts -refseq $::refseqdir/hg19 tmp/ali_star.bam
	cg select -q {$counts-qorts-ali_star > 0} tmp/gene_counts-qorts-ali_star.tsv
} {geneid	counts-qorts-ali_star	count_cds-qorts-ali_star	count_utr-qorts-ali_star	count_ambig_gene-qorts-ali_star
ENSG00000100412.11	2	2	0	14
ENSG00000100413.12	1	0	1	14
ENSG00000183486.8	78	61	17	0
_total_no_feature	5	0	0	0
_total_ambiguous	14	0	0	0}

testsummarize
