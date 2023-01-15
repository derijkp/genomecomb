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

test counters {process_project rnaseqc and qorts} {

	test_cleantmp
	mkdir tmp/samples/s1
	cg_bam_sort data/star-p.sam tmp/samples/s1/map-sstar-s1.bam
	analysisinfo_write {} tmp/samples/s1/map-sstar-s1.bam sample s1
	mkdir tmp/samples/s2
	file copy tmp/samples/s1/map-sstar-s1.bam tmp/samples/s2/map-sstar-s2.bam
	analysisinfo_write {} tmp/samples/s2/map-sstar-s2.bam sample s2
	cg process_project -stack 1 -v 2 \
		-refdir $::refseqdir/hg19 \
		-realign 0 \
		-removeduplicates 0 \
		-varcallers {} \
		-aligners star \
		-counters {rnaseqc qorts} tmp \
		>& tmp/startuplog.txt
	cg select -overwrite 1 -q {$counts-qorts-sstar-s1 > 0 || $counts-qorts-sstar-s2 > 0} tmp/compar/gene_counts-qorts-sstar-tmp.tsv tmp/temp.tsv
	file_write tmp/expected.tsv [deindent {
		#filetype tsv/countfile
		#fileversion    0.99
		#fields field	number	type	description
		#fields geneid	1	String	id field
		#fields counts	1	Integer	gene counts per sample
		#fields count_cds	1	Integer	gene counts per sample
		#fields count_utr	1	Integer	gene counts per sample
		#fields count_ambig_gene	1	Integer	gene counts per sample
		geneid	counts-qorts-sstar-s2	count_cds-qorts-sstar-s2	count_utr-qorts-sstar-s2	count_ambig_gene-qorts-sstar-s2	counts-qorts-sstar-s1	count_cds-qorts-sstar-s1	count_utr-qorts-sstar-s1	count_ambig_gene-qorts-sstar-s1
		ENSG00000100412.11	2	2	0	14	2	2	0	14
		ENSG00000100413.12	1	0	1	14	1	0	1	14
		ENSG00000183486.8	78	61	17	0	78	61	17	0
		_total_no_feature	5	0	0	0	5	0	0	0
		_total_ambiguous	14	0	0	0	14	0	0	0
	}]\n
	exec diff tmp/temp.tsv tmp/expected.tsv
	cg select -overwrite 1 -q {$counts-rnaseqc-sstar-s1 > 0 || $counts-rnaseqc-sstar-s2 > 0} tmp/compar/gene_counts-rnaseqc-sstar-tmp.tsv tmp/temp.tsv
	file_write tmp/expected.tsv [deindent {
		#filetype tsv/countfile
		#fileversion    0.99
		#fields field	number	type	description
		#fields geneid	1	String	id field
		#fields genename	1	String	name field
		#fields counts	1	Integer	gene counts per sample
		geneid	genename	counts-rnaseqc-sstar-s2	counts-rnaseqc-sstar-s1
		ENSG00000183486.8	MX2	38	38
		ENSG00000100412.11	ACO2	1	1
		ENSG00000100413.12	POLR3H	8	8
	}]\n
	exec diff tmp/temp.tsv tmp/expected.tsv

} {}

testsummarize
