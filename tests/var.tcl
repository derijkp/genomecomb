#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test var_gatk {var_gatk basic} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg var_gatk -stack 1 tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas 2> /dev/null
	cg tsvdiff tmp/var-gatk-bwa.tsv.lz4 data/var-gatk-bwa.tsv.lz4
} {}

testsummarize

