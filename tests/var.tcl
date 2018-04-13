#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test var {var_gatk basic} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg var_gatk -stack 1 tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas
	cg tsvdiff tmp/var-gatk-bwa.tsv.lz4 data/var-gatk-bwa.tsv.lz4
} {}

test var {var_freebayes basic} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg var_freebayes -stack 1 tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas 2> /dev/null
	cg tsvdiff tmp/var-freebayes-bwa.tsv.lz4 data/var-freebayes-bwa.tsv.lz4
} {}

test var {var_gatkh basic} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg var_gatkh -stack 1 tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas
	cg tsvdiff tmp/var-gatkh-bwa.tsv.lz4 data/var-gatkh-bwa.tsv.lz4
	exec gunzip -d tmp/varall-gatkh-bwa.gvcf.gz
	exec gunzip -d -c data/varall-gatkh-bwa.gvcf.gz > tmp/expected.gvcf
	catch {exec diff tmp/varall-gatkh-bwa.gvcf tmp/expected.gvcf | grep -v GATKCommandLine | grep -v {^12c12$}} e
	if {$e ne "---\nchild process exited abnormally"} {error $e}
} {}

test var {var_gatkh -distrchr 1} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg var_gatkh -stack 1 -distrchr 1 tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas
	cg tsvdiff tmp/var-gatkh-bwa.tsv.lz4 data/var-gatkh-bwa.tsv.lz4
	exec gunzip -d tmp/varall-gatkh-bwa.gvcf.gz
	exec gunzip -d -c data/varall-gatkh-bwa.gvcf.gz > tmp/expected.gvcf
	catch {exec diff tmp/varall-gatkh-bwa.gvcf tmp/expected.gvcf | grep -v GATKCommandLine | grep -v {^12c12$}} e
	if {$e ne "---\nchild process exited abnormally"} {error $e}
} {}

testsummarize

