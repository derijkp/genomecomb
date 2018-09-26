#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test var {var_gatk basic} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg var_gatk -stack 1 tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas
	cg tsvdiff tmp/var-gatk-bwa.tsv.lz4 data/var-gatk-bwa.tsv.lz4
	string_change [cg covered tmp/sreg-gatk-bwa.tsv.lz4] [list \n\n \n]
} {chromosome	bases
chr21	1047
chr22	156
total	1203}

test var {var_distrreg gatk} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg var_distrreg -stack 1 -method gatk tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas
	cg tsvdiff tmp/var-gatk-bwa.tsv.lz4 data/var-gatk-bwa.tsv.lz4
	string_change [cg covered tmp/sreg-gatk-bwa.tsv.lz4] [list \n\n \n]
} {chromosome	bases
chr21	1047
chr22	156
total	1203}

test var {var_distrreg gatk -d 3} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg var_distrreg -stack 1 -d 3 -method gatk tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas
	cg tsvdiff tmp/var-gatk-bwa.tsv.lz4 data/var-gatk-bwa.tsv.lz4
	string_change [cg covered tmp/sreg-gatk-bwa.tsv.lz4] [list \n\n \n]
} {chromosome	bases
chr21	1047
chr22	156
total	1203}

test var {var_distrreg gatk sequencedgenome (part)} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg select -q {$chromosome regexp "chr2."} $::refseqdir/hg19/extra/reg_hg19_sequencedgenome.tsv.lz4 tmp/distrreg.tsv
	cg var_distrreg -stack 1 -method gatk -distrreg tmp/distrreg.tsv tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas
	cg tsvdiff tmp/var-gatk-bwa.tsv.lz4 data/var-gatk-bwa.tsv.lz4
	string_change [cg covered tmp/sreg-gatk-bwa.tsv.lz4] [list \n\n \n]
} {chromosome	bases
chr21	1047
chr22	156
total	1203}

test var {var_distrreg gatk size} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg var_distrreg -stack 1 -method gatk -distrreg 30000000 tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas
	cg tsvdiff tmp/var-gatk-bwa.tsv.lz4 data/var-gatk-bwa.tsv.lz4
	string_change [cg covered tmp/sreg-gatk-bwa.tsv.lz4] [list \n\n \n]
} {chromosome	bases
chr21	1047
chr22	156
total	1203}

test var {var_distrreg gatk result exists already} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	foreach file {
		tmp/sreg-cov3-bwa.tsv.lz4
		tmp/varall-gatk-bwa.tsv.lz4 tmp/varall-gatk-bwa.tsv.lz4.lz4i tmp/varall-gatk-bwa.tsv.analysisinfo
		tmp/var-gatk-bwa.tsv.lz4 tmp/var-gatk-bwa.tsv.lz4.lz4i tmp/var-gatk-bwa.tsv.analysisinfo
		tmp/sreg-gatk-bwa.tsv.lz4 tmp/sreg-gatk-bwa.tsv.lz4.lz4i tmp/sreg-gatk-bwa.tsv.analysisinfo
		tmp/reg_cluster-gatk-bwa.tsv.lz4
	} {file_write $file {}}
	cg var_distrreg -stack 1 -v 2 -distrreg chr -method gatk tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas > tmp/log 2> tmp/logerror
	catch {exec grep "targets already completed or running" tmp/logerror} temp
	if {$temp eq ""} {error "not skipped"}
	llength [glob -nocomplain tmp/*.old]
} 0

test var {var_sam basic} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg var_sam -stack 1 tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas
	cg tsvdiff tmp/var-sam-bwa.tsv.lz4 data/var-sam-bwa.tsv.lz4
	string_change [cg covered tmp/sreg-sam-bwa.tsv.lz4] [list \n\n \n]
} {chromosome	bases
chr21	1047
chr22	156
total	1203}

test var {var_distrreg sam} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg var_distrreg -stack 1 -method sam tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas
	cg tsvdiff tmp/var-sam-bwa.tsv.lz4 data/var-sam-bwa.tsv.lz4
	string_change [cg covered tmp/sreg-sam-bwa.tsv.lz4] [list \n\n \n]
} {chromosome	bases
chr21	1047
chr22	156
total	1203}

test var {var_distrreg sam -d 3} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg var_distrreg -stack 1 -d 3 -method sam tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas
	cg tsvdiff tmp/var-sam-bwa.tsv.lz4 data/var-sam-bwa.tsv.lz4
	string_change [cg covered tmp/sreg-sam-bwa.tsv.lz4] [list \n\n \n]
} {chromosome	bases
chr21	1047
chr22	156
total	1203}

test var {var_distrreg sam result exists already} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg var_distrreg -stack 1 -method sam tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas
	cg var_distrreg -stack 1 -method sam tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas > tmp/log 2> tmp/logerror
	catch {exec grep "targets already completed or running" tmp/logerror} temp
	if {$temp eq ""} {error "not skipped"}
	llength [glob -nocomplain tmp/*.old]
} 0

test var {var_freebayes basic} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg var_freebayes -stack 1 tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas 2> /dev/null
	cg tsvdiff tmp/var-freebayes-bwa.tsv.lz4 data/var-freebayes-bwa.tsv.lz4
} {}

test var {var_distrreg freebayes} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg var_distrreg -stack 1 -v 2 -method freebayes tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas
	cg tsvdiff tmp/var-freebayes-bwa.tsv.lz4 data/var-freebayes-bwa.tsv.lz4
} {}

test var {var_gatkh basic} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	# use low quality settings because test bwa.bam has low coverage
	cg var_gatkh -stack 1 -mincoverage 5 -mingenoqual 12 tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas
	cg tsvdiff tmp/var-gatkh-bwa.tsv.lz4 data/var-gatkh-bwa.tsv.lz4
	cg tsvdiff tmp/varall-gatkh-bwa.gvcf.gz data/varall-gatkh-bwa.gvcf.gz
	string_change [cg covered tmp/sreg-gatkh-bwa.tsv.lz4] [list \n\n \n]
} {chromosome	bases
chr21	1018
chr22	142
total	1160}

test var {var_distrreg gatkh} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg var_distrreg -stack 1 -cleanup 0 -distrreg 1 -method gatkh -mincoverage 5 -mingenoqual 12 tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas
	cg tsvdiff tmp/var-gatkh-bwa.tsv.lz4 data/var-gatkh-bwa.tsv.lz4
	cg tsvdiff tmp/varall-gatkh-bwa.gvcf.gz data/varall-gatkh-bwa.gvcf.gz
	string_change [cg covered tmp/sreg-gatkh-bwa.tsv.lz4] [list \n\n \n]
} {chromosome	bases
chr21	1018
chr22	142
total	1160}

testsummarize

