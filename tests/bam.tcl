#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test bam_histo {bam_histo} {
	test_cleantmp
	file copy data/reg_hg19_smallpartexome.tsv tmp/regfile.tsv
	set namecol info
	set regionfile tmp/regfile.tsv
	set bamfile genomecomb.testdata/exomes_yri_chr2122.reference/samples/NA19240chr2122/map-rdsbwa-NA19240chr2122.bam
	set intervals {1 5 10 20 50 100 200 500 1000}
	cg bam_histo -n $namecol $regionfile $bamfile $intervals > tmp/result.tsv
	exec diff tmp/result.tsv genomecomb.testdata/expected/bam_histo-NA19240_smallpartchr2122.tsv
} {}

test bam_histo {bam_histo chr_clipped} {
	test_cleantmp
	cg select -f {chromosome=chr_clip($chromosome) begin end info} data/reg_hg19_smallpartexome.tsv tmp/regfile.tsv
	set regionfile tmp/regfile.tsv
	set bamfile genomecomb.testdata/exomes_yri_chr2122.reference/samples/NA19240chr2122/map-rdsbwa-NA19240chr2122.bam
	set intervals {1 5 10 20 50 100 200 500 1000}
	cg bam_histo -n info $regionfile $bamfile $intervals > tmp/result.tsv
	exec diff tmp/result.tsv genomecomb.testdata/expected/bam_histo-NA19240_smallpartchr2122.tsv
} {}

test bam2fastq {bam2fastq} {
	test_cleantmp
	set bamfile genomecomb.testdata/exomes_yri_chr2122.reference/samples/NA19240chr2122/map-rdsbwa-NA19240chr2122.bam
	cg bam2fastq $bamfile tmp/out_1.fq tmp/out_2.fq
	exec gzip tmp/out_1.fq tmp/out_2.fq
	exec diff tmp/out_1.fq.gz genomecomb.testdata/expected/NA19240chr2122_1.fq.gz
	exec diff tmp/out_2.fq.gz genomecomb.testdata/expected/NA19240chr2122_2.fq.gz
} {}

test bam2fastq {bam2fastq} {
	test_cleantmp
	set bamfile genomecomb.testdata/exomes_yri_chr2122.reference/samples/NA19240chr2122/map-rdsbwa-NA19240chr2122.bam
	cg bam2fastq $bamfile tmp/out_1.fq.gz tmp/out_2.fq.gz
	exec diff tmp/out_1.fq.gz genomecomb.testdata/expected/NA19240chr2122_1.fq.gz
	exec diff tmp/out_2.fq.gz genomecomb.testdata/expected/NA19240chr2122_2.fq.gz
} {}

testsummarize
