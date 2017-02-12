#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test bam_histo {bam_histo} {
	test_cleantmp
	file copy data/reg_hg19_smallpartexome.tsv tmp/regfile.tsv
	set namecol info
	set regionfile tmp/regfile.tsv
	set bamfile genomecomb.testdata/ori/test-map-rdsbwa-NA19240chr2122.bam
	set intervals {1 5 10 20 50 100 200 500 1000}
	cg bam_histo -n $namecol $regionfile $bamfile $intervals > tmp/result.tsv
	exec diff tmp/result.tsv genomecomb.testdata/expected/bam_histo-NA19240_smallpartchr2122.tsv
} {}

test bam_histo {bam_histo chr_clipped} {
	test_cleantmp
	cg select -f {chromosome=chr_clip($chromosome) begin end info} data/reg_hg19_smallpartexome.tsv tmp/regfile.tsv
	set regionfile tmp/regfile.tsv
	set bamfile genomecomb.testdata/ori/test-map-rdsbwa-NA19240chr2122.bam
	set intervals {1 5 10 20 50 100 200 500 1000}
	cg bam_histo -n info $regionfile $bamfile $intervals > tmp/result.tsv
	exec diff tmp/result.tsv genomecomb.testdata/expected/bam_histo-NA19240_smallpartchr2122.tsv
} {}

test bam2fastq {bam2fastq} {
	test_cleantmp
	set bamfile genomecomb.testdata/ori/test-map-rdsbwa-NA19240chr2122.bam
	cg bam2fastq $bamfile tmp/out_1.fq tmp/out_2.fq
	exec gunzip -c genomecomb.testdata/expected/test-NA19240chr2122_1.fq.gz > tmp/expected1.fq
	exec gunzip -c genomecomb.testdata/expected/test-NA19240chr2122_2.fq.gz > tmp/expected2.fq
	exec diff tmp/out_1.fq tmp/expected1.fq
	exec diff tmp/out_2.fq tmp/expected2.fq
} {}

test bam2fastq {bam2fastq gz} {
	test_cleantmp
	set bamfile genomecomb.testdata/ori/test-map-rdsbwa-NA19240chr2122.bam
	cg bam2fastq $bamfile tmp/out_1.fq.gz tmp/out_2.fq.gz
	exec gunzip tmp/out_1.fq.gz tmp/out_2.fq.gz
	exec gunzip -c genomecomb.testdata/expected/test-NA19240chr2122_1.fq.gz > tmp/expected1.fq
	exec gunzip -c genomecomb.testdata/expected/test-NA19240chr2122_2.fq.gz > tmp/expected2.fq
	exec diff tmp/out_1.fq tmp/expected1.fq
	exec diff tmp/out_2.fq tmp/expected2.fq
} {}

test cg_regextract {regextract} {
	test_cleantmp
	set bamfile genomecomb.testdata/ori/test-map-rdsbwa-NA19240chr2122.bam
	cg regextract -min 20 $bamfile > tmp/reg-cov20.tsv
	string_change [cg covered tmp/reg-cov20.tsv] [list \n\n \n]
} {chromosome	bases
chr21	495709
chr22	1168356
total	1664065}

test cg_regextract {regextract -filtered 1 } {
	test_cleantmp
	set bamfile genomecomb.testdata/ori/test-map-rdsbwa-NA19240chr2122.bam
	cg regextract --filtered 1 -min 20 $bamfile > tmp/reg-cov20.tsv
	string_change [cg covered tmp/reg-cov20.tsv] [list \n\n \n]
} {chromosome	bases
chr21	464015
chr22	1004396
total	1468411}

test cg_regextract {--filtered 1 -q 1 -Q 0} {
	test_cleantmp
	cg select -q {$chromosome in "chr21 chr22"} $::refseqdir/hg19/extra/reg_hg19_exome_SeqCap_EZ_v3.tsv tmp/regfile.tsv
	set bamfile genomecomb.testdata/ori/test-map-rdsbwa-NA19240chr2122.bam
	# TARGET_TERRITORY 2126556
	# PCT_TARGET_BASES_20X	0.56393
	# reg: 1199228.72508
	# next settings come closest to hsmetrics?
	cg regextract --filtered 1 -q 1 -Q 0 -min 20 $bamfile > tmp/reg-cov20.tsv
	set regionfile tmp/regfile.tsv
	file delete tmp/multireg.tsv
	cg multireg tmp/multireg.tsv tmp/regfile.tsv tmp/reg-cov20.tsv
	cg select -q {$regfile == 1} tmp/multireg.tsv tmp/selmultireg.tsv
	set total [lindex [cg covered tmp/selmultireg.tsv] end]
	set cov20 [lindex [exec cg select -q {$reg-cov20 == 1} tmp/selmultireg.tsv | cg covered] end]
	format %.4f [expr 100.0*$cov20/$total]
} {56.1604}

test cg_regextract {-q 20 -Q 20} {
	test_cleantmp
	cg select -q {$chromosome in "chr21 chr22"} $::refseqdir/hg19/extra/reg_hg19_exome_SeqCap_EZ_v3.tsv tmp/regfile.tsv
	set bamfile genomecomb.testdata/ori/test-map-rdsbwa-NA19240chr2122.bam
	# TARGET_TERRITORY 2126556
	# PCT_TARGET_BASES_20X	0.56393
	# reg: 1199228.72508
	# next settings come closest to hsmetrics?
	cg regextract -q 20 -Q 20 -min 20 $bamfile > tmp/reg-cov20.tsv
	set regionfile tmp/regfile.tsv
	file delete tmp/multireg.tsv
	cg multireg tmp/multireg.tsv tmp/regfile.tsv tmp/reg-cov20.tsv
	cg select -q {$regfile == 1} tmp/multireg.tsv tmp/selmultireg.tsv
	set total [lindex [cg covered tmp/selmultireg.tsv] end]
	set cov20 [lindex [exec cg select -q {$reg-cov20 == 1} tmp/selmultireg.tsv | cg covered] end]
	format %.4f [expr 100.0*$cov20/$total]
} {55.0057}

testsummarize
