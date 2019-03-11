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

test bam_histo {bam_histo different order regfile} {
	test_cleantmp
	cg select -f {info chromosome begin end} data/reg_hg19_smallpartexome.tsv tmp/regfile.tsv
	set namecol info
	set regionfile tmp/regfile.tsv
	set bamfile genomecomb.testdata/ori/test-map-rdsbwa-NA19240chr2122.bam
	set intervals {1 5 10 20 50 100 200 500 1000}
	cg bam_histo -n $namecol $regionfile $bamfile $intervals > tmp/result.tsv
	exec diff tmp/result.tsv genomecomb.testdata/expected/bam_histo-NA19240_smallpartchr2122.tsv
} {}

test depth_histo {depth_histo} {
	test_cleantmp
	file copy data/reg_hg19_smallpartexome.tsv tmp/regfile.tsv
	set regionfile tmp/regfile.tsv
	set bamfile genomecomb.testdata/ori/test-map-rdsbwa-NA19240chr2122.bam
	set max 1000
	cg depth_histo -max $max $bamfile $regionfile > tmp/result.tsv
	exec diff tmp/result.tsv genomecomb.testdata/expected/depth_histo-NA19240_smallpartchr2122.tsv
} {}

test depth_histo {depth_histo -Q -max} {
	test_cleantmp
	file copy data/reg_hg19_smallpartexome.tsv tmp/regfile.tsv
	set regionfile tmp/regfile.tsv
	set bamfile genomecomb.testdata/ori/test-map-rdsbwa-NA19240chr2122.bam
	set max 200
	cg depth_histo -max $max -Q 20 $bamfile $regionfile > tmp/result.tsv
	exec diff tmp/result.tsv genomecomb.testdata/expected/depth_histo_Q20-NA19240_smallpartchr2122.tsv
} {}

test depth_histo {depth_histo no regfile} {
	test_cleantmp
	file copy data/reg_hg19_smallpartexome.tsv tmp/regfile.tsv
	set regionfile {}
	set bamfile genomecomb.testdata/ori/test-map-rdsbwa-NA19240chr2122.bam
	set max 1000
	cg depth_histo -max $max $bamfile $regionfile > tmp/result.tsv
	cg select -f {depth {ontarget=0} {offtarget=$ontarget + $offtarget}} genomecomb.testdata/expected/depth_histo-NA19240_smallpartchr2122.tsv tmp/expected.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test depth_histo {depth_histo compressed regfile} {
	test_cleantmp
	compress data/reg_hg19_smallpartexome.tsv tmp/regfile.tsv.zst
	set regionfile tmp/regfile.tsv.zst
	set bamfile genomecomb.testdata/ori/test-map-rdsbwa-NA19240chr2122.bam
	set max 1000
	cg depth_histo -max $max $bamfile $regionfile > tmp/result.tsv
	exec diff tmp/result.tsv genomecomb.testdata/expected/depth_histo-NA19240_smallpartchr2122.tsv
} {}

test bam2fastq {bam2fastq} {
	test_cleantmp
	set bamfile genomecomb.testdata/ori/test-map-rdsbwa-NA19240chr2122.bam
	cg bam2fastq $bamfile tmp/out_1.fq tmp/out_2.fq
	exec gunzip -c genomecomb.testdata/expected/test-NA19240chr2122_1.fq.gz > tmp/expected1.fq
	exec gunzip -c genomecomb.testdata/expected/test-NA19240chr2122_2.fq.gz > tmp/expected2.fq
	cg sortfastq tmp/out_1.fq tmp/sout_1.fq
	cg sortfastq tmp/out_2.fq tmp/sout_2.fq
	exec diff --brief tmp/sout_1.fq tmp/expected1.fq
	exec diff --brief tmp/sout_2.fq tmp/expected2.fq
} {}

test bam2fastq {bam2fastq gz} {
	test_cleantmp
	set bamfile genomecomb.testdata/ori/test-map-rdsbwa-NA19240chr2122.bam
	cg bam2fastq $bamfile tmp/out_1.fq.gz tmp/out_2.fq.gz
	exec gunzip tmp/out_1.fq.gz tmp/out_2.fq.gz
	exec gunzip -c genomecomb.testdata/expected/test-NA19240chr2122_1.fq.gz > tmp/expected1.fq
	exec gunzip -c genomecomb.testdata/expected/test-NA19240chr2122_2.fq.gz > tmp/expected2.fq
	cg sortfastq tmp/out_1.fq tmp/sout_1.fq
	cg sortfastq tmp/out_2.fq tmp/sout_2.fq
	exec diff --brief tmp/sout_1.fq tmp/expected1.fq
	exec diff --brief tmp/sout_2.fq tmp/expected2.fq
} {}

test cg_regextract {regextract} {
	test_cleantmp
	set bamfile genomecomb.testdata/ori/test-map-rdsbwa-NA19240chr2122.bam
	cg regextract -min 20 $bamfile > tmp/reg-cov20.tsv 2>@ stderr
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
	cg select -q {$chromosome in "chr21 chr22"} [gzfile $::refseqdir/hg19/extra/reg_hg19_exome_SeqCap_EZ_v3.tsv] tmp/regfile.tsv
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
} {56.1521}

test cg_regextract {-q 20 -Q 20} {
	test_cleantmp
	cg select -q {$chromosome in "chr21 chr22"} [gzfile $::refseqdir/hg19/extra/reg_hg19_exome_SeqCap_EZ_v3.tsv] tmp/regfile.tsv
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
} {54.9979}

test cg_regextract {small sam -min 1} {
	write_sam tmp/temp.sam {
		chr1	100	20M	20	chr1	121	20M	20
		chr1	100	2M1X17M	20	chr1	121	17M1X2M	20
		chr1	100	2S18M	20	chr1	121	2S16M2S	20
		chr1	100	2S2M2D16M	20	chr1	120	2S16M2D2M	20
		chr1	100	2S2M2I14M	20	chr1	123	2S14M2I2M	20
		chr1	102	4S2M2I12M	20	chr1	111	30M	30
		chr1	107	20M	20	chr1	112	20M	20
		chr1	108	20M	20	chr1	113	20M	20
		chr2	50	20M	20	chr2	60	20M	20
		chr2	100	50M	50	chr2	100	50M	50
	}
	cg regextract -min 1 tmp/temp.sam
} {chromosome	begin	end
chr1	99	140
chr2	49	79
chr2	99	149}

test cg_bam2reg {small sam -min 2} {
	write_sam tmp/temp.sam {
		chr1	100	20M	20	chr1	121	20M	20
		chr1	102	4S2M2I12M	20	chr1	111	30M	30
		chr2	50	20M	20	chr2	60	20M	20
		chr2	100	50M	50	chr2	100	50M	50
	}
	cg bam2reg -mincoverage 2 tmp/temp.sam
	exec cg zcat tmp/sreg-cov2-temp.tsv.zst
} {chromosome	begin	end
chr1	101	119
chr1	120	140
chr2	59	69
chr2	99	149}

test cg_bam2reg {small sam -min 2 with target} {
	write_sam tmp/temp.sam {
		chr1	100	20M	20	chr1	121	20M	20
		chr1	102	4S2M2I12M	20	chr1	111	30M	30
		chr2	50	20M	20	chr2	60	20M	20
		chr2	100	50M	50	chr2	100	50M	50
	}
	cg bam2reg -mincoverage 2 tmp/temp.sam tmp/reg.tsv
	file_read tmp/reg.tsv
} {chromosome	begin	end
chr1	101	119
chr1	120	140
chr2	59	69
chr2	99	149
}

test cg_bam2reg {small sam 2} {
	write_sam tmp/temp.sam {
		chr1	100	20M	20	chr1	121	20M	20
		chr1	102	4S2M2I12M	20	chr1	111	30M	30
		chr2	50	20M	20	chr2	60	20M	20
		chr2	100	50M	50	chr2	100	50M	50
	}
	cg bam2reg tmp/temp.sam 2
	exec cg zcat tmp/sreg-cov2-temp.tsv.zst
} {chromosome	begin	end
chr1	101	119
chr1	120	140
chr2	59	69
chr2	99	149}

test cg_bam2reg {small sam 2 target} {
	write_sam tmp/temp.sam {
		chr1	100	20M	20	chr1	121	20M	20
		chr1	102	4S2M2I12M	20	chr1	111	30M	30
		chr2	50	20M	20	chr2	60	20M	20
		chr2	100	50M	50	chr2	100	50M	50
	}
	cg bam2reg tmp/temp.sam 2 tmp/reg.tsv
	file_read tmp/reg.tsv
} {chromosome	begin	end
chr1	101	119
chr1	120	140
chr2	59	69
chr2	99	149
}

test cg_bam2reg {basic} {
	mklink genomecomb.testdata/ori/test-map-rdsbwa-NA19240chr2122.bam tmp/test.bam
	cg bam2reg -stack 1 -mincoverage 20 tmp/test.bam tmp/result.tsv
	catch {exec cg zcat tmp/result.tsv | head -4} temp
	append temp \n[exec cg zcat tmp/result.tsv | tail -3]
} {chromosome	begin	end
chr21	9439337	9439639
chr21	9439640	9439747
chr21	9444036	9444482
chr22	51229714	51229749
chr22	51229752	51229787
chr22	51237064	51237741}

testsummarize
