#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test bam_histo {bam_histo} {
	test_cleantmp
	file copy data/reg_hg19_smallpartexome.tsv tmp/regfile.tsv
	set namecol info
	set regionfile tmp/regfile.tsv
	set bamfile $::smalltestdir/ori/test-map-rdsbwa-NA19240chr2122.bam
	set intervals {1 5 10 20 50 100 200 500 1000}
	cg bam_histo -n $namecol $regionfile $bamfile $intervals > tmp/result.tsv
	exec diff tmp/result.tsv $::smalltestdir/expected/bam_histo-NA19240_smallpartchr2122.tsv
} {}

test bam_histo {bam_histo chr_clipped} {
	test_cleantmp
	cg select -f {chromosome=chr_clip($chromosome) begin end info} data/reg_hg19_smallpartexome.tsv tmp/regfile.tsv
	set regionfile tmp/regfile.tsv
	set bamfile $::smalltestdir/ori/test-map-rdsbwa-NA19240chr2122.bam
	set intervals {1 5 10 20 50 100 200 500 1000}
	cg bam_histo -n info $regionfile $bamfile $intervals > tmp/result.tsv
	exec diff tmp/result.tsv $::smalltestdir/expected/bam_histo-NA19240_smallpartchr2122.tsv
} {}

test bam_histo {bam_histo different order regfile} {
	test_cleantmp
	cg select -f {info chromosome begin end} data/reg_hg19_smallpartexome.tsv tmp/regfile.tsv
	set namecol info
	set regionfile tmp/regfile.tsv
	set bamfile $::smalltestdir/ori/test-map-rdsbwa-NA19240chr2122.bam
	set intervals {1 5 10 20 50 100 200 500 1000}
	cg bam_histo -n $namecol $regionfile $bamfile $intervals > tmp/result.tsv
	exec diff tmp/result.tsv $::smalltestdir/expected/bam_histo-NA19240_smallpartchr2122.tsv
} {}

test depth_histo {depth_histo} {
	test_cleantmp
	file copy data/reg_hg19_smallpartexome.tsv tmp/regfile.tsv
	set regionfile tmp/regfile.tsv
	set bamfile $::smalltestdir/ori/test-map-rdsbwa-NA19240chr2122.bam
	set max 1000
	cg depth_histo -max $max $bamfile $regionfile > tmp/result.tsv
	exec diff tmp/result.tsv $::smalltestdir/expected/depth_histo-NA19240_smallpartchr2122.tsv
} {}

test depth_histo {depth_histo -Q -max} {
	test_cleantmp
	file copy data/reg_hg19_smallpartexome.tsv tmp/regfile.tsv
	set regionfile tmp/regfile.tsv
	set bamfile $::smalltestdir/ori/test-map-rdsbwa-NA19240chr2122.bam
	set max 200
	cg depth_histo -max $max -Q 20 $bamfile $regionfile > tmp/result.tsv
	exec diff tmp/result.tsv $::smalltestdir/expected/depth_histo_Q20-NA19240_smallpartchr2122.tsv
} {}

test depth_histo {depth_histo no regfile} {
	test_cleantmp
	file copy data/reg_hg19_smallpartexome.tsv tmp/regfile.tsv
	set regionfile {}
	set bamfile $::smalltestdir/ori/test-map-rdsbwa-NA19240chr2122.bam
	set max 1000
	cg depth_histo -max $max $bamfile $regionfile > tmp/result.tsv
	cg select -f {depth {ontarget=0} {offtarget=$ontarget + $offtarget}} $::smalltestdir/expected/depth_histo-NA19240_smallpartchr2122.tsv tmp/expected.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test depth_histo {depth_histo compressed regfile} {
	test_cleantmp
	compress data/reg_hg19_smallpartexome.tsv tmp/regfile.tsv.zst
	set regionfile tmp/regfile.tsv.zst
	set bamfile $::smalltestdir/ori/test-map-rdsbwa-NA19240chr2122.bam
	set max 1000
	cg depth_histo -max $max $bamfile $regionfile > tmp/result.tsv
	exec diff tmp/result.tsv $::smalltestdir/expected/depth_histo-NA19240_smallpartchr2122.tsv
} {}

test bam2fastq {bam2fastq} {
	test_cleantmp
	set bamfile $::smalltestdir/ori/test-map-rdsbwa-NA19240chr2122.bam
	cg bam2fastq $bamfile tmp/out_1.fq tmp/out_2.fq
	cg zcat $::smalltestdir/expected/bam/test-NA19240chr2122_1.fq.zst > tmp/expected1.fq
	cg zcat $::smalltestdir/expected/bam/test-NA19240chr2122_2.fq.zst > tmp/expected2.fq
	cg sortfastq tmp/out_1.fq tmp/sout_1.fq
	cg sortfastq tmp/out_2.fq tmp/sout_2.fq
	exec diff --brief tmp/sout_1.fq tmp/expected1.fq
	exec diff --brief tmp/sout_2.fq tmp/expected2.fq
} {}

test bam2fastq {bam2fastq gz} {
	test_cleantmp
	set bamfile $::smalltestdir/ori/test-map-rdsbwa-NA19240chr2122.bam
	cg bam2fastq $bamfile tmp/out_1.fq.gz tmp/out_2.fq.gz
	exec gunzip tmp/out_1.fq.gz tmp/out_2.fq.gz
	cg zcat $::smalltestdir/expected/bam/test-NA19240chr2122_1.fq.zst > tmp/expected1.fq
	cg zcat $::smalltestdir/expected/bam/test-NA19240chr2122_2.fq.zst > tmp/expected2.fq
	cg sortfastq tmp/out_1.fq tmp/sout_1.fq
	cg sortfastq tmp/out_2.fq tmp/sout_2.fq
	exec diff --brief tmp/sout_1.fq tmp/expected1.fq
	exec diff --brief tmp/sout_2.fq tmp/expected2.fq
} {}

test bam2fastq {bam2fastq -sortmethod collate gz} {
	test_cleantmp
	set bamfile $::smalltestdir/ori/test-map-rdsbwa-NA19240chr2122.bam
	cg bam2fastq -sortmethod collate $bamfile tmp/out_1.fq.gz tmp/out_2.fq.gz
	exec gunzip tmp/out_1.fq.gz tmp/out_2.fq.gz
	cg zcat $::smalltestdir/expected/bam/test-NA19240chr2122_1.fq.zst > tmp/expected1.fq
	cg zcat $::smalltestdir/expected/bam/test-NA19240chr2122_2.fq.zst > tmp/expected2.fq
	cg sortfastq tmp/out_1.fq tmp/sout_1.fq
	cg sortfastq tmp/out_2.fq tmp/sout_2.fq
	exec diff --brief tmp/sout_1.fq tmp/expected1.fq
	exec diff --brief tmp/sout_2.fq tmp/expected2.fq
} {}

test cg_regextract {regextract basic} {
	test_cleantmp
	set bamfile data/bwa.sam
	cg regextract -min 2 $bamfile > tmp/reg-cov20.tsv 2>@ stderr
	string_change [cg covered tmp/reg-cov20.tsv] [list \n\n \n]
} {chromosome	bases
chr21	3376
chr22	1017
total	4393}

test cg_regextract {regextract -region} {
	test_cleantmp
	exec samtools view --no-PG -h -b data/bwa.sam > tmp/bwa.bam
	exec samtools index tmp/bwa.bam
	set bamfile tmp/bwa.bam
	cg regextract -min 2 -region chr22 $bamfile > tmp/reg-cov20.tsv 2>@ stderr
	string_change [cg covered tmp/reg-cov20.tsv] [list \n\n \n]
} {chromosome	bases
chr22	1017
total	1017}

test cg_regextract {regextract -region chr21_random1} {
	test_cleantmp
	set c [file_read data/bwa.sam]
	set c [string_change $c {chr21 chr21_random1 chr22 chr21_random2}]
	file_write tmp/bwa.sam $c
	exec samtools view --no-PG -h -b tmp/bwa.sam > tmp/bwa.bam
	exec samtools index tmp/bwa.bam
	set bamfile tmp/bwa.bam
	# dummy reference sequence, only needed here for the fai
	file_write tmp/ref.fas ""
	file_write tmp/ref.fas.fai [deindent {
		chr21_random1	48129895	2784007886	48129895	48129896
		chr21_random2	51304566	2832165523	51304566	51304567
		chr22	51304566	2832165523	51304566	51304567
	}]
	cg regextract -min 2 -refseq tmp/ref.fas -region chr21_ $bamfile > tmp/reg-cov20.tsv 2>@ stderr
	string_change [cg covered tmp/reg-cov20.tsv] [list \n\n \n]
} {chromosome	bases
chr21_random1	3376
chr21_random2	1017
total	4393}

test cg_regextract {regextract} {
	test_cleantmp
	set bamfile $::smalltestdir/ori/test-map-rdsbwa-NA19240chr2122.bam
	cg regextract -min 20 $bamfile > tmp/reg-cov20.tsv 2>@ stderr
	string_change [cg covered tmp/reg-cov20.tsv] [list \n\n \n]
} {chromosome	bases
chr21	495709
chr22	1168356
total	1664065}

test cg_regextract {regextract -region 2} {
	test_cleantmp
	set bamfile $::smalltestdir/ori/test-map-rdsbwa-NA19240chr2122.bam
	cg regextract -min 20 -region chr21 $bamfile > tmp/reg-cov20.tsv 2>@ stderr
	string_change [cg covered tmp/reg-cov20.tsv] [list \n\n \n]
} {chromosome	bases
chr21	495709
total	495709}

test cg_regextract {regextract -filtered 1 } {
	test_cleantmp
	set bamfile $::smalltestdir/ori/test-map-rdsbwa-NA19240chr2122.bam
	cg regextract --filtered 1 -min 20 $bamfile > tmp/reg-cov20.tsv
	string_change [cg covered tmp/reg-cov20.tsv] [list \n\n \n]
} {chromosome	bases
chr21	464015
chr22	1004396
total	1468411}

test cg_regextract {--filtered 1 -q 1 -Q 0} {
	test_cleantmp
	cg select -q {$chromosome in "chr21 chr22"} [gzfile $::refseqdir/hg19/extra/reg_hg19_exome_SeqCap_EZ_v3.tsv] tmp/regfile.tsv
	set bamfile $::smalltestdir/ori/test-map-rdsbwa-NA19240chr2122.bam
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
	set bamfile $::smalltestdir/ori/test-map-rdsbwa-NA19240chr2122.bam
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

test cg_bam2reg {small bam -min 2 distrreg} {
	write_sam tmp/temp.sam {
		chr1	100	20M	20	chr1	121	20M	20
		chr1	102	4S2M2I12M	20	chr1	111	30M	30
		chr2	50	20M	20	chr2	60	20M	20
		chr2	100	50M	50	chr2	100	50M	50
	}
	exec samtools view --no-PG -b tmp/temp.sam > tmp/temp.bam
	exec samtools index tmp/temp.bam
	cg bam2reg -distrreg chr -refseq $::refseqdir/hg19/genome_hg19.ifas -mincoverage 2 tmp/temp.bam
	exec cg zcat tmp/sreg-cov2-temp.tsv.zst
} {chromosome	begin	end
chr1	101	119
chr1	120	140
chr2	59	69
chr2	99	149}

test cg_bam2reg {basic} {
	mklink $::smalltestdir/ori/test-map-rdsbwa-NA19240chr2122.bam tmp/test.bam
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

test cg_bamreorder {basic genomecomb} {
	file_write tmp/ref.fai [deindent {
		chr22	51304566	2832165523	51304566	51304567
		chr21_gl000210_random	27682	2832137827	27682	27683
		chr21	48129895	2784007886	48129895	48129896
	}]
	file delete tmp/result.bam  tmp/result.bam.bai
	catch_exec cg bamreorder data/bwa.bam tmp/result.bam tmp/ref
	list [exec md5sum tmp/result.bam] [exec samtools view --no-PG -H tmp/result.bam]
} {{7bcf2370e55bb6f10a9794b9e76fc567  tmp/result.bam} {@HD	VN:1.3	SO:coordinate
@SQ	SN:chr22	LN:51304566
@SQ	SN:chr21_gl000210_random	LN:27682
@SQ	SN:chr21	LN:48129895
@RG	ID:NA19240m	SM:NA19240m	PL:illumina	PU:NA19240m	LB:solexa-123
@PG	ID:bwa	PN:bwa	VN:0.7.15-r1140	CL:bwa mem -t 2 -M -R @RG\tID:NA19240m\tSM:NA19240m\tPL:illumina\tPU:NA19240m\tLB:solexa-123 /home/peter/dev/genomecomb/tests/genomecomb.testdata/refseqtest/hg19/genome_hg19.ifas.bwa/genome_hg19.ifas /home/peter/dev/genomecomb/tests/tmp/seq_1.fastq /home/peter/dev/genomecomb/tests/tmp/seq_2.fastq}}

test cg_sam_catmerge {basic cg_sam_catmerge} {
	write_sam tmp/temp1.sam {
		chr1	102	4S2M2I12M	20	chr1	111	30M	30
		chr1	107	20M	20	chr1	112	20M	20
		chr1	108	20M	20	chr1	113	20M	20
		chr1	100	20M	20	chr1	121	20M	20
	}
	write_sam tmp/temp2.sam {
		chr2	50	20M	20	chr2	60	20M	20
		chr2	100	50M	50	chr2	100	50M	50
	}
	file_write tmp/expected.sam [deindent $::sam_header]\n[deindent {
		A4	99	chr1	100	60	20M	=	121	41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A1	99	chr1	102	60	4S2M2I12M	=	111	39	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A2	99	chr1	107	60	20M	=	112	25	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A3	99	chr1	108	60	20M	=	113	25	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A1	147	chr1	111	60	30M	=	102	-39	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A2	147	chr1	112	60	20M	=	107	-25	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A3	147	chr1	113	60	20M	=	108	-25	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A4	147	chr1	121	60	20M	=	100	-41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A1	99	chr2	50	60	20M	=	60	30	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A1	147	chr2	60	60	20M	=	50	-30	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A2	99	chr2	100	60	50M	=	100	50	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	--------------------------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A2	147	chr2	100	60	50M	=	100	-50	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	--------------------------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
	}]\n
	file delete tmp/result.sam
	cg sam_catmerge -stack 1 -v 2 tmp/result.sam tmp/temp1.sam tmp/temp2.sam
	exec diff -I {@HD	} tmp/result.sam tmp/expected.sam
} {}

test cg_sam_catmerge {cg_sam_catmerge nosort} {
	file_write tmp/temp1.sam [deindent $::sam_header]\n[deindent {
		A4	147	chr1	121	60	20M	=	100	-41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A4	99	chr1	100	60	20M	=	121	41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
	}]\n
	write_sam tmp/temp2.sam {
		chr2	50	20M	20	chr2	60	20M	20
		chr2	100	50M	50	chr2	100	50M	50
	}
	file_write tmp/expected.sam [deindent $::sam_header]\n[deindent {
		A4	147	chr1	121	60	20M	=	100	-41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A4	99	chr1	100	60	20M	=	121	41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A1	99	chr2	50	60	20M	=	60	30	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A1	147	chr2	60	60	20M	=	50	-30	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A2	99	chr2	100	60	50M	=	100	50	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	--------------------------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A2	147	chr2	100	60	50M	=	100	-50	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	--------------------------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
	}]\n
	cg sam_catmerge -stack 1 -sort nosort tmp/result.sam tmp/temp1.sam tmp/temp2.sam
	exec diff tmp/result.sam tmp/expected.sam
} {}

test cg_sam_catmerge {cg_sam_catmerge nosort compressed sams} {
	file_write tmp/temp1.sam [deindent $::sam_header]\n[deindent {
		A4	147	chr1	121	60	20M	=	100	-41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A4	99	chr1	100	60	20M	=	121	41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
	}]\n
	write_sam tmp/temp2.sam {
		chr2	50	20M	20	chr2	60	20M	20
		chr2	100	50M	50	chr2	100	50M	50
	}
	file_write tmp/expected.sam [deindent $::sam_header]\n[deindent {
		A4	147	chr1	121	60	20M	=	100	-41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A4	99	chr1	100	60	20M	=	121	41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A1	99	chr2	50	60	20M	=	60	30	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A1	147	chr2	60	60	20M	=	50	-30	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A2	99	chr2	100	60	50M	=	100	50	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	--------------------------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A2	147	chr2	100	60	50M	=	100	-50	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	--------------------------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
	}]\n
	cg zst tmp/temp1.sam
	cg zst tmp/temp2.sam
	cg sam_catmerge -stack 1 -sort nosort tmp/result.sam tmp/temp1.sam.zst tmp/temp2.sam.zst
	exec diff tmp/result.sam tmp/expected.sam
} {}

test cg_sam_catmerge {cg_sam_catmerge to cram} {
	file_write tmp/temp1.sam [deindent $::sam_header]\n[deindent {
		A4	147	chr1	121	60	20M	=	100	-41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A4	99	chr1	100	60	20M	=	121	41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
	}]\n
	write_sam tmp/temp2.sam {
		chr2	50	20M	20	chr2	60	20M	20
		chr2	100	50M	50	chr2	100	50M	50
	}
	file_write tmp/expected.sam [deindent $::sam_header]\n[deindent {
		A4	99	chr1	100	60	20M	=	121	41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A4	147	chr1	121	60	20M	=	100	-41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A1	99	chr2	50	60	20M	=	60	30	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A1	147	chr2	60	60	20M	=	50	-30	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A2	99	chr2	100	60	50M	=	100	50	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	--------------------------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A2	147	chr2	100	60	50M	=	100	-50	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	--------------------------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
	}]\n
	cg sam_catmerge -stack 1 -refseq $::refseqdir/hg19/genome_hg19.ifas tmp/result.cram tmp/temp1.sam tmp/temp2.sam
	cramdiff tmp/result.cram tmp/expected.sam
} {}

test cg_sam_catmerge {cg_sam_catmerge to compressed sam} {
	file_write tmp/temp1.sam [deindent $::sam_header]\n[deindent {
		A4	147	chr1	121	60	20M	=	100	-41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A4	99	chr1	100	60	20M	=	121	41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
	}]\n
	write_sam tmp/temp2.sam {
		chr2	50	20M	20	chr2	60	20M	20
		chr2	100	50M	50	chr2	100	50M	50
	}
	file_write tmp/expected.sam [deindent $::sam_header]\n[deindent {
		A4	99	chr1	100	60	20M	=	121	41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A4	147	chr1	121	60	20M	=	100	-41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A1	99	chr2	50	60	20M	=	60	30	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A1	147	chr2	60	60	20M	=	50	-30	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A2	99	chr2	100	60	50M	=	100	50	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	--------------------------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A2	147	chr2	100	60	50M	=	100	-50	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	--------------------------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
	}]\n
	cg sam_catmerge -stack 1 tmp/result.sam.zst tmp/temp1.sam tmp/temp2.sam
	cramdiff tmp/result.sam.zst tmp/expected.sam
} {}

test cg_sam_catmerge {cg_sam_catmerge nosort bam} {
	file_write tmp/temp1.sam [deindent $::sam_header]\n[deindent {
		A4	147	chr1	121	60	20M	=	100	-41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A4	99	chr1	100	60	20M	=	121	41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
	}]\n
	write_sam tmp/temp2.sam {
		chr2	50	20M	20	chr2	60	20M	20
		chr2	100	50M	50	chr2	100	50M	50
	}
	file_write tmp/expected.sam [deindent $::sam_header]\n[deindent {
		A4	147	chr1	121	60	20M	=	100	-41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A4	99	chr1	100	60	20M	=	121	41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A1	99	chr2	50	60	20M	=	60	30	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A1	147	chr2	60	60	20M	=	50	-30	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A2	99	chr2	100	60	50M	=	100	50	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	--------------------------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A2	147	chr2	100	60	50M	=	100	-50	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	--------------------------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
	}]\n
	cg sam_catmerge -stack 1 -sort nosort -index 0 tmp/result.bam tmp/temp1.sam tmp/temp2.sam
	exec samtools view --no-PG -h tmp/result.bam -o tmp/result.sam
	exec diff tmp/result.sam tmp/expected.sam
} {}

test cg_sam_catmerge {-mergesort 1 -distrreg 1 with unaligned} {
	write_sam tmp/utemp1.sam {
		chr1	1000100	20M	20	chr1	1000121	20M	20
		chr1	1000107	20M	20	chr1	1000112	20M	20
		chr2	1000100	50M	50	chr2	1000100	50M	50
	}
	cg sam_sort tmp/utemp1.sam tmp/temp1.sam
	write_sam tmp/utemp2.sam {
		chr1	1000102	4S2M2I12M	20	chr1	1000111	30M	30
		chr1	1000108	20M	20	chr1	1000113	20M	20
		chr2	1000050	20M	20	chr2	1000060	20M	20
		*	0	*	20	*	0	*	20
	} B
	cg sam_sort tmp/utemp2.sam tmp/temp2.sam
	file_write tmp/expected-chr1.sam [deindent $::sam_header]\n[deindent {
		A1	99	chr1	1000100	60	20M	=	1000121	41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		B1	99	chr1	1000102	60	4S2M2I12M	=	1000111	39	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A2	99	chr1	1000107	60	20M	=	1000112	25	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		B2	99	chr1	1000108	60	20M	=	1000113	25	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		B1	147	chr1	1000111	60	30M	=	1000102	-39	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A2	147	chr1	1000112	60	20M	=	1000107	-25	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		B2	147	chr1	1000113	60	20M	=	1000108	-25	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A1	147	chr1	1000121	60	20M	=	1000100	-41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
	}]\n
	file_write tmp/expected-chr2.sam [deindent $::sam_header]\n[deindent {
		B3	99	chr2	1000050	60	20M	=	1000060	30	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		B3	147	chr2	1000060	60	20M	=	1000050	-30	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A3	99	chr2	1000100	60	50M	=	1000100	50	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	--------------------------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A3	147	chr2	1000100	60	50M	=	1000100	-50	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	--------------------------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
	}]\n
	file_write tmp/expected-unaligned.sam [deindent $::sam_header]\n[deindent {
		B4	103	*	0	60	*	*	0	20	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		B4	151	*	0	60	*	*	0	-20	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
	}]\n
	cg sam_catmerge -stack 1 -mergesort 1 -distrreg 1 \
		-refseq $::refseqdir/hg19/genome_hg19.ifas \
		tmp/result.sam tmp/temp1.sam tmp/temp2.sam
	exec diff -I {@HD	} tmp/result-chr1-0-29953082.sam tmp/expected-chr1.sam
	exec diff -I {@HD	} tmp/result-chr2-0-33142192.sam tmp/expected-chr2.sam
	exec diff -I {@HD	} tmp/result-unaligned.sam tmp/expected-unaligned.sam
} {}

test cg_sam_catmerge {-mergesort 1 -distrreg 1 .sam.zst} {
	write_sam tmp/utemp1.sam {
		chr1	1000107	20M	20	chr1	1000112	20M	20
		chr1	1000100	20M	20	chr1	1000121	20M	20
		chr2	1000100	50M	50	chr2	1000100	50M	50
	}
	cg sam_sort tmp/utemp1.sam tmp/temp1.sam
	write_sam tmp/utemp2.sam {
		chr1	1000102	4S2M2I12M	20	chr1	1000111	30M	30
		chr1	1000108	20M	20	chr1	1000113	20M	20
		chr2	1000050	20M	20	chr2	1000060	20M	20
	}
	cg sam_sort tmp/utemp2.sam tmp/temp2.sam
	file_write tmp/expected-chr1_start.sam [deindent $::sam_header]\n[deindent {
		A2	99	chr1	1000100	60	20M	=	1000121	41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A1	99	chr1	1000102	60	4S2M2I12M	=	1000111	39	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A1	99	chr1	1000107	60	20M	=	1000112	25	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A2	99	chr1	1000108	60	20M	=	1000113	25	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A1	147	chr1	1000111	60	30M	=	1000102	-39	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A1	147	chr1	1000112	60	20M	=	1000107	-25	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A2	147	chr1	1000113	60	20M	=	1000108	-25	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A2	147	chr1	1000121	60	20M	=	1000100	-41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
	}]\n
	file_write tmp/expected-chr2_start.sam [deindent $::sam_header]\n[deindent {
		A3	99	chr2	1000050	60	20M	=	1000060	30	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A3	147	chr2	1000060	60	20M	=	1000050	-30	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A3	99	chr2	1000100	60	50M	=	1000100	50	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	--------------------------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A3	147	chr2	1000100	60	50M	=	1000100	-50	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	--------------------------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
	}]\n
	cg sam_catmerge -stack 1 -sort coordinate -mergesort 1 -distrreg 1 \
		-refseq $::refseqdir/hg19/genome_hg19.ifas \
		tmp/result.sam.zst tmp/temp1.sam tmp/temp2.sam
	cg zcat tmp/result-chr1-0-29953082.sam.zst > tmp/test-chr1_start.sam
	cg zcat tmp/result-chr2-0-33142192.sam.zst > tmp/test-chr2_start.sam
	exec diff -I {@HD	} tmp/test-chr1_start.sam tmp/expected-chr1_start.sam
	exec diff -I {@HD	} tmp/test-chr2_start.sam tmp/expected-chr2_start.sam
} {}

test cg_sam_catmerge {-mergesort 1 -distrreg 1 .sam 4 files} {
	write_sam tmp/utemp1.sam {
		chr1	1000107	20M	20	chr1	1000112	20M	20
		chr2	1000100	50M	50	chr2	1000100	50M	50
	} A
	cg sam_sort tmp/utemp1.sam tmp/temp1.sam
	write_sam tmp/utemp2.sam {
		chr1	1000102	4S2M2I12M	20	chr1	1000111	30M	30
		chr1	1000108	20M	20	chr1	1000113	20M	20
	} B
	cg sam_sort tmp/utemp2.sam tmp/temp2.sam
	write_sam tmp/utemp3.sam {
		chr1	1000100	20M	20	chr1	1000121	20M	20
	} C
	cg sam_sort tmp/utemp3.sam tmp/temp3.sam
	write_sam tmp/utemp4.sam {
		chr2	1000050	20M	20	chr2	1000060	20M	20
	} D
	cg sam_sort tmp/utemp4.sam tmp/temp4.sam
	file_write tmp/expected-chr1.sam [deindent $::sam_header]\n[deindent {
		C1	99	chr1	1000100	60	20M	=	1000121	41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		B1	99	chr1	1000102	60	4S2M2I12M	=	1000111	39	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A1	99	chr1	1000107	60	20M	=	1000112	25	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		B2	99	chr1	1000108	60	20M	=	1000113	25	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		B1	147	chr1	1000111	60	30M	=	1000102	-39	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A1	147	chr1	1000112	60	20M	=	1000107	-25	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		B2	147	chr1	1000113	60	20M	=	1000108	-25	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		C1	147	chr1	1000121	60	20M	=	1000100	-41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
	}]\n
	file_write tmp/expected-chr2.sam [deindent $::sam_header]\n[deindent {
		D1	99	chr2	1000050	60	20M	=	1000060	30	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		D1	147	chr2	1000060	60	20M	=	1000050	-30	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A2	99	chr2	1000100	60	50M	=	1000100	50	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	--------------------------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A2	147	chr2	1000100	60	50M	=	1000100	-50	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	--------------------------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
	}]\n
	file delete tmp/test-chr1.sam tmp/test-chr2.sam
	cg sam_catmerge -stack 1 -mergesort 1 -distrreg chr \
		-refseq $::refseqdir/hg19/genome_hg19.ifas \
		tmp/result.sam tmp/temp1.sam tmp/temp2.sam tmp/temp3.sam tmp/temp4.sam
	exec diff -I {@HD	} tmp/result-chr1.sam tmp/expected-chr1.sam
	exec diff -I {@HD	} tmp/result-chr2.sam tmp/expected-chr2.sam
} {}

test cg_sam_catmerge {-mergesort 1 -maxopenfiles 2 .sam 5 files } {
	write_sam tmp/utemp1.sam {
		chr1	1000107	20M	20	chr1	1000112	20M	20
		chr2	1000100	50M	50	chr2	1000100	50M	50
	} A
	cg sam_sort tmp/utemp1.sam tmp/temp1.sam
	write_sam tmp/utemp2.sam {
		chr1	1000102	4S2M2I12M	20	chr1	1000111	30M	30
		chr1	1000108	20M	20	chr1	1000113	20M	20
	} B
	cg sam_sort tmp/utemp2.sam tmp/temp2.sam
	write_sam tmp/utemp3.sam {
		chr1	1000100	20M	20	chr1	1000121	20M	20
	} C
	cg sam_sort tmp/utemp3.sam tmp/temp3.sam
	write_sam tmp/utemp4.sam {
		chr2	1000050	20M	20	chr2	1000060	20M	20
	} D
	cg sam_sort tmp/utemp4.sam tmp/temp4.sam
	write_sam tmp/utemp5.sam {
		chr2	1000050	20M	20	chr2	1000060	20M	20
	} E
	cg sam_sort tmp/utemp5.sam tmp/temp5.sam
	file_write tmp/expected.sam [deindent $::sam_header]\n[deindent {
		C1	99	chr1	1000100	60	20M	=	1000121	41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		B1	99	chr1	1000102	60	4S2M2I12M	=	1000111	39	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A1	99	chr1	1000107	60	20M	=	1000112	25	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		B2	99	chr1	1000108	60	20M	=	1000113	25	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		B1	147	chr1	1000111	60	30M	=	1000102	-39	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A1	147	chr1	1000112	60	20M	=	1000107	-25	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		B2	147	chr1	1000113	60	20M	=	1000108	-25	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		C1	147	chr1	1000121	60	20M	=	1000100	-41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		D1	99	chr2	1000050	60	20M	=	1000060	30	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		E1	99	chr2	1000050	60	20M	=	1000060	30	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		D1	147	chr2	1000060	60	20M	=	1000050	-30	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		E1	147	chr2	1000060	60	20M	=	1000050	-30	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A2	99	chr2	1000100	60	50M	=	1000100	50	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	--------------------------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A2	147	chr2	1000100	60	50M	=	1000100	-50	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	--------------------------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
	}]\n
	file delete tmp/result.sam
	cg sam_catmerge -stack 1 -mergesort 1 -maxopenfiles 2 -deletesams 1 \
		-refseq $::refseqdir/hg19/genome_hg19.ifas \
		tmp/result.sam tmp/temp1.sam tmp/temp2.sam tmp/temp3.sam tmp/temp4.sam tmp/temp5.sam
	exec diff -I {@HD	} tmp/result.sam tmp/expected.sam
} {}

test cg_bam2sam {basic sortchr} {
	file_write tmp/temp.sam [deindent {
		@HD	VN:1.4	GO:none	SO:coordinate
		@SQ	SN:chr2	LN:243199373
		@SQ	SN:chr1	LN:249250621
		@RG	ID:sample1	PL:illumina	PU:sample1	LB:solexa-123	SM:sample1
		A3	99	chr2	50	60	20M	=	60	30	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A3	147	chr2	60	60	20M	=	50	-30	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A4	99	chr2	100	60	50M	=	100	50	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	--------------------------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A4	147	chr2	100	60	50M	=	100	-50	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	--------------------------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A1	99	chr1	100	60	20M	=	121	41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A2	99	chr1	102	60	4S2M2I12M	=	111	39	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A2	147	chr1	111	60	30M	=	102	-39	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A1	147	chr1	121	60	20M	=	100	-41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
	}]\n
	file_write tmp/expected.sam [deindent {
		@HD	VN:1.4	GO:none	SO:coordinate
		@SQ	SN:chr1	LN:249250621
		@SQ	SN:chr2	LN:243199373
		@RG	ID:sample1	PL:illumina	PU:sample1	LB:solexa-123	SM:sample1
		A1	99	chr1	100	60	20M	=	121	41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A2	99	chr1	102	60	4S2M2I12M	=	111	39	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A2	147	chr1	111	60	30M	=	102	-39	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A1	147	chr1	121	60	20M	=	100	-41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A3	99	chr2	50	60	20M	=	60	30	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A3	147	chr2	60	60	20M	=	50	-30	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A4	99	chr2	100	60	50M	=	100	50	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	--------------------------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
		A4	147	chr2	100	60	50M	=	100	-50	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	--------------------------------------------------	RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25
	}]\n
	exec samtools view --no-PG -b -1 tmp/temp.sam > tmp/temp.bam
	exec cg bam2sam tmp/temp.bam tmp/result.sam >@ stdout 2>@ stderr
	exec diff tmp/result.sam tmp/expected.sam
} {}

test cg_bam2sam {basic -sortchr 0} {
	file_write tmp/temp.sam [deindent {
		@HD	VN:1.4	GO:none	SO:coordinate
		@SQ	SN:chr2	LN:243199373
		@SQ	SN:chr1	LN:249250621
		@RG	ID:sample1	PL:illumina	PU:sample1	LB:solexa-123	SM:sample1
		A3	99	chr2	50	60	20M	=	60	30	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	NM:i:4	MQ:i:60	AS:i:241	XS:i:25
		A3	147	chr2	60	60	20M	=	50	-30	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	NM:i:4	MQ:i:60	AS:i:241	XS:i:25
		A4	99	chr2	100	60	50M	=	100	50	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	--------------------------------------------------	RG:Z:sample1	NM:i:4	MQ:i:60	AS:i:241	XS:i:25
		A4	147	chr2	100	60	50M	=	100	-50	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	--------------------------------------------------	RG:Z:sample1	NM:i:4	MQ:i:60	AS:i:241	XS:i:25
		A1	99	chr1	100	60	20M	=	121	41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	NM:i:4	MQ:i:60	AS:i:241	XS:i:25
		A2	99	chr1	102	60	4S2M2I12M	=	111	39	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	NM:i:4	MQ:i:60	AS:i:241	XS:i:25
		A2	147	chr1	111	60	30M	=	102	-39	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	------------------------------	RG:Z:sample1	NM:i:4	MQ:i:60	AS:i:241	XS:i:25
		A1	147	chr1	121	60	20M	=	100	-41	AAAAAAAAAAAAAAAAAAAA	--------------------	RG:Z:sample1	NM:i:4	MQ:i:60	AS:i:241	XS:i:25
	}]\n
	exec samtools view --no-PG -b -1 tmp/temp.sam > tmp/temp.bam
	exec cg bam2sam -sortchr 0 tmp/temp.bam > tmp/result.sam
	exec diff tmp/result.sam tmp/temp.sam
} {}

test bam_sort {bam_sort basic} {
	test_cleantmp
	exec samtools sort --no-PG -n -o tmp/test.bam data/bwa.bam
	cg bam_sort tmp/test.bam tmp/sorted.bam
	exec diff tmp/sorted.bam data/bwa.bam
} {}

test bam_sort {bam_sort pipe} {
	test_cleantmp
	exec samtools sort --no-PG -n -o tmp/test.bam data/bwa.bam
	cg bam_sort < tmp/test.bam > tmp/sorted.bam
	exec diff tmp/sorted.bam data/bwa.bam
} {}

test bam_sort {bam_sort pipe sam} {
	test_cleantmp
	exec samtools view --no-PG -h data/bwa.bam > tmp/expected.sam
	exec samtools sort --no-PG -n data/bwa.bam | samtools view --no-PG -h > tmp/test.sam
	cg bam_sort -inputformat sam -outputformat sam < tmp/test.sam > tmp/sorted.sam
	cg tsvdiff tmp/sorted.sam tmp/expected.sam
} {}

test bam_sort {bam_sort to cram} {
	test_cleantmp
	exec samtools sort --no-PG -n -o tmp/test.bam data/bwa.bam
	cg bam_sort -refseq $::refseqdir/hg19 tmp/test.bam tmp/sorted.cram
	cramdiff tmp/sorted.cram data/bwa.sam
} {}

test bam_sort {bam_sort -method gnusort} {
	test_cleantmp
	exec samtools sort --no-PG -n -o tmp/test.bam data/bwa.bam
	cg bam_sort -method gnusort tmp/test.bam tmp/sorted.bam
	cg tsvdiff tmp/sorted.bam data/bwa.bam
} {}

test bam_sort {bam_sort -sort name} {
	test_cleantmp
	exec samtools sort --no-PG -n -o tmp/expected.bam data/bwa.bam
	cg bam_sort -sort name data/bwa.bam tmp/sorted.bam
	exec diff tmp/sorted.bam tmp/expected.bam
} {}

# biobambam no longer in default distro
#test bam_sort {bam_sort -method biobambam} {
#	test_cleantmp
#	exec samtools sort --no-PG -n -o tmp/test.bam data/bwa.bam
#	cg bam_sort -method biobambam tmp/test.bam tmp/sorted.bam
#	cg tsvdiff tmp/sorted.bam data/bwa.bam
#} {}
#
#test bam_sort {bam_sort -method biobambam pipe} {
#	test_cleantmp
#	exec samtools sort --no-PG -n data/bwa.bam > tmp/test.bam
#	cg bam_sort -method biobambam < tmp/test.bam > tmp/sorted.bam
#	cg tsvdiff tmp/sorted.bam data/bwa.bam
#} {}
#
#test bam_sort {bam_sort -method biobambam sam} {
#	test_cleantmp
#	exec samtools view --no-PG -h data/bwa.bam > tmp/expected.sam
#	exec samtools sort --no-PG -n data/bwa.bam | samtools view --no-PG -h > tmp/test.sam
#	cg bam_sort -inputformat sam -outputformat sam -method biobambam tmp/test.sam tmp/sorted.sam
#	cg tsvdiff tmp/sorted.sam tmp/expected.sam
#} {}
#
#test bam_sort {bam_sort -method biobambam -sort name} {
#	test_cleantmp
#	exec samtools sort --no-PG -n -o tmp/expected.bam data/bwa.bam
#	cg bam_sort -method biobambam -sort name data/bwa.bam tmp/sorted.bam
#	cg tsvdiff tmp/sorted.bam tmp/expected.bam
#} {}

test bam_sort {(internal) cg _sam_sort_gnusort, without @HD line} {
	test_cleantmp
	exec samtools sort --no-PG -o tmp/expected.sam data/bwa.bam
	exec samtools sort --no-PG -n -O sam data/bwa.bam | grep -v @HD > tmp/namesorted.sam
	exec cg _sam_sort_gnusort coordinate < tmp/namesorted.sam > tmp/sorted.sam
	exec diff tmp/sorted.sam tmp/expected.sam
} {1c1
< @HD	VN:1.6	SO:coordinate
---
> @HD	VN:1.3	SO:coordinate
child process exited abnormally} error

test noise {pretest} {
	file_write tmp/temp.pup [deindent {
		chr21	9439154	A	3	.gG	CDI
		chr21	9439155	C	3	.,.	CDI
		chr21	9439389	G	40	.....,.-1C..,,.......,,,....,.....,.,,..,..	@B9@DI>BDIGD<DDF=;DDDDFHJDIGGJIDGD?HDDF@
		chr21	9439515	C	53	,$....T*+3gct.,,.,+3gct.T.,T.t..,T.,+3gct.+3GCT.+3GCT..,+3gct.+3GCTT,t,,+3gct.t.+3GCT.tT,,,t,T,+3gct,+3gctt,.	CC>ADD<DIID<CDCFDC;DCJDD<<<7J<<IDDD/JC<GDJ>DBCDJ<<CD@
		chr21	9439626	C	20	,,,,.,,,,,,,.,..,.,,	FFFFDHGHGFFD;<IFDH?8
		chr21	9439632	C	20	,,,.,,,,,,,,..,.,,,^],	FFFDFHJIJAJDJHDHD0@>
		chr21	9439651	T	22	,,g.,..,.,,,g,G,.,,..,	@DIBIEDDFDB<DD;DD?DHDB
		chr21	9439680	A	26	.$,..,.,,,,,.,.,,..,.,,,,..	0FBDJ>J<DBDCDHDDHIDJD0D?GF
		chr21	9439702	C	26	,.,,,g,G,A,,..,.g,,,..,.g,	FCJCIHI>IEFJHIFHD<<5JIBF<>
		chr21	9439709	C	27	.,,,a,A,.,,..,.a,,,..,.aa,,	<HDHJI9JCFIDFJIGII8JGDGCED>
		chr21	9444119	c	57	.$.,..........................,....,.......,.............^!.	CEDEEDE@EHH=IJIJIHJHJIIJJIJJJCJJJ<IHJIIFIAFIEJGHHHHHDC:@?
		chr21	43539135	G	26	...A........,.,..,....,,,^],	DDFEHGJEJJJI@JDGIDJJJJDDCD
		chr21	43539136	C	26	..-2NN...T......,.,..,....,,,,	DDD@H<JHGJJBDJDAJCIJIJDDAD
		chr21	43539137	A	27	..*..-1N....Aa.,.,..,....,,,,^].	DEEDHDJ>HJJBCJCHJCJIJJDCCDC
		chr21	43539138	G	27	..*..*.......,.,..,....,,,,.	DDDFHHIFJJJGCJ@BJCJJIJDCACC
		chr21	43539139	A	27	............,.,..,....,,,,.	?DCFFCIGIJJG>ICGJDJGJJDC>AC
	}]\n
	file_write tmp/expected.tsv [deindent {
		depth	nrdiff	bin	chromosome	pos	ref	depth	bases	qual
		40	0	0.00	chr21	9439389	G	40	.....,.-1C..,,.......,,,....,.....,.,,..,..	@B9@DI>BDIGD<DDF=;DDDDFHJDIGGJIDGD?HDDF@
		53	25	50.00	chr21	9439515	C	53	,$....T*+3gct.,,.,+3gct.T.,T.t..,T.,+3gct.+3GCT.+3GCT..,+3gct.+3GCTT,t,,+3gct.t.+3GCT.tT,,,t,T,+3gct,+3gctt,.	CC>ADD<DIID<CDCFDC;DCJDD<<<7J<<IDDD/JC<GDJ>DBCDJ<<CD@
		20	0	0.00	chr21	9439626	C	20	,,,,.,,,,,,,.,..,.,,	FFFFDHGHGFFD;<IFDH?8
		20	0	0.00	chr21	9439632	C	20	,,,.,,,,,,,,..,.,,,^],	FFFDFHJIJAJDJHDHD0@>
		22	3	15.00	chr21	9439651	T	22	,,g.,..,.,,,g,G,.,,..,	@DIBIEDDFDB<DD;DD?DHDB
		26	0	0.00	chr21	9439680	A	26	.$,..,.,,,,,.,.,,..,.,,,,..	0FBDJ>J<DBDCDHDDHIDJD0D?GF
		26	5	20.00	chr21	9439702	C	26	,.,,,g,G,A,,..,.g,,,..,.g,	FCJCIHI>IEFJHIFHD<<5JIBF<>
		27	5	20.00	chr21	9439709	C	27	.,,,a,A,.,,..,.a,,,..,.aa,,	<HDHJI9JCFIDFJIGII8JGDGCED>
		57	0	0.00	chr21	9444119	c	57	.$.,..........................,....,.......,.............^!.	CEDEEDE@EHH=IJIJIHJHJIIJJIJJJCJJJ<IHJIIFIAFIEJGHHHHHDC:@?
		26	1	5.00	chr21	43539135	G	26	...A........,.,..,....,,,^],	DDFEHGJEJJJI@JDGIDJJJJDDCD
		26	1	5.00	chr21	43539136	C	26	..-2NN...T......,.,..,....,,,,	DDD@H<JHGJJBDJDAJCIJIJDDAD
		27	3	15.00	chr21	43539137	A	27	..*..-1N....Aa.,.,..,....,,,,^].	DEEDHDJ>HJJBCJCHJCJIJJDCCDC
		27	2	10.00	chr21	43539138	G	27	..*..*.......,.,..,....,,,,.	DDDFHHIFJJJGCJ@BJCJJIJDCACC
		27	0	0.00	chr21	43539139	A	27	............,.,..,....,,,,.	?DCFFCIGIJJG>ICGJDJGJJDC>AC
	}]\n
	exec noise 20 0 0 < tmp/temp.pup > tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
	# cg select -overwrite 1 -f {chromosome {begin=$pos - 1} depth nrdiff pct_alt="${bin}"} tmp/expected.tsv tmp/expected2.tsv
	file_write tmp/expected2.tsv [deindent {
		chromosome	begin	depth	ref	nr_diff	nr_mismatch	nr_del	nr_ins	nr_A	nr_C	nr_G	nr_T	pct_diff	pct_mismatch	pct_del	pct_ins	pct_A	pct_C	pct_G	pct_T	pct_maxmismatch
		chr21	9439388	40	G	0	0	0	0	0	40	0	0	0.00	0.00	0.00	0.00	0.00	100.00	0.00	0.00	0.00
		chr21	9439514	53	C	25	13	1	11	0	39	0	13	47.17	24.53	1.89	20.75	0.00	73.58	0.00	24.53	24.53
		chr21	9439625	20	C	0	0	0	0	0	20	0	0	0.00	0.00	0.00	0.00	0.00	100.00	0.00	0.00	0.00
		chr21	9439631	20	C	0	0	0	0	0	20	0	0	0.00	0.00	0.00	0.00	0.00	100.00	0.00	0.00	0.00
		chr21	9439650	22	T	3	3	0	0	0	0	3	19	13.64	13.64	0.00	0.00	0.00	0.00	13.64	86.36	13.64
		chr21	9439679	26	A	0	0	0	0	26	0	0	0	0.00	0.00	0.00	0.00	100.00	0.00	0.00	0.00	0.00
		chr21	9439701	26	C	5	5	0	0	1	21	4	0	19.23	19.23	0.00	0.00	3.85	80.77	15.38	0.00	15.38
		chr21	9439708	27	C	5	5	0	0	5	22	0	0	18.52	18.52	0.00	0.00	18.52	81.48	0.00	0.00	18.52
		chr21	9444118	57	c	0	0	0	0	0	57	0	0	0.00	0.00	0.00	0.00	0.00	100.00	0.00	0.00	0.00
		chr21	43539134	26	G	1	1	0	0	1	25	0	0	3.85	3.85	0.00	0.00	3.85	96.15	0.00	0.00	3.85
		chr21	43539135	26	C	1	1	0	0	0	25	0	1	3.85	3.85	0.00	0.00	0.00	96.15	0.00	3.85	3.85
		chr21	43539136	27	A	3	2	1	0	24	0	0	0	11.11	7.41	3.70	0.00	88.89	0.00	0.00	0.00	7.41
		chr21	43539137	27	G	2	0	2	0	0	26	0	0	7.41	0.00	7.41	0.00	0.00	96.30	0.00	0.00	0.00
		chr21	43539138	27	A	0	0	0	0	27	0	0	0	0.00	0.00	0.00	0.00	100.00	0.00	0.00	0.00	0.00
	}]\n
	exec noise 20 1 < tmp/temp.pup > tmp/result2.tsv
	exec diff tmp/result2.tsv tmp/expected2.tsv
	file_write tmp/expected3.tsv [deindent {
		#filetype	tsv/noisefile
		#fileversion	0.99
		#description	a noise (=percentage alternative allele) histogram for all positions with depth >= mindepth 
		#mindepth	20
		#general bam stats:
		#numalignedbases	427
		#nummismatch	32
		#numdel	4
		#numins	11
		percentage	count
		0.00	6
		5.00	2
		10.00	1
		15.00	2
		20.00	2
		25.00	0
		30.00	0
		35.00	0
		40.00	0
		45.00	0
		50.00	1
		55.00	0
		60.00	0
		65.00	0
		70.00	0
		75.00	0
		80.00	0
		85.00	0
		90.00	0
		95.00	0
		100.00	0
	}]\n
	# cg select -f {{pct=100.0*$nrdiff/$depth} depth nrdiff} tmp/result.tsv
	# cg select -f {{pct=100.0*$nrdiff/$depth}} -g pct tmp/result.tsv
	exec noise 20 < tmp/temp.pup > tmp/result3.tsv
	exec diff tmp/result3.tsv tmp/expected3.tsv
	file_write tmp/expected4.tsv [deindent {
		#filetype	tsv/noisefile
		#fileversion	0.99
		#description	a noise (=percentage alternative allele) histogram for all positions with depth >= mindepth 
		#mindepth	20
		#general bam stats:
		#numalignedbases	427
		#nummismatch	32
		#numdel	4
		#numins	11
		percentage	count
		0.00	7
		5.00	2
		10.00	1
		15.00	1
		20.00	2
		25.00	1
		30.00	0
		35.00	0
		40.00	0
		45.00	0
		50.00	0
		55.00	0
		60.00	0
		65.00	0
		70.00	0
		75.00	0
		80.00	0
		85.00	0
		90.00	0
		95.00	0
		100.00	0
	}]\n
	# cg select -f {{pct=100.0*$nrdiff/$depth} depth nrdiff} tmp/result.tsv
	# cg select -f {{pct=100.0*$nrdiff/$depth}} -g pct tmp/result.tsv
	exec noise 20 0 -1 0 b < tmp/temp.pup > tmp/result4.tsv
	exec diff tmp/result4.tsv tmp/expected4.tsv
} {}

test noise {noise -typediffs b} {
	test_cleantmp
	exec samtools view -b data/map-rdsbwa-NA19240_small.sam > tmp/temp.bam
	exec samtools index tmp/temp.bam
	set bamfile tmp/temp.bam
	set refseq $::refseqdir/hg19
	# exec samtools mpileup -f $refseq/genome_hg19.ifas -A --ignore-overlaps tmp/temp.bam > tmp/temp.pup
	file_write tmp/expected.tsv [deindent {
		#filetype	tsv/noisefile
		#fileversion	0.99
		#description	a noise (=percentage alternative allele) histogram for all positions with depth >= mindepth 
		#mindepth	20
		#general bam stats:
		#numalignedbases	8709
		#nummismatch	103
		#numdel	73
		#numins	3
		percentage	count
		0.00	128
		5.00	10
		10.00	1
		15.00	0
		20.00	0
		25.00	0
		30.00	0
		35.00	0
		40.00	0
		45.00	0
		50.00	0
		55.00	0
		60.00	0
		65.00	0
		70.00	0
		75.00	0
		80.00	0
		85.00	0
		90.00	0
		95.00	0
		100.00	1
	}]\n
	cg noise -typediffs base -refseq $refseq $bamfile > tmp/result.tsv
	catch {
		exec diff tmp/result.tsv tmp/expected.tsv
	}
} 0

test noise {noise -perpos 1} {
	test_cleantmp
	exec samtools view -b data/map-rdsbwa-NA19240_small.sam > tmp/test.bam
	exec samtools index tmp/test.bam
	set refseq $::refseqdir/hg19
	file_write tmp/expected.tsv [deindent {
		chromosome	begin	depth	ref	nr_diff	nr_mismatch	nr_del	nr_ins	nr_A	nr_C	nr_G	nr_T	pct_diff	pct_mismatch	pct_del	pct_ins	pct_A	pct_C	pct_G	pct_T	pct_maxmismatch
		chr21	9825515	80	G	0	0	0	0	0	80	0	0	0.00	0.00	0.00	0.00	0.00	100.00	0.00	0.00	0.00
		chr21	9825516	80	G	0	0	0	0	0	80	0	0	0.00	0.00	0.00	0.00	0.00	100.00	0.00	0.00	0.00
		chr21	9825517	81	T	0	0	0	0	0	0	0	81	0.00	0.00	0.00	0.00	0.00	0.00	0.00	100.00	0.00
		chr21	9825518	81	G	0	0	0	0	0	81	0	0	0.00	0.00	0.00	0.00	0.00	100.00	0.00	0.00	0.00
		chr21	9825519	82	C	0	0	0	0	0	82	0	0	0.00	0.00	0.00	0.00	0.00	100.00	0.00	0.00	0.00
		chr21	9825520	82	C	0	0	0	0	0	82	0	0	0.00	0.00	0.00	0.00	0.00	100.00	0.00	0.00	0.00
		chr21	9825521	83	C	0	0	0	0	0	83	0	0	0.00	0.00	0.00	0.00	0.00	100.00	0.00	0.00	0.00
		chr21	9825522	85	T	0	0	0	0	0	0	0	85	0.00	0.00	0.00	0.00	0.00	0.00	0.00	100.00	0.00
		chr21	9825523	85	T	8	8	0	0	0	8	0	77	9.41	9.41	0.00	0.00	0.00	9.41	0.00	90.59	9.41
		chr21	9825524	85	G	3	0	3	0	0	82	0	0	3.53	0.00	3.53	0.00	0.00	96.47	0.00	0.00	0.00
		chr21	9825525	85	C	3	0	3	0	0	82	0	0	3.53	0.00	3.53	0.00	0.00	96.47	0.00	0.00	0.00
		chr21	9825526	85	C	3	0	3	0	0	82	0	0	3.53	0.00	3.53	0.00	0.00	96.47	0.00	0.00	0.00
		chr21	9825527	85	C	3	0	3	0	0	82	0	0	3.53	0.00	3.53	0.00	0.00	96.47	0.00	0.00	0.00
		chr21	9825528	85	T	3	0	3	0	0	0	0	82	3.53	0.00	3.53	0.00	0.00	0.00	0.00	96.47	0.00
		chr21	9825529	86	C	3	0	3	0	0	83	0	0	3.49	0.00	3.49	0.00	0.00	96.51	0.00	0.00	0.00
		chr21	9825530	86	A	86	83	3	0	0	0	83	0	100.00	96.51	3.49	0.00	0.00	0.00	96.51	0.00	96.51
		chr21	9825531	86	C	3	0	3	0	0	83	0	0	3.49	0.00	3.49	0.00	0.00	96.51	0.00	0.00	0.00
		chr21	9825532	86	G	3	0	3	0	0	83	0	0	3.49	0.00	3.49	0.00	0.00	96.51	0.00	0.00	0.00
		chr21	9825533	87	G	3	0	3	0	0	84	0	0	3.45	0.00	3.45	0.00	0.00	96.55	0.00	0.00	0.00
		chr21	9825534	87	T	3	0	3	0	0	0	0	84	3.45	0.00	3.45	0.00	0.00	0.00	0.00	96.55	0.00
		chr21	9825535	87	C	3	0	3	0	0	84	0	0	3.45	0.00	3.45	0.00	0.00	96.55	0.00	0.00	0.00
		chr21	9825536	87	C	3	0	3	0	0	84	0	0	3.45	0.00	3.45	0.00	0.00	96.55	0.00	0.00	0.00
		chr21	9825537	87	C	3	0	3	0	0	84	0	0	3.45	0.00	3.45	0.00	0.00	96.55	0.00	0.00	0.00
		chr21	9825538	87	C	3	0	3	0	0	84	0	0	3.45	0.00	3.45	0.00	0.00	96.55	0.00	0.00	0.00
		chr21	9825539	87	G	4	1	3	0	1	83	0	0	4.60	1.15	3.45	0.00	1.15	95.40	0.00	0.00	1.15
		chr21	9825540	85	G	3	0	3	0	0	82	0	0	3.53	0.00	3.53	0.00	0.00	96.47	0.00	0.00	0.00
		chr21	9825541	83	C	3	0	3	0	0	80	0	0	3.61	0.00	3.61	0.00	0.00	96.39	0.00	0.00	0.00
		chr21	9825542	83	C	3	0	3	0	0	80	0	0	3.61	0.00	3.61	0.00	0.00	96.39	0.00	0.00	0.00
		chr21	9825543	82	C	3	0	3	0	0	79	0	0	3.66	0.00	3.66	0.00	0.00	96.34	0.00	0.00	0.00
		chr21	9825544	81	T	6	0	3	3	0	0	0	78	7.41	0.00	3.70	3.70	0.00	0.00	0.00	96.30	0.00
	}]\n
	# exec samtools mpileup -f $refseq/genome_hg19.ifas -A --ignore-overlaps tmp/temp.bam > tmp/temp.pup
	cg noise -mindepth 80 -refseq $refseq -perpos 1 data/map-rdsbwa-NA19240_small.sam > tmp/result.tsv
	catch {
		exec diff tmp/result.tsv tmp/expected.tsv
	}
} 0

test noise {noise -perpos 1} {
	test_cleantmp
	exec samtools view -b data/map-rdsbwa-NA19240_small.sam > tmp/test.bam
	exec samtools index tmp/test.bam
	set refseq $::refseqdir/hg19
	file_write tmp/expected.tsv [deindent {
		chromosome	begin	depth	ref	nr_diff	nr_mismatch	nr_del	nr_ins	nr_A	nr_C	nr_G	nr_T	pct_diff	pct_mismatch	pct_del	pct_ins	pct_A	pct_C	pct_G	pct_T	pct_maxmismatch
		chr21	9825515	80	G	0	0	0	0	0	80	0	0	0.00	0.00	0.00	0.00	0.00	100.00	0.00	0.00	0.00
		chr21	9825516	80	G	0	0	0	0	0	80	0	0	0.00	0.00	0.00	0.00	0.00	100.00	0.00	0.00	0.00
		chr21	9825517	81	T	0	0	0	0	0	0	0	81	0.00	0.00	0.00	0.00	0.00	0.00	0.00	100.00	0.00
		chr21	9825518	81	G	0	0	0	0	0	81	0	0	0.00	0.00	0.00	0.00	0.00	100.00	0.00	0.00	0.00
		chr21	9825519	82	C	0	0	0	0	0	82	0	0	0.00	0.00	0.00	0.00	0.00	100.00	0.00	0.00	0.00
		chr21	9825520	82	C	0	0	0	0	0	82	0	0	0.00	0.00	0.00	0.00	0.00	100.00	0.00	0.00	0.00
		chr21	9825521	83	C	0	0	0	0	0	83	0	0	0.00	0.00	0.00	0.00	0.00	100.00	0.00	0.00	0.00
		chr21	9825522	85	T	0	0	0	0	0	0	0	85	0.00	0.00	0.00	0.00	0.00	0.00	0.00	100.00	0.00
		chr21	9825523	85	T	8	8	0	0	0	8	0	77	9.41	9.41	0.00	0.00	0.00	9.41	0.00	90.59	9.41
		chr21	9825524	85	G	3	0	3	0	0	82	0	0	3.53	0.00	3.53	0.00	0.00	96.47	0.00	0.00	0.00
		chr21	9825525	85	C	3	0	3	0	0	82	0	0	3.53	0.00	3.53	0.00	0.00	96.47	0.00	0.00	0.00
		chr21	9825526	85	C	3	0	3	0	0	82	0	0	3.53	0.00	3.53	0.00	0.00	96.47	0.00	0.00	0.00
		chr21	9825527	85	C	3	0	3	0	0	82	0	0	3.53	0.00	3.53	0.00	0.00	96.47	0.00	0.00	0.00
		chr21	9825528	85	T	3	0	3	0	0	0	0	82	3.53	0.00	3.53	0.00	0.00	0.00	0.00	96.47	0.00
		chr21	9825529	86	C	3	0	3	0	0	83	0	0	3.49	0.00	3.49	0.00	0.00	96.51	0.00	0.00	0.00
		chr21	9825530	86	A	86	83	3	0	0	0	83	0	100.00	96.51	3.49	0.00	0.00	0.00	96.51	0.00	96.51
		chr21	9825531	86	C	3	0	3	0	0	83	0	0	3.49	0.00	3.49	0.00	0.00	96.51	0.00	0.00	0.00
		chr21	9825532	86	G	3	0	3	0	0	83	0	0	3.49	0.00	3.49	0.00	0.00	96.51	0.00	0.00	0.00
		chr21	9825533	87	G	3	0	3	0	0	84	0	0	3.45	0.00	3.45	0.00	0.00	96.55	0.00	0.00	0.00
		chr21	9825534	87	T	3	0	3	0	0	0	0	84	3.45	0.00	3.45	0.00	0.00	0.00	0.00	96.55	0.00
		chr21	9825535	87	C	3	0	3	0	0	84	0	0	3.45	0.00	3.45	0.00	0.00	96.55	0.00	0.00	0.00
		chr21	9825536	87	C	3	0	3	0	0	84	0	0	3.45	0.00	3.45	0.00	0.00	96.55	0.00	0.00	0.00
		chr21	9825537	87	C	3	0	3	0	0	84	0	0	3.45	0.00	3.45	0.00	0.00	96.55	0.00	0.00	0.00
		chr21	9825538	87	C	3	0	3	0	0	84	0	0	3.45	0.00	3.45	0.00	0.00	96.55	0.00	0.00	0.00
		chr21	9825539	87	G	4	1	3	0	1	83	0	0	4.60	1.15	3.45	0.00	1.15	95.40	0.00	0.00	1.15
		chr21	9825540	85	G	3	0	3	0	0	82	0	0	3.53	0.00	3.53	0.00	0.00	96.47	0.00	0.00	0.00
		chr21	9825541	83	C	3	0	3	0	0	80	0	0	3.61	0.00	3.61	0.00	0.00	96.39	0.00	0.00	0.00
		chr21	9825542	83	C	3	0	3	0	0	80	0	0	3.61	0.00	3.61	0.00	0.00	96.39	0.00	0.00	0.00
		chr21	9825543	82	C	3	0	3	0	0	79	0	0	3.66	0.00	3.66	0.00	0.00	96.34	0.00	0.00	0.00
		chr21	9825544	81	T	6	0	3	3	0	0	0	78	7.41	0.00	3.70	3.70	0.00	0.00	0.00	96.30	0.00
	}]\n
	# exec samtools mpileup -f $refseq/genome_hg19.ifas -A --ignore-overlaps tmp/temp.bam > tmp/temp.pup
	cg noise -mindepth 80 -refseq $refseq -perpos 1 data/map-rdsbwa-NA19240_small.sam > tmp/result.tsv
	catch {
		exec diff tmp/result.tsv tmp/expected.tsv
	}
} 0

test bam_varia {longshot_replacebam} {
	test_cleantmp
	file copy -force data/minimap2.bam tmp/minimap2.bam
	file mtime tmp/minimap2.bam 1607468400
	file copy -force data/minimap2-p.bam tmp/minimap2-p.bam
	file mtime tmp/minimap2-p.bam 1607554800
	longshot_replacebam tmp/minimap2-p.bam tmp/minimap2.bam
	file lstat tmp/minimap2.bam a
	list [file link tmp/minimap2.bam] $a(mtime) [file size tmp/minimap2.bam] [file mtime tmp/minimap2-p.bam] [file size tmp/minimap2-p.bam]
} {minimap2-p.bam 1607468400 19446 1607554800 19446}

testsummarize
