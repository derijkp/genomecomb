#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl
package require genomecomb

test svmulticompar {basic} {
	test_cleantmp
	cg svmulticompar tmp/temp.tsv data/cgsv1.tsv data/cgsv2.tsv data/cgsv3.tsv
	exec diff tmp/temp.tsv data/expected-svmulticompar.tsv
} {} 

test svmulticompar {name} {
	test_cleantmp
	file copy data/cgsv1.tsv tmp/sv-cgsv-s1.tsv
	file copy data/cgsv2.tsv tmp/sv-cgsv-s2.tsv
	file copy data/cgsv3.tsv tmp/sv-test-s2.tsv
	cg svmulticompar tmp/temp.tsv tmp/sv-cgsv-s1.tsv tmp/sv-cgsv-s2.tsv tmp/sv-test-s2.tsv
	cg select -sh /dev/null tmp/temp.tsv tmp/shtemp.tsv
	cg select -sh /dev/null data/expected-svmulticompar.tsv tmp/expected.tsv
	exec diff tmp/shtemp.tsv tmp/expected.tsv
	cg select -a tmp/temp.tsv
} {cgsv-s1
cgsv-s2
test-s2} 

test svmulticompar {trans} {
	test_cleantmp
	write_tab tmp/s1.tsv {
		chromosome	begin	end	type	ref	alt
		chr1	10001	10001	trans	{}	[chr2:100[
	}
	write_tab tmp/s2.tsv {
		chromosome	begin	end	type	ref	alt
		chr1	10001	10001	trans	{}	[chr2:500[
		chr1	10021	10021	trans	{}	[chr2:200[
	}
	file_write tmp/expected.tsv [deindent {
		chromosome	begin	end	type	ref	alt	lbegin-s1	lend-s1	lalt-s1	lbegin-s2	lend-s2	lalt-s2
		1	10001	10001	trans		[2:100[	10001	10001	[2:100[	10021	10021	[2:200[
		1	10001	10001	trans		[2:500[				10001	10001	[2:500[
	}]\n
	cg svmulticompar tmp/temp.tsv tmp/s1.tsv tmp/s2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test svmulticompar {bnd} {
	test_cleantmp
	write_tab tmp/s1.tsv {
		chromosome	begin	end	type	ref	alt
		chr1	10001	10001	bnd	{}	[chr2:100[
	}
	write_tab tmp/s2.tsv {
		chromosome	begin	end	type	ref	alt
		chr1	10002	10002	bnd	{}	]chr2:100]
		chr1	10021	10021	bnd	{}	[chr2:200[
	}
	file_write tmp/expected.tsv [deindent {
		chromosome	begin	end	type	ref	alt	lbegin-s1	lend-s1	lalt-s1	lbegin-s2	lend-s2	lalt-s2
		1	10001	10001	bnd		[2:100[	10001	10001	[2:100[	10021	10021	[2:200[
		1	10002	10002	bnd		]2:100]				10002	10002	]2:100]
	}]\n
	cg svmulticompar tmp/temp.tsv tmp/s1.tsv tmp/s2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test svmulticompar {ins} {
	test_cleantmp
	write_tab tmp/s1.tsv {
		chromosome	begin	end	type	ref	alt
		chr1	10001	10001	ins	{}	AGCTAGCT
	}
	write_tab tmp/s2.tsv {
		chromosome	begin	end	type	ref	alt
		chr1	10001	10001	ins	{}	TGCTAGCTA
		chr1	10021	10021	ins	{}	20
		chr1	20021	20021	ins	{}	20
	}
	file_write tmp/expected.tsv [deindent {
		chromosome	begin	end	type	ref	alt	lbegin-s1	lend-s1	lalt-s1	lbegin-s2	lend-s2	lalt-s2
		1	10001	10001	ins		AGCTAGCT	10001	10001	AGCTAGCT	10001	10001	TGCTAGCTA
		1	10021	10021	ins		20				10021	10021	20
		1	20021	20021	ins		20				20021	20021	20
	}]\n
	cg svmulticompar tmp/temp.tsv tmp/s1.tsv tmp/s2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test svmulticompar {cnv} {
	test_cleantmp
	write_tab tmp/s1.tsv {
		chr	begin	end	type	avgNormalizedCvg	relativeCvg	calledPloidy	ploidyScore
		chr1	1720	2721	amp	84.8	1.95	4	9
		chr1	9650	10000	del	76.5	1.76	4	9
	}
	write_tab tmp/s2.tsv {
		chr	begin	end	type	avgNormalizedCvg	relativeCvg	calledPloidy	ploidyScore
		chr1	1700	2721	amp	84.8	1.95	3	9
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	lbegin-s1	lend-s1	lalt-s1	avgNormalizedCvg-s1	relativeCvg-s1	calledPloidy-s1	ploidyScore-s1	lbegin-s2	lend-s2	lalt-s2	avgNormalizedCvg-s2	relativeCvg-s2	calledPloidy-s2	ploidyScore-s2
		1	1720	2721	amp	1001	?	1720	2721	?	84.8	1.95	4	9	1700	2721	?	84.8	1.95	3	9
		1	9650	10000	del	350	{}	9650	10000	{}	76.5	1.76	4	9	{}	{}	{}	{}	{}	{}	{}
	}
	cg svmulticompar tmp/temp.tsv tmp/s1.tsv tmp/s2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test svmulticompar {sv del duplicated error} {
	test_cleantmp
	write_tab tmp/s1.tsv {
		chromosome	begin	end	type
		chr1	5727767	5728301	del
	}
	write_tab tmp/s2.tsv {
		chromosome	begin	end	type
		chr1	5727878	5728315	del
	}
	write_tab tmp/s3.tsv {
		chromosome	begin	end	type
		chr1	5727847	5728315	del
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	lbegin-s1	lend-s1	lalt-s1	lbegin-s2	lend-s2	lalt-s2	lbegin-s3	lend-s3	lalt-s3
		1	5727767	5728301	del	534	{}	5727767	5728301	{}	{}	{}	{}	{}	{}	{}
		1	5727878	5728315	del	437	{}	{}	{}	{}	5727878	5728315	{}	5727847	5728315	{}
	}
	file delete tmp/temp.tsv
	cg svmulticompar -overlap 80 tmp/temp.tsv tmp/s1.tsv tmp/s2.tsv tmp/s3.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test svmulticompar {sv repeated del} {
	test_cleantmp
	write_tab tmp/s1.tsv {
		chromosome	begin	end	type
		chr1	5727767	5728301	del
	}
	write_tab tmp/s2.tsv {
		chromosome	begin	end	type	name
		chr1	5727767	5728301	del	1
		chr1	5727767	5728301	del	2
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	lbegin-s1	lend-s1	lalt-s1	lbegin-s2	lend-s2	lalt-s2	name-s2
		1	5727767	5728301	del	534	{}	{}	{}	{}	5727767	5728301	{}	2
		1	5727767	5728301	del	534	{}	5727767	5728301	{}	5727767	5728301	{}	1
	}
	file delete tmp/temp.tsv
	cg svmulticompar -overlap 80 tmp/temp.tsv tmp/s1.tsv tmp/s2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test svmulticompar {sv trans, diff chr} {
	test_cleantmp
	write_tab tmp/s1.tsv {
		chromosome	begin	end	type	ref	alt
		1	5727760	5727760	trans	{}	[2:1000[
	}
	write_tab tmp/s2.tsv {
		chromosome	begin	end	type	ref	alt
		1	5727767	5728301	trans	534	[chr2:1200[
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	lbegin-s1	lend-s1	lalt-s1	lbegin-s2	lend-s2	lalt-s2
		1	5727760	5727760	trans	{}	[2:1000[	5727760	5727760	[2:1000[	5727767	5728301	[2:1200[
	}
	file delete tmp/temp.tsv
	cg svmulticompar -overlap 80 tmp/temp.tsv tmp/s1.tsv tmp/s2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test svmulticompar {bugcheck wrong length (incorrectly added 2 in output)} {
	test_cleantmp
	file_write tmp/s1.tsv [deindent {
		chromosome	begin	end	type	ref	alt	quality	alleleSeq1	alleleSeq2	zyg	phased	genotypes	filter	genoqual	PL	PR	SR	IMPRECISE	CIPOS	CIEND	CIGAR	MATEID	EVENT	HOMLEN	HOMSEQ	SVINSLEN	SVINSSEQ	LEFT_SVINSSEQ	RIGHT_SVINSSEQ	INV3	INV5	BND_DEPTH	MATE_BND_DEPTH	JUNCTION_QUAL
	}]\n
	file_write tmp/s2.tsv [deindent {
		chromosome	begin	end	type	ref	alt	quality	alleleSeq1	alleleSeq2	zyg	phased	genotypes	filter	genoqual	PL	PR	SR	IMPRECISE	CIPOS	CIEND	CIGAR	MATEID	EVENT	HOMLEN	HOMSEQ	SVINSLEN	SVINSSEQ	LEFT_SVINSSEQ	RIGHT_SVINSSEQ	INV3	INV5	BND_DEPTH	MATE_BND_DEPTH	JUNCTION_QUAL
		chr21	22781178	22781332	inv	154	i	694	154	i	t	0	0;1	PASS	694	744,0,999	42,11	57,16		0,1	-1,0				1	T						1			
	}]\n
	file_write tmp/expected.tsv [deindent {
		chromosome	begin	end	type	ref	alt	lbegin-s1	lend-s1	lalt-s1	quality-s1	alleleSeq1-s1	alleleSeq2-s1	zyg-s1	phased-s1	genotypes-s1	filter-s1	genoqual-s1	PL-s1	PR-s1	SR-s1	IMPRECISE-s1	CIPOS-s1	CIEND-s1	CIGAR-s1	MATEID-s1	EVENT-s1	HOMLEN-s1	HOMSEQ-s1	SVINSLEN-s1	SVINSSEQ-s1	LEFT_SVINSSEQ-s1	RIGHT_SVINSSEQ-s1	INV3-s1	INV5-s1	BND_DEPTH-s1	MATE_BND_DEPTH-s1	JUNCTION_QUAL-s1	lbegin-s2	lend-s2	lalt-s2	quality-s2	alleleSeq1-s2	alleleSeq2-s2	zyg-s2	phased-s2	genotypes-s2	filter-s2	genoqual-s2	PL-s2	PR-s2	SR-s2	IMPRECISE-s2	CIPOS-s2	CIEND-s2	CIGAR-s2	MATEID-s2	EVENT-s2	HOMLEN-s2	HOMSEQ-s2	SVINSLEN-s2	SVINSSEQ-s2	LEFT_SVINSSEQ-s2	RIGHT_SVINSSEQ-s2	INV3-s2	INV5-s2	BND_DEPTH-s2	MATE_BND_DEPTH-s2	JUNCTION_QUAL-s2
		21	22781178	22781332	inv	154	i																																22781178	22781332	i	694	154	i	t	0	0;1	PASS	694	744,0,999	42,11	57,16		0,1	-1,0				1	T						1			
	}]\n
	file delete tmp/result.tsv
	cg svmulticompar tmp/result.tsv tmp/s1.tsv tmp/s2.tsv
	cg checktsv tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {} 

test svmulticompar {bugcheck sort error} {
	test_cleantmp
	file_write tmp/s1.tsv [deindent {
		chromosome	begin	end	type	ref	alt	quality	alleleSeq1	alleleSeq2	zyg	phased	genotypes	filter	genoqual	PL	PR	SR	IMPRECISE	CIPOS	CIEND	CIGAR	MATEID	EVENT	HOMLEN	HOMSEQ	SVINSLEN	SVINSSEQ	LEFT_SVINSSEQ	RIGHT_SVINSSEQ	INV3	INV5	BND_DEPTH	MATE_BND_DEPTH	JUNCTION_QUAL
		chr2	16750000	16750000	bnd		]CHR20:43204700]G	166		]CHR19:43204700]G	t	0	0;1	PASS	166	216,0,955	35,4	39,10					MantaBND:1312:0:6:0:0:0:1				1	G					71	53	
	}]\n
	file_write tmp/s2.tsv [deindent {
		chromosome	begin	end	type	ref	alt	quality	alleleSeq1	alleleSeq2	zyg	phased	genotypes	filter	genoqual	PL	PR	SR	IMPRECISE	CIPOS	CIEND	CIGAR	MATEID	EVENT	HOMLEN	HOMSEQ	SVINSLEN	SVINSSEQ	LEFT_SVINSSEQ	RIGHT_SVINSSEQ	INV3	INV5	BND_DEPTH	MATE_BND_DEPTH	JUNCTION_QUAL
		chr2	16750000	16750000	bnd		[CHR20:57210257[GG	200		[CHR20:57210257[GG	t	0	0;1	PASS	200	250,0,911	49,3	28,9					MantaBND:1209:0:2:0:0:0:1				1	G					35	43	
	}]\n
	file delete tmp/result.tsv
	cg svmulticompar tmp/result.tsv tmp/s1.tsv tmp/s2.tsv
	cg checksort tmp/result.tsv
} {} 

test_cleantmp

set ::env(PATH) $keeppath

testsummarize
