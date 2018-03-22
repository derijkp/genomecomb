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

test svmulticompar {trans} {
	test_cleantmp
	write_tab tmp/s1.tsv {
		chromosome	begin	end	type	ref	alt
		chr1	10001	10001	trans	{}	[chr2:100[
	}
	write_tab tmp/s2.tsv {
		chromosome	begin	end	type	ref	alt
		chr1	10021	10021	trans	{}	[chr2:200[
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	lbegin-s1	lend-s1	lalt-s1	lbegin-s2	lend-s2	lalt-s2
		1	10001	10001	trans	{}	[chr2:100[	10001	10001	[chr2:100[	10021	10021	[chr2:200[
	}
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
		1	5727767	5728301	del	534	{}	5727767	5728301	{}	5727767	5728301	{}	1
		1	5727767	5728301	del	534	{}	{}	{}	{}	5727767	5728301	{}	2
	}
	file delete tmp/temp.tsv
	cg svmulticompar -overlap 80 tmp/temp.tsv tmp/s1.tsv tmp/s2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test_cleantmp

set ::env(PATH) $keeppath

testsummarize
