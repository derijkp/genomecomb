#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test bed2tsv {bed2tsv} {
	exec cg bed2tsv data/sample.bed
} {#browser position chr7:127471196-127495720
#browser hide all
#track name="ColorByStrandDemo" description="Color by strand demonstration" visibility=2 colorByStrand="255,0,0 0,0,255"
chromosome	begin	end	name	score	strand
chr7	127471196  127472363  Pos1  0  +
chr7	127472363  127473530  Pos2  0  +
chr7	127473530  127474697  Pos3  0  +
chr7	127474697  127475864  Pos4  0  +
chr7	127475864  127477031  Neg1  0  -
chr7	127477031  127478198  Neg2  0  -
chr7	127478198  127479365  Neg3  0  -
chr7	127479365  127480532  Pos5  0  +
chr7	127480532  127481699  Neg4  0  -}

test bed2tsv {bed2tsv compress} {
	exec cg bed2tsv data/sample.bed tmp/test.tsv.zst
	cg unzip tmp/test.tsv.zst
	file_read tmp/test.tsv
} {#browser position chr7:127471196-127495720
#browser hide all
#track name="ColorByStrandDemo" description="Color by strand demonstration" visibility=2 colorByStrand="255,0,0 0,0,255"
chromosome	begin	end	name	score	strand
chr7	127471196  127472363  Pos1  0  +
chr7	127472363  127473530  Pos2  0  +
chr7	127473530  127474697  Pos3  0  +
chr7	127474697  127475864  Pos4  0  +
chr7	127475864  127477031  Neg1  0  -
chr7	127477031  127478198  Neg2  0  -
chr7	127478198  127479365  Neg3  0  -
chr7	127479365  127480532  Pos5  0  +
chr7	127480532  127481699  Neg4  0  -
}

test_cleantmp

test splitalleles {splitalleles} {
	write_tab tmp/test.tsv {
		chromosome begin end type ref alt freq-sample
	 	chr1 4000 4001 snp G A 0.5
	 	chr1 4001 4002 snp A G,C 0.5,0.1
	 	chr1 4100 4101 del A {} 0.5
	 	chr1 4200 4200 ins {} A 0.5
	 	chr1 4208 4208 ins {} A,AA 0.8,0.9
	}
	exec cg splitalleles tmp/test.tsv > tmp/temp.tsv
	write_tab tmp/expected.tsv {
		chromosome begin end type ref alt freq-sample
	 	chr1 4000 4001 snp G A 0.5
	 	chr1 4001 4002 snp A C 0.1
	 	chr1 4001 4002 snp A G 0.5
	 	chr1 4100 4101 del A {} 0.5
	 	chr1 4200 4200 ins {} A 0.5
	 	chr1 4208 4208 ins {} A 0.8
	 	chr1 4208 4208 ins {} AA 0.9
	}
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test splitalleles {splitalleles2} {
	write_tab tmp/test.tsv {
		chromosome begin end type ref alt freq-sample
	 	chr1 4001 4002 snp A G,C 0.5,0.1
	 	chr1 4008 4009 snp G C 0.9
	}
	exec cg splitalleles tmp/test.tsv > tmp/temp.tsv
	write_tab tmp/expected.tsv {
		chromosome begin end type ref alt freq-sample
	 	chr1 4001 4002 snp A C 0.1
	 	chr1 4001 4002 snp A G 0.5
	 	chr1 4008 4009 snp G C 0.9
	}
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test splitalleles {splitalleles comment} {
	write_tab tmp/test.tsv {
		chromosome begin end type ref alt freq-sample
	 	chr1 4000 4001 snp G A 0.5
	 	chr1 4001 4002 snp A G,C 0.5,0.1
	} {#comment line}
	exec cg splitalleles tmp/test.tsv > tmp/temp.tsv
	write_tab tmp/expected.tsv {
		chromosome begin end type ref alt freq-sample
	 	chr1 4000 4001 snp G A 0.5
	 	chr1 4001 4002 snp A C 0.1
	 	chr1 4001 4002 snp A G 0.5
	} {#comment line}
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test collapsealleles {collapsealleles} {
	write_tab tmp/test.tsv {
		chromosome begin end type ref alt freq-sample
	 	chr1 4000 4001 snp G A 0.5
	 	chr1 4001 4002 snp A C 0.1
	 	chr1 4001 4002 snp A G 0.5
	 	chr1 4100 4101 del A {} 0.5
	 	chr1 4200 4200 ins {} A 0.5
	 	chr1 4208 4208 ins {} A 0.8
	 	chr1 4208 4208 ins {} AA 0.9
	}
	exec cg collapsealleles tmp/test.tsv > tmp/temp.tsv
	write_tab tmp/expected.tsv {
		chromosome begin end type ref alt freq-sample
	 	chr1 4000 4001 snp G A 0.5
	 	chr1 4001 4002 snp A C,G 0.1,0.5
	 	chr1 4100 4101 del A {} 0.5
	 	chr1 4200 4200 ins {} A 0.5
	 	chr1 4208 4208 ins {} A,AA 0.8,0.9
	}
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test collapsealleles {collapsealleles 2} {
	write_tab tmp/test.tsv {
		chromosome begin end type ref alt freq-sample
	 	chr1 4200 4200 ins {} A 0.5
	 	chr1 4200 4200 ins {} A 0.8
	 	chr1 5000 5001 snp G A 0.5
	}
	exec cg collapsealleles < tmp/test.tsv > tmp/temp.tsv
	write_tab tmp/expected.tsv {
		chromosome begin end type ref alt freq-sample
	 	chr1 4200 4200 ins {} A,A 0.5,0.8
	 	chr1 5000 5001 snp G A 0.5
	}
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test collapsealleles {collapsealleles 3} {
	write_tab tmp/test.tsv {
		chromosome begin end type ref alt zyg-sample1 alleleSeq1-sample1 alleleSeq2-sample1 genotypes-sample1 freq-sample1 zyg-sample2 alleleSeq1-sample2 alleleSeq2-sample2 genotypes-sample2 freq-sample2
	 	chr1 4000 4001 snp G A   t G A 0;1 0.5   t A A 0;0 0.5
	 	chr1 4001 4002 snp A C   c C G 1;2 0.1   o A G 0;2 0.1
	 	chr1 4001 4002 snp A G   c C G 2;1 0.5   t A G 0;1 0.5
	}
	exec cg collapsealleles tmp/test.tsv > tmp/temp.tsv
	write_tab tmp/expected.tsv {
		chromosome begin end type ref alt zyg-sample1 alleleSeq1-sample1 alleleSeq2-sample1 genotypes-sample1 freq-sample1 zyg-sample2 alleleSeq1-sample2 alleleSeq2-sample2 genotypes-sample2 freq-sample2
	 	chr1 4000 4001 snp G A     t G A 0;1 0.5       t A A 0;0 0.5
	 	chr1 4001 4002 snp A C,G   c C G 1;2 0.1,0.5   t A G 0;2 0.1,0.5
	}
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test collapsealleles {splitalleles and collapsealleles without all fields} {
	cg splitalleles data/var-annot.tsv tmp/split.tsv
	cg collapsealleles tmp/split.tsv > tmp/collapsed.tsv
	exec diff data/var-annot.tsv tmp/collapsed.tsv
} {}

test collapsealleles {collapsealleles -duplicates merge} {
	write_tab tmp/test.tsv {
		chromosome begin end type ref alt freq-sample
	 	chr1 4200 4200 snp G A 0.5
	 	chr1 4200 4200 snp G A 0.8
	}
	exec cg collapsealleles -duplicates merge < tmp/test.tsv > tmp/temp.tsv
	write_tab tmp/expected.tsv {
		chromosome begin end type ref alt freq-sample
	 	chr1 4200 4200 snp G A 0.5,0.8
	}
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test collapsealleles {collapsealleles -duplicates keep} {
	write_tab tmp/test.tsv {
		chromosome begin end type ref alt freq-sample
	 	chr1 4200 4200 snp G A 0.5
	 	chr1 4200 4200 snp G A 0.8
	}
	exec cg collapsealleles -duplicates keep < tmp/test.tsv > tmp/temp.tsv
	write_tab tmp/expected.tsv {
		chromosome begin end type ref alt freq-sample
	 	chr1 4200 4200 snp G A,A 0.5,0.8
	}
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test collapsealleles {collapsealleles -duplicates first} {
	write_tab tmp/test.tsv {
		chromosome begin end type ref alt freq-sample
	 	chr1 4200 4200 snp G A 0.5
	 	chr1 4200 4200 snp G A 0.8
	}
	exec cg collapsealleles -duplicates first < tmp/test.tsv > tmp/temp.tsv
	write_tab tmp/expected.tsv {
		chromosome begin end type ref alt freq-sample
	 	chr1 4200 4200 snp G A 0.5
	}
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test collapsealleles {collapsealleles -duplicates max} {
	write_tab tmp/test.tsv {
		chromosome begin end type ref alt freq-sample
	 	chr1 4200 4200 snp G A 0.5
	 	chr1 4200 4200 snp G A 0.8
	}
	exec cg collapsealleles -duplicates {max freq-sample} < tmp/test.tsv > tmp/temp.tsv
	write_tab tmp/expected.tsv {
		chromosome begin end type ref alt freq-sample
	 	chr1 4200 4200 snp G A 0.8
	}
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test collapsealleles {collapsealleles -duplicates min} {
	write_tab tmp/test.tsv {
		chromosome begin end type ref alt freq-sample
	 	chr1 4200 4200 snp G A 0.8
	 	chr1 4200 4200 snp G A 0.5
	}
	exec cg collapsealleles -duplicates {min freq-sample} < tmp/test.tsv > tmp/temp.tsv
	write_tab tmp/expected.tsv {
		chromosome begin end type ref alt freq-sample
	 	chr1 4200 4200 snp G A 0.5
	}
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test collapsealleles {collapsealleles -duplicates error} {
	write_tab tmp/test.tsv {
		chromosome begin end type ref alt freq-sample
	 	chr1 4200 4200 snp G A 0.5
	 	chr1 4200 4200 snp G A 0.8
	}
	exec cg collapsealleles -duplicates error < tmp/test.tsv > tmp/temp.tsv
} {duplicates found:
chr1 4200 4200 snp G A 0.5
chr1 4200 4200 snp G A 0.8} error

test format {long -type analysis} {
	write_tab tmp/wide.tsv {
		var	zyg-gatk-bwa-sample1	zyg-sam-bwa-sample1	zyg-gatk-bwa-sample2
		test	v	u	r
	}
	exec cg long -type analysis tmp/wide.tsv
} {sample	var	zyg
gatk-bwa-sample1	test	v
sam-bwa-sample1	test	u
gatk-bwa-sample2	test	r}

test format {long -type sample} {
	write_tab tmp/wide.tsv {
		var	zyg-gatk-bwa-sample1	zyg-sam-bwa-sample1	zyg-gatk-bwa-sample2
		test	v	u	r
	}
	exec cg long -type sample tmp/wide.tsv
} {sample	var	zyg-gatk-bwa	zyg-sam-bwa
sample1	test	v	u
sample2	test	r	}

test format {long} {
	write_tab tmp/wide.tsv {
		chromosome begin end type ref alt freq-sample1 sequenced-sample1 alleleSeq1-sample1 alleleSeq2-sample1 zyg-sample1 freq-sample2 sequenced-sample2 alleleSeq1-sample2 alleleSeq2-sample2 zyg-sample2
	 	chr1 4200 4200 snp G A 0.5 v G A t 0.8 v A A m
	 	chr1 4200 4200 ins {} A 0.8 v {} A t 0.1 r {} {} r
	 	chr1 5000 5001 snp G T 0.9 v T T m 0.0 r G G r
	}
	exec cg long tmp/wide.tsv tmp/long.tsv
	write_tab tmp/expected.tsv {
		sample chromosome begin end type ref alt sequenced zyg alleleSeq1 alleleSeq2 freq
	 	sample1 chr1 4200 4200 snp G A v t G A 0.5
	 	sample2 chr1 4200 4200 snp G A v m A A 0.8
	 	sample1 chr1 4200 4200 ins {} A v t {} A 0.8
	 	sample2 chr1 4200 4200 ins {} A r r {} {} 0.1
	 	sample1 chr1 5000 5001 snp G T v m T T 0.9
	 	sample2 chr1 5000 5001 snp G T r r G G 0.0
	}
	exec diff tmp/long.tsv tmp/expected.tsv
} {}

test format {long} {
	write_tab tmp/wide.tsv {
		chromosome begin end type ref alt freq-gatk-bwa-sample1 sequenced-gatk-bwa-sample1 alleleSeq1-gatk-bwa-sample1 alleleSeq2-gatk-bwa-sample1 zyg-gatk-bwa-sample1 freq-gatk-bwa-sample2 sequenced-gatk-bwa-sample2 alleleSeq1-gatk-bwa-sample2 alleleSeq2-gatk-bwa-sample2 zyg-gatk-bwa-sample2
	 	chr1 4200 4200 snp G A 0.5 v G A t 0.8 v A A m
	 	chr1 4200 4200 ins {} A 0.8 v {} A t 0.1 r {} {} r
	 	chr1 5000 5001 snp G T 0.9 v T T m 0.0 r G G r
	}
	exec cg long tmp/wide.tsv tmp/long.tsv
	write_tab tmp/expected.tsv {
		sample chromosome begin end type ref alt sequenced zyg alleleSeq1 alleleSeq2 freq
	 	gatk-bwa-sample1 chr1 4200 4200 snp G A v t G A 0.5
	 	gatk-bwa-sample2 chr1 4200 4200 snp G A v m A A 0.8
	 	gatk-bwa-sample1 chr1 4200 4200 ins {} A v t {} A 0.8
	 	gatk-bwa-sample2 chr1 4200 4200 ins {} A r r {} {} 0.1
	 	gatk-bwa-sample1 chr1 5000 5001 snp G T v m T T 0.9
	 	gatk-bwa-sample2 chr1 5000 5001 snp G T r r G G 0.0
	}
	exec diff tmp/long.tsv tmp/expected.tsv
} {}

test format {long from temp} {
	write_tab tmp/test.tsv.temp {
		chromosome begin end type ref alt freq-sample1 sequenced-sample1 alleleSeq1-sample1 alleleSeq2-sample1 zyg-sample1 freq-sample2 sequenced-sample2 alleleSeq1-sample2 alleleSeq2-sample2 zyg-sample2
	 	chr1 4200 4200 snp G A 0.5 v G A t 0.8 v A A m
	 	chr1 4200 4200 ins {} A 0.8 v {} A t 0.1 r {} {} r
	 	chr1 5000 5001 snp G T 0.9 v T T m 0.0 r G G r
	}
	exec cg long tmp/test.tsv.temp tmp/test.tsv
	write_tab tmp/expected.tsv {
		sample chromosome begin end type ref alt sequenced zyg alleleSeq1 alleleSeq2 freq
	 	sample1 chr1 4200 4200 snp G A v t G A 0.5
	 	sample2 chr1 4200 4200 snp G A v m A A 0.8
	 	sample1 chr1 4200 4200 ins {} A v t {} A 0.8
	 	sample2 chr1 4200 4200 ins {} A r r {} {} 0.1
	 	sample1 chr1 5000 5001 snp G T v m T T 0.9
	 	sample2 chr1 5000 5001 snp G T r r G G 0.0
	}
	exec diff tmp/test.tsv tmp/expected.tsv
} {}

test format {long -norm} {
	write_tab tmp/wide.tsv {
		chromosome begin end type ref alt freq-sample1 sequenced-sample1 alleleSeq1-sample1 alleleSeq2-sample1 zyg-sample1 freq-sample2 sequenced-sample2 alleleSeq1-sample2 alleleSeq2-sample2 zyg-sample2
	 	chr1 4200 4200 snp G A 0.5 v G A t 0.8 v A A m
	 	chr1 4200 4200 ins {} A 0.8 v {} A t 0.1 r {} {} r
	 	chr1 5000 5001 snp G T 0.9 v T T m 0.0 r G G r
	}
	exec cg long -norm 1 tmp/wide.tsv tmp/long.tsv
	write_tab tmp/expected.tsv {
		id chromosome begin end type ref alt
	 	1 chr1 4200 4200 snp G A
	 	2 chr1 4200 4200 ins {} A
	 	3 chr1 5000 5001 snp G T
	}
	write_tab tmp/expected.tsv.sampledata.tsv {
		id sample sequenced zyg alleleSeq1 alleleSeq2 freq
	 	1 sample1 v t G A 0.5
	 	1 sample2 v m A A 0.8
	 	2 sample1 v t {} A 0.8
	 	2 sample2 r r {} {} 0.1
	 	3 sample1 v m T T 0.9
	 	3 sample2 r r G G 0.0
	}
	exec diff tmp/long.tsv tmp/expected.tsv
	exec diff tmp/long.tsv.sampledata.tsv tmp/expected.tsv.sampledata.tsv
} {}

test format {long with post, multialt} {
	write_tab tmp/wide.tsv {
		chromosome begin end type ref alt freq-sample1 sequenced-sample1 alleleSeq1-sample1 alleleSeq2-sample1 zyg-sample1 freq-sample2 sequenced-sample2 alleleSeq1-sample2 alleleSeq2-sample2 zyg-sample2 post
	 	chr1 4200 4200 snp G A 0.5 v G A t 0.8 v A A m 1
	 	chr1 4200 4200 ins {} A 0.8 v {} A t 0.1 r {} {} r 2
	 	chr1 5000 5001 snp G T,C 0.9 v T T m 0.0 v C C m 3,4
	}
	write_tab tmp/expected.tsv {
		sample chromosome begin end type ref alt sequenced zyg alleleSeq1 alleleSeq2 freq post
	 	sample1 chr1 4200 4200 snp G A v t G A 0.5 1
	 	sample2 chr1 4200 4200 snp G A v m A A 0.8 1
	 	sample1 chr1 4200 4200 ins {} A v t {} A 0.8 2
	 	sample2 chr1 4200 4200 ins {} A r r {} {} 0.1 2
	 	sample1 chr1 5000 5001 snp G T,C v m T T 0.9 3,4
	 	sample2 chr1 5000 5001 snp G T,C v m C C 0.0 3,4
	}
	exec cg long tmp/wide.tsv tmp/long.tsv
	exec diff tmp/long.tsv tmp/expected.tsv
} {}

test format {long with overlapping field} {
	write_tab tmp/wide.tsv {
		id	validated validated-s1 validated-s2 validated-s3
		1	2	1	0	1
		2	0	0	0	0
	}
	write_tab tmp/expected.tsv {
		sample	id	validated_global	validated
		s1	1	2	1
		s2	1	2	0
		s3	1	2	1
		s1	2	0	0
		s2	2	0	0
		s3	2	0	0
	}
	exec cg long tmp/wide.tsv tmp/long.tsv
	exec diff tmp/long.tsv tmp/expected.tsv
} {}

test format {wide} {
	write_tab tmp/long.tsv {
		chromosome sample begin end type ref alt freq sequenced alleleSeq1 alleleseq2 zyg
	 	chr1 sample1 4200 4200 ins {} A 0.8 v {} A t
	 	chr1 sample2 4200 4200 ins {} A 0.1 r {} {} r
	 	chr1 sample1 4200 4200 snp G A 0.5 v G A t
	 	chr1 sample2 4200 4200 snp G A 0.8 v A A m
	 	chr1 sample1 5000 5001 snp G T 0.9 v T T m
	 	chr1 sample2 5000 5001 snp G T 0.0 v G G r
	}
	exec cg wide tmp/long.tsv tmp/wide.tsv
	write_tab tmp/expected.tsv {
		chromosome begin end type ref alt freq-sample1 sequenced-sample1 alleleSeq1-sample1 alleleseq2-sample1 zyg-sample1 freq-sample2 sequenced-sample2 alleleSeq1-sample2 alleleseq2-sample2 zyg-sample2
	 	chr1 4200 4200 ins {} A 0.8 v {} A t 0.1 r {} {} r
	 	chr1 4200 4200 snp G A 0.5 v G A t 0.8 v A A m
	 	chr1 5000 5001 snp G T 0.9 v T T m 0.0 v G G r
	} {#type	wide tsv
#samplefields	sample}
	exec diff tmp/wide.tsv tmp/expected.tsv
} {}

test format {wide other} {
	write_tab tmp/long.tsv {
		id sample samplea data
	 	1 s1 sa1 i1d11
	 	1 s1 sa2 i1d12
	 	1 s2 sa1 i1d21
	 	2 s1 sa1 i2d11
	 	2 s2 sa1 i2d21
	}
	exec cg wide -s {sample samplea} -f id tmp/long.tsv tmp/wide.tsv
	write_tab tmp/expected.tsv {
		id data-s1-sa1 data-s1-sa2 data-s2-sa1
	 	1 i1d11 i1d12 i1d21
	 	2 i2d11 ? i2d21 
	} {#type	wide tsv
#samplefields	sample samplea}
	exec diff tmp/wide.tsv tmp/expected.tsv
} {}

test format {wide other unsorted} {
	write_tab tmp/long.tsv {
		id sample samplea data
	 	1 s1 sa1 i1d11
	 	2 s1 sa1 i2d11
	 	1 s1 sa2 i1d12
	 	2 s2 sa1 i2d21
	 	1 s2 sa1 i1d21
	}
	exec cg wide -s {sample samplea} -f id tmp/long.tsv tmp/wide.tsv
	write_tab tmp/expected.tsv {
		id data-s1-sa1 data-s1-sa2 data-s2-sa1
	 	1 i1d11 i1d12 i1d21
	 	2 i2d11 ? i2d21 
	} {#type	wide tsv
#samplefields	sample samplea}
	exec diff tmp/wide.tsv tmp/expected.tsv
} {}

test format {keyvalue} {
	write_tab tmp/test.tsv {
		id	v1	v2
	 	id1	1	2
	 	id2	3	3
	}
	exec cg keyvalue tmp/test.tsv tmp/result.tsv
	write_tab tmp/expected.tsv {
		id	key	value
	 	id1	v1	1
	 	id1	v2	2
	 	id2	v1	3
	 	id2	v2	3
	}
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test format {keyvalue} {
	write_tab tmp/test.tsv {
		test	v1	v2
	 	id1	1	2
	 	id2	3	3
	}
	exec cg keyvalue -idfields test tmp/test.tsv tmp/result.tsv
	write_tab tmp/expected.tsv {
		test	key	value
	 	id1	v1	1
	 	id1	v2	2
	 	id2	v1	3
	 	id2	v2	3
	}
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test format {keyvalue -keyname and -valuename} {
	write_tab tmp/test.tsv {
		test	v1	v2
	 	id1	1	2
	 	id2	3	3
	}
	exec cg keyvalue -idfields test -keyname version -valuename count tmp/test.tsv tmp/result.tsv
	write_tab tmp/expected.tsv {
		test	version	count
	 	id1	v1	1
	 	id1	v2	2
	 	id2	v1	3
	 	id2	v2	3
	}
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test format {keyvalue} {
	write_tab tmp/test.tsv {
		test	extra	v1	v2
	 	id1	a	a1	a2
	 	id1	b	b1	b2
	 	id2	a	3	3
	}
	exec cg keyvalue -idfields {test extra} tmp/test.tsv tmp/result.tsv
	write_tab tmp/expected.tsv {
		test	extra	key	value
	 	id1	a	v1	a1
	 	id1	a	v2	a2
	 	id1	b	v1	b1
	 	id1	b	v2	b2
	 	id2	a	v1	3
	 	id2	a	v2	3
	}
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test format {keyvalue} {
	write_tab tmp/wide.tsv {
		chromosome begin end type ref alt freq-test-sample1 sequenced-test-sample1 freq-try-sample1 sequenced-try-sample1
	 	chr1 4200 4200 snp G A 0.5 v 0.8 v
	}
	exec cg keyvalue -samplefields {group sample} tmp/wide.tsv tmp/keyvalue.tsv
	write_tab tmp/expected.tsv {
		group	sample	chromosome	begin	end	type	ref	alt	key	value
		test	sample1	chr1	4200	4200	snp	G	A	freq	0.5
		test	sample1	chr1	4200	4200	snp	G	A	sequenced	v
		try	sample1	chr1	4200	4200	snp	G	A	freq	0.8
		try	sample1	chr1	4200	4200	snp	G	A	sequenced	v
	}
	exec diff tmp/keyvalue.tsv tmp/expected.tsv
} {}

test format {keyvalue} {
	write_tab tmp/wide.tsv {
		chromosome begin end type ref alt freq-sample1 sequenced-sample1 alleleSeq1-sample1 alleleSeq2-sample1 zyg-sample1 freq-sample2 sequenced-sample2 alleleSeq1-sample2 alleleSeq2-sample2 zyg-sample2
	 	chr1 4200 4200 snp G A 0.5 v G A t 0.8 v A A m
	 	chr1 4200 4200 ins {} A 0.8 v {} A t 0.1 r {} {} r
	}
	exec cg keyvalue tmp/wide.tsv tmp/keyvalue.tsv
	write_tab tmp/expected.tsv {
		sample	chromosome	begin	end	type	ref	alt	key	value
		sample1	chr1	4200	4200	snp	G	A	freq	0.5
		sample1	chr1	4200	4200	snp	G	A	sequenced	v
		sample1	chr1	4200	4200	snp	G	A	alleleSeq1	G
		sample1	chr1	4200	4200	snp	G	A	alleleSeq2	A
		sample1	chr1	4200	4200	snp	G	A	zyg	t
		sample2	chr1	4200	4200	snp	G	A	freq	0.8
		sample2	chr1	4200	4200	snp	G	A	sequenced	v
		sample2	chr1	4200	4200	snp	G	A	alleleSeq1	A
		sample2	chr1	4200	4200	snp	G	A	alleleSeq2	A
		sample2	chr1	4200	4200	snp	G	A	zyg	m
		sample1	chr1	4200	4200	ins	{}	A	freq	0.8
		sample1	chr1	4200	4200	ins	{}	A	sequenced	v
		sample1	chr1	4200	4200	ins	{}	A	alleleSeq1	{}
		sample1	chr1	4200	4200	ins	{}	A	alleleSeq2	A
		sample1	chr1	4200	4200	ins	{}	A	zyg	t
		sample2	chr1	4200	4200	ins	{}	A	freq	0.1
		sample2	chr1	4200	4200	ins	{}	A	sequenced	r
		sample2	chr1	4200	4200	ins	{}	A	alleleSeq1	{}
		sample2	chr1	4200	4200	ins	{}	A	alleleSeq2	{}
		sample2	chr1	4200	4200	ins	{}	A	zyg	r
	}
	exec diff tmp/keyvalue.tsv tmp/expected.tsv
} {}

test format {colvalue} {
	write_tab tmp/keyvalue.tsv {
		sample	key	value
		sample1	a	1
		sample1	b	2	
		sample2	a	3
		sample2	b	4	
	}
	exec cg colvalue tmp/keyvalue.tsv tmp/wide.tsv
	write_tab tmp/expected.tsv {
		sample	a	b
		sample1	1	2
		sample2	3	4
	}
	cg tsvdiff tmp/wide.tsv tmp/expected.tsv
} {}

test format {colvalue} {
	write_tab tmp/keyvalue.tsv {
		sample	x	test	nr
		sample1	n	a	1
		sample1	n	b	2	
		sample2	n	a	3
		sample2	n	b	4	
	}
	exec cg colvalue -idfields sample -keyfields test -valuefield nr tmp/keyvalue.tsv tmp/wide.tsv
	write_tab tmp/expected.tsv {
		sample	a	b
		sample1	1	2
		sample2	3	4
	}
	cg tsvdiff tmp/wide.tsv tmp/expected.tsv
} {}

test format {colvalue} {
	write_tab tmp/keyvalue.tsv {
		group	sample	chromosome	begin	end	type	ref	alt	key	value
		test	sample1	chr1	4200	4200	snp	G	A	freq	0.5
		test	sample1	chr1	4200	4200	snp	G	A	sequenced	v
		try	sample1	chr1	4200	4200	snp	G	A	freq	0.8
		try	sample1	chr1	4200	4200	snp	G	A	sequenced	v
	}
	exec cg colvalue -keyfields {key group sample} tmp/keyvalue.tsv tmp/wide.tsv
	write_tab tmp/expected.tsv {
		chromosome begin end type ref alt freq-test-sample1 sequenced-test-sample1 freq-try-sample1 sequenced-try-sample1
	 	chr1 4200 4200 snp G A 0.5 v 0.8 v
	}
	cg tsvdiff tmp/wide.tsv tmp/expected.tsv
} {}

test gene2reg {gene2reg basic} {
	cg select -overwrite 1 -q {$name2 == "HES4"} data/gene_test.tsv tmp/genetmp.tsv
	cg gene2reg tmp/genetmp.tsv tmp/result.tsv
	exec diff tmp/result.tsv data/expected-gene2reg.tsv
} {}

test gene2reg {gene2reg -upstream 100} {
	cg select -overwrite 1 -q {$name2 == "HES4"} data/gene_test.tsv tmp/genetmp.tsv
	cg gene2reg -upstream 100 tmp/genetmp.tsv tmp/result.tsv
	exec diff tmp/result.tsv data/expected-gene2reg.tsv
} {2d1
< chr1	934243	934343	-	downstream	-1	-1	-1	-1	-1	HES4	NM_021170	0	592	HES4	cmpl	cmpl	1,0,0,0,	49
12,13d10
< chr1	935477	935577	-	upstream	-1	-1	-1	-1	-1	HES4	NM_021170	0	592	HES4	cmpl	cmpl	1,0,0,0,	49
< chr1	934243	934343	-	downstream	-1	-1	-1	-1	-1	HES4	NM_001142467	0	592	HES4	cmpl	cmpl	1,0,0,	50
21d17
< chr1	935477	935577	-	upstream	-1	-1	-1	-1	-1	HES4	NM_001142467	0	592	HES4	cmpl	cmpl	1,0,0,	50
child process exited abnormally} error

test gene2reg {gene2reg nocds} {
	cg select -q {$name2 == "HES4"} data/gene_test.tsv tmp/genetmp.tsv
	cg gene2reg -nocds 1 tmp/genetmp.tsv tmp/result.tsv
	cg select -overwrite 1 -f {
		chromosome begin end strand
		{type=if($type in "UTR CDS","RNA",$type)} 
		element rna_start rna_end
		protein_start=-1 protein_end=-1
		*
	} data/expected-gene2reg.tsv tmp/expected-gene2reg.tsv
	exec diff tmp/result.tsv tmp/expected-gene2reg.tsv
} {2c2,3
< chr1	934343	934812	-	RNA	exon4	416	884	-1	-1	HES4	NM_021170	0	592	HES4	cmpl	cmpl	1,0,0,0,	49
---
> chr1	934343	934438	-	RNA	exon4	790	884	-1	-1	HES4	NM_021170	0	592	HES4	cmpl	cmpl	1,0,0,0,	49
> chr1	934438	934812	-	RNA	exon4	416	789	-1	-1	HES4	NM_021170	0	592	HES4	cmpl	cmpl	1,0,0,0,	49
8,9c9,12
< chr1	935245	935477	-	RNA	exon1	0	231	-1	-1	HES4	NM_021170	0	592	HES4	cmpl	cmpl	1,0,0,0,	49
< chr1	934343	934812	-	RNA	exon3	494	962	-1	-1	HES4	NM_001142467	0	592	HES4	cmpl	cmpl	1,0,0,	50
---
> chr1	935245	935353	-	RNA	exon1	124	231	-1	-1	HES4	NM_021170	0	592	HES4	cmpl	cmpl	1,0,0,0,	49
> chr1	935353	935477	-	RNA	exon1	0	123	-1	-1	HES4	NM_021170	0	592	HES4	cmpl	cmpl	1,0,0,0,	49
> chr1	934343	934438	-	RNA	exon3	868	962	-1	-1	HES4	NM_001142467	0	592	HES4	cmpl	cmpl	1,0,0,	50
> chr1	934438	934812	-	RNA	exon3	494	867	-1	-1	HES4	NM_001142467	0	592	HES4	cmpl	cmpl	1,0,0,	50
13c16,17
< chr1	935071	935477	-	RNA	exon1	0	405	-1	-1	HES4	NM_001142467	0	592	HES4	cmpl	cmpl	1,0,0,	50
---
> chr1	935071	935353	-	RNA	exon1	124	405	-1	-1	HES4	NM_001142467	0	592	HES4	cmpl	cmpl	1,0,0,	50
> chr1	935353	935477	-	RNA	exon1	0	123	-1	-1	HES4	NM_001142467	0	592	HES4	cmpl	cmpl	1,0,0,	50
child process exited abnormally} error

test correctvariants {basic} {
	test_cleantmp
	cg select -f {chromosome begin end type ref alt sequenced=$sequenced-sample1 zyg=$zyg-sample1 alleleSeq1=$alleleSeq1-sample1 alleleSeq2=$alleleSeq2-sample1} data/var_h19referrors.tsv tmp/vars.tsv
	exec cg correctvariants -f 1 -s 1 tmp/vars.tsv tmp/temp.tsv $::refseqdir/hg19
	cg select -f {chromosome begin end type ref alt sequenced zyg} tmp/temp.tsv tmp/part.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced	zyg
		1	0	1	snp	N	A	v	c
		1	700	701	snp	N	A	r	o
		1	702	702	ins	{}	A	v	t
		1	1609710	1609711	snp	C	A	v	m
		1	1609720	1609720	ins	{}	A	v	m
		1	2482141	2482142	snp	A	C	r	o
		1	2482141	2482142	snp	A	G	v	m
		15	26517559	26517560	snp	A	T	v	t
	}
	exec diff tmp/part.tsv tmp/expected.tsv
} {}

test correctvariants {basic with doubles} {
	test_cleantmp
	cg select -f {chromosome begin end type ref alt sequenced=$sequenced-sample1 zyg=$zyg-sample1 alleleSeq1=$alleleSeq1-sample1 alleleSeq2=$alleleSeq2-sample1} data/var_h19referrors.tsv tmp/vars.tsv
	set f [open tmp/vars.tsv a]
	puts $f [join {15 26517559 26517560 snp T A v t A T} \t]
	close $f
	exec cg correctvariants -f 1 -s 1 tmp/vars.tsv tmp/temp.tsv $::refseqdir/hg19
	cg select -f {chromosome begin end type ref alt sequenced zyg} tmp/temp.tsv tmp/part.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced	zyg
		1	0	1	snp	N	A	v	c
		1	700	701	snp	N	A	r	o
		1	702	702	ins	{}	A	v	t
		1	1609710	1609711	snp	C	A	v	m
		1	1609720	1609720	ins	{}	A	v	m
		1	2482141	2482142	snp	A	C	r	o
		1	2482141	2482142	snp	A	G	v	m
		15	26517559	26517560	snp	A	T	v	t
	}
	exec diff tmp/part.tsv tmp/expected.tsv
} {}

test correctvariants {basic missing zyg col} {
	test_cleantmp
	cg select -f {chromosome begin end type ref alt sequenced=$sequenced-sample1 alleleSeq1=$alleleSeq1-sample1 alleleSeq2=$alleleSeq2-sample1} data/var_h19referrors.tsv tmp/vars.tsv
	exec cg correctvariants -f 1 -s 1 tmp/vars.tsv tmp/temp.tsv $::refseqdir/hg19
	cg select -f {chromosome begin end type ref alt sequenced} tmp/temp.tsv tmp/part.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced
		1	0	1	snp	N	A	v
		1	700	701	snp	N	A	r
		1	702	702	ins	{}	A	v
		1	1609710	1609711	snp	C	A	v
		1	1609720	1609720	ins	{}	A	v
		1	2482141	2482142	snp	A	C	r
		1	2482141	2482142	snp	A	G	v
		15	26517559	26517560	snp	A	T	v
	}
	exec diff tmp/part.tsv tmp/expected.tsv
} {}

test correctvariants {basic multicompar} {
	test_cleantmp
	exec cg correctvariants -f 1 -s 1 data/var_h19referrors.tsv tmp/temp.tsv $::refseqdir/hg19
	cg select -f {chromosome begin end type ref alt sequenced-* zyg-*} tmp/temp.tsv tmp/part.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	sequenced-sample2	zyg-sample1	zyg-sample2
		1	0	1	snp	N	A	v	r	c	o
		1	700	701	snp	N	A	r	v	o	m
		1	702	702	ins	{}	A	v	r	t	r
		1	1609710	1609711	snp	C	A	v	v	m	c
		1	1609720	1609720	ins	{}	A	v	v	m	t
		1	2482141	2482142	snp	A	C	r	v	o	t
		1	2482141	2482142	snp	A	G	v	r	m	o
		15	26517559	26517560	snp	A	T	v	v	t	t
	}
	exec diff tmp/part.tsv tmp/expected.tsv
} {}

test correctvariants {basic multicompar doubles} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	alleleSeq1-sample1 alleleSeq2-sample1 sequenced-sample2	zyg-sample2	alleleSeq1-sample2 alleleSeq2-sample2
		1	2482141	2482142	snp	A	G	v	m	G	G	v	r	A	A
		1	2482141	2482142	snp	A	G	v	m	G	G	v	r	A	A
		15	26517559	26517560	snp	G	T,A	v	t	T	A	v	o	A	A
		15	26517559	26517560	snp	G	T,A	v	m	T	T	v	o	A	A
	}
	exec cg correctvariants -f 1 -s 0 tmp/temp.tsv tmp/result.tsv $::refseqdir/hg19
	cg select -f {chromosome begin end type ref alt sequenced-* zyg-*} tmp/temp.tsv tmp/part.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	alleleSeq1-sample1 alleleSeq2-sample1 sequenced-sample2	zyg-sample2	alleleSeq1-sample2 alleleSeq2-sample2
		1	2482141	2482142	snp	A	G	v	m	G	G	v	r	A	A
		15	26517559	26517560	snp	A	T	v	t	T	A	r	r	A	A
	}
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test correctvariants {-c 1} {
	test_cleantmp
	exec cg correctvariants -split 0 -c 1 data/updatavartest.tsv tmp/temp.tsv $::refseqdir/hg19
	exec diff tmp/temp.tsv data/expected-updatavartest.tsv
} {}

test correctvariants {-f} {
	test_cleantmp
	exec cg select -f {chromosome begin end type ref alt} data/updatavartest.tsv tmp/temp.tsv
	exec cg correctvariants -split 0 -c 1 tmp/temp.tsv tmp/temp2.tsv $::refseqdir/hg19
	exec diff tmp/temp2.tsv data/expected-updatavartest2.tsv
} {}

test renamesamples {basic} {
	test_cleantmp
	test_genomecombdir
	file copy tmp/test/compar/compar-test.tsv tmp/test/compar/annot_compar-test.tsv
	cg razip tmp/test/compar/compar-test.tsv
	cg zst tmp/test/samples/sample1/var-sample1.tsv
	cg renamesamples tmp/test sample1 new1 sample2 new2 sample3 new3
	file delete -force tmp/test/rename.old
	join [bsort [glob tmp/test/* tmp/test/*/* tmp/test/*/*/*]] \n
} {tmp/test/compar
tmp/test/compar/annot_compar-test.tsv
tmp/test/compar/compar-test.tsv.rz
tmp/test/samples
tmp/test/samples/new1
tmp/test/samples/new1/sreg-new1.tsv
tmp/test/samples/new1/var-new1.tsv.zst
tmp/test/samples/new2
tmp/test/samples/new2/sreg-new2.tsv
tmp/test/samples/new2/var-new2.tsv
tmp/test/samples/new3
tmp/test/samples/new3/prevar-new3.tsv
tmp/test/samples/new3/sreg-new3.tsv
tmp/test/samples/new3/var-new3.tsv}

test renamesamples {rename dir itself} {
	test_cleantmp
	test_genomecombdir
	file copy tmp/test/compar/compar-test.tsv tmp/test/compar/annot_compar-test.tsv
	cg razip tmp/test/compar/compar-test.tsv
	cg zst tmp/test/samples/sample1/var-sample1.tsv
	file delete -force tmp/test/samples/sample2 tmp/test/samples/sample3 tmp/test/compar
	cg renamesamples tmp/test/samples/sample1 sample1 new1
	file delete -force tmp/test/rename.old
	join [bsort [glob tmp/test/samples/* tmp/test/samples/*/*]] \n
} {tmp/test/samples/new1
tmp/test/samples/new1/sreg-new1.tsv
tmp/test/samples/new1/var-new1.tsv.zst
tmp/test/samples/sample1.old
tmp/test/samples/sample1.old/sreg-sample1.tsv
tmp/test/samples/sample1.old/var-sample1.tsv.zst}

test renamesamples {do not follow softlink to dir} {
	test_cleantmp
	test_genomecombdir
	file mkdir tmp/linkeddir
	file copy -force tmp/test/compar/compar-test.tsv tmp/linkeddir/linkeddirfile-test.tsv
	file copy -force tmp/test/samples/sample1 tmp/linkeddir/linkedsample1
	mklink tmp/linkeddir tmp/test/linkeddir
	cg renamesamples tmp/test sample1 new1 sample2 new2 sample3 new3
	file delete -force tmp/test/rename.old
	join [bsort [glob tmp/linkeddir/*/* tmp/test/* tmp/test/*/* tmp/test/*/*/*]] \n
} {tmp/linkeddir/linkedsample1/sreg-sample1.tsv
tmp/linkeddir/linkedsample1/var-sample1.tsv
tmp/test/compar
tmp/test/compar/compar-test.tsv
tmp/test/linkeddir
tmp/test/linkeddir/linkeddirfile-test.tsv
tmp/test/linkeddir/linkedsample1
tmp/test/linkeddir/linkedsample1/sreg-sample1.tsv
tmp/test/linkeddir/linkedsample1/var-sample1.tsv
tmp/test/samples
tmp/test/samples/new1
tmp/test/samples/new1/sreg-new1.tsv
tmp/test/samples/new1/var-new1.tsv
tmp/test/samples/new2
tmp/test/samples/new2/sreg-new2.tsv
tmp/test/samples/new2/var-new2.tsv
tmp/test/samples/new3
tmp/test/samples/new3/prevar-new3.tsv
tmp/test/samples/new3/sreg-new3.tsv
tmp/test/samples/new3/var-new3.tsv}

test tsv2bed {basic} {
	write_tab tmp/test.tsv {
		chromosome begin test end name
	 	chr1 4000 t1 4100 a
	 	chr1 5000 t2 5500 b
	}
	exec cg tsv2bed tmp/test.tsv tmp/temp.bed
	write_tab tmp/expected.bed {
	 	chr1 4000 4100
	 	chr1 5000 5500
	}
	exec diff tmp/temp.bed tmp/expected.bed
} {}

test tsv2bed {pipe 1} {
	write_tab tmp/test.tsv {
		chromosome begin test end name
	 	chr1 4000 t1 4100 a
	 	chr1 5000 t2 5500 b
	}
	exec cg tsv2bed tmp/test.tsv > tmp/temp.bed
	write_tab tmp/expected.bed {
	 	chr1 4000 4100
	 	chr1 5000 5500
	}
	exec diff tmp/temp.bed tmp/expected.bed
} {}

test tsv2bed {pipe 2} {
	write_tab tmp/test.tsv {
		chromosome begin test end name
	 	chr1 4000 t1 4100 a
	 	chr1 5000 t2 5500 b
	}
	exec cg tsv2bed < tmp/test.tsv > tmp/temp.bed
	write_tab tmp/expected.bed {
	 	chr1 4000 4100
	 	chr1 5000 5500
	}
	exec diff tmp/temp.bed tmp/expected.bed
} {}

test tsv2bed {fields} {
	write_tab tmp/test.tsv {
		con begin test end name
	 	chr1 4000 t1 4100 a
	 	chr1 5000 t2 5500 b
	}
	exec cg tsv2bed --fields {con begin end name} < tmp/test.tsv > tmp/temp.bed
	write_tab tmp/expected.bed {
	 	chr1 4000 4100 a
	 	chr1 5000 5500 b
	}
	exec diff tmp/temp.bed tmp/expected.bed
} {}

test tsv2bed {fields poss defaults} {
	write_tab tmp/test.tsv {
		chromosome begin test end name
	 	chr1 4000 t1 4100 a
	 	chr1 5000 t2 5500 b
	}
	exec cg tsv2bed --fields {{} {} {} name} < tmp/test.tsv > tmp/temp.bed
	write_tab tmp/expected.bed {
	 	chr1 4000 4100 a
	 	chr1 5000 5500 b
	}
	exec diff tmp/temp.bed tmp/expected.bed
} {}

test tsv2bed {fields poss other defaults} {
	write_tab tmp/test.tsv {
		chromosome begin test end score name
	 	chr1 4000 t1 4100	1	a
	 	chr1 5000 t2 5500	2	b
	}
	exec cg tsv2bed --fields {{} {} {} {} {}} < tmp/test.tsv > tmp/temp.bed
	write_tab tmp/expected.bed {
	 	chr1 4000 4100 a 1
	 	chr1 5000 5500 b	2
	}
	exec diff tmp/temp.bed tmp/expected.bed
} {}

test tsv2bed {fixed chr} {
	write_tab tmp/test.tsv {
		begin test end name
	 	4000 t1 4100 a
	 	5000 t2 5500 b
	}
	exec cg tsv2bed --fields {chr1 {} {} name} < tmp/test.tsv > tmp/temp.bed
	write_tab tmp/expected.bed {
	 	chr1 4000 4100 a
	 	chr1 5000 5500 b
	}
	exec diff tmp/temp.bed tmp/expected.bed
} {}

test format {tsvjoin} {
	write_tab tmp/file1.tsv {
		id	v1
		id1	a
		id2	b
	}
	write_tab tmp/file2.tsv {
		id	v2
		id1	1
		id2	2
	}
	write_tab tmp/expected.tsv {
		id	v1	v2
		id1	a	1
		id2	b	2
	}
	exec cg tsvjoin tmp/file1.tsv tmp/file2.tsv tmp/result.tsv
	cg tsvdiff tmp/result.tsv tmp/expected.tsv
} {}

test format {tsvjoin -pre2} {
	write_tab tmp/file1.tsv {
		id	v1
		id1	a
		id2	b
	}
	write_tab tmp/file2.tsv {
		id	v2
		id1	1
		id2	2
	}
	write_tab tmp/expected.tsv {
		id	v1	p2_v2
		id1	a	1
		id2	b	2
	}
	exec cg tsvjoin -pre2 p2_ tmp/file1.tsv tmp/file2.tsv tmp/result.tsv
	cg tsvdiff tmp/result.tsv tmp/expected.tsv
} {}

test format {tsvjoin -comments 1} {
	file_write tmp/file1.tsv [deindent {
		# test 1
		id	v1
		id1	a
		id2	b
	}]
	write_tab tmp/file2.tsv {
		# test 2
		id	v2
		id1	1
		id2	2
	}
	write_tab tmp/expected.tsv {
		# test 1
		id	v1	v2
		id1	a	1
		id2	b	2
	}
	exec cg tsvjoin -comments 1 tmp/file1.tsv tmp/file2.tsv tmp/result.tsv
	cg tsvdiff tmp/result.tsv tmp/expected.tsv
} {}

test format {tsvjoin -comments 2} {
	file_write tmp/file1.tsv [deindent {
		# test 1
		id	v1
		id1	a
		id2	b
	}]
	write_tab tmp/file2.tsv {
		# test 2
		id	v2
		id1	1
		id2	2
	}
	write_tab tmp/expected.tsv {
		# test 2
		id	v1	v2
		id1	a	1
		id2	b	2
	}
	exec cg tsvjoin -comments 2 tmp/file1.tsv tmp/file2.tsv tmp/result.tsv
	cg tsvdiff tmp/result.tsv tmp/expected.tsv
} {}

test format {tsvjoin l -type left} {
	file_write tmp/file1.tsv [deindent {
		# test 1
		id	v1
		id1	a
		id2	b
		id3	c
	}]
	file_write tmp/file2.tsv [deindent {
		# test 2
		id	v2
		id1	1
		id3	3
		id4	4
	}]
	file_write tmp/expected.tsv [deindent {
		# test 1
		id	v1	v2
		id1	a	1
		id2	b	
		id3	c	3
	}]
	exec cg tsvjoin -stack 1 -type left tmp/file1.tsv tmp/file2.tsv tmp/result.tsv
	cg tsvdiff tmp/result.tsv tmp/expected.tsv
} {}

test format {tsvjoin l -type right} {
	file_write tmp/file1.tsv [deindent {
		# test 1
		id	v1
		id1	a
		id2	b
		id3	c
	}]
	file_write tmp/file2.tsv [deindent {
		# test 2
		id	v2
		id1	1
		id3	3
		id4	4
	}]
	file_write tmp/expected.tsv [deindent {
		# test 1
		id	v1	v2
		id1	a	1
		id3	c	3
		id4		4
	}]
	exec cg tsvjoin -stack 1 -type right tmp/file1.tsv tmp/file2.tsv tmp/result.tsv
	cg tsvdiff tmp/result.tsv tmp/expected.tsv
} {}

test format {tsvjoin l -type inner} {
	file_write tmp/file1.tsv [deindent {
		# test 1
		id	v1
		id1	a
		id2	b
		id3	c
	}]
	file_write tmp/file2.tsv [deindent {
		# test 2
		id	v2
		id1	1
		id3	3
		id4	4
	}]
	file_write tmp/expected.tsv [deindent {
		# test 1
		id	v1	v2
		id1	a	1
		id3	c	3
	}]
	exec cg tsvjoin -stack 1 -type inner tmp/file1.tsv tmp/file2.tsv tmp/result.tsv
	cg tsvdiff tmp/result.tsv tmp/expected.tsv
} {}

test format {tsvjoin -type full with duplicates} {
	file_write tmp/file1.tsv [deindent {
		# test 1
		id	v1
		id1	a
		id2	b
		id3	c
		id4	d
		id4	db
		id6	f
	}]
	file_write tmp/file2.tsv [deindent {
		# test 2
		id	v2
		id1	1
		id2	2
		id2	2b
		id3	3
		id4	4
		id5	5
	}]
	file_write tmp/expected.tsv [deindent {
		# test 1
		id	v1	v2
		id1	a	1
		id2	b	2
		id2	b	2b
		id3	c	3
		id4	d	4
		id4	db	4
		id5		5
		id6	f	
	}]
	exec cg tsvjoin -stack 1 tmp/file1.tsv tmp/file2.tsv tmp/result.tsv
	cg tsvdiff tmp/result.tsv tmp/expected.tsv
} {}

test format {tsvjoin -type full with duplicates, file2 longer} {
	file_write tmp/file1.tsv [deindent {
		# test 1
		id	v1
		id1	a
		id2	b
		id3	c
		id4	d
		id4	db
	}]
	file_write tmp/file2.tsv [deindent {
		# test 2
		id	v2
		id1	1
		id2	2
		id2	2b
		id3	3
		id4	4
		id5	5
		id6	6
	}]
	file_write tmp/expected.tsv [deindent {
		# test 1
		id	v1	v2
		id1	a	1
		id2	b	2
		id2	b	2b
		id3	c	3
		id4	d	4
		id4	db	4
		id5		5
		id6		6	
	}]
	exec cg tsvjoin -stack 1 tmp/file1.tsv tmp/file2.tsv tmp/result.tsv
	cg tsvdiff tmp/result.tsv tmp/expected.tsv
} {}

test format {tsvjoin unsorted} {
	write_tab tmp/file1.tsv {
		id	v1
		id1	a
		id3	c
		id2	b
	}
	write_tab tmp/file2.tsv {
		id	v2
		id4	4
		id1	1
		id2	2
	}
	write_tab tmp/expected.tsv {
		id	v1	v2
		id1	a	1
		id2	b	2
		id3	c	{}
		id4	{}	4
	}
	exec cg tsvjoin tmp/file1.tsv tmp/file2.tsv tmp/result.tsv
	cg tsvdiff tmp/result.tsv tmp/expected.tsv
} {}

test format {tsvjoin unsorted -sorted 1 error} {
	write_tab tmp/file1.tsv {
		id	v1
		id1	a
		id3	c
		id2	b
	}
	write_tab tmp/file2.tsv {
		id	v2
		id4	4
		id1	1
		id2	2
	}
	write_tab tmp/expected.tsv {
		id	v1	v2
		id1	a	1
		id2	b	2
		id3	c	{}
		id4	{}	4
	}
	exec cg tsvjoin -sorted 1 tmp/file1.tsv tmp/file2.tsv tmp/result.tsv
	cg tsvdiff tmp/result.tsv tmp/expected.tsv
} {file tmp/file1.tsv not properly sorted on id (id3 comes before id2)} error

test csv2tsv {csv2tsv} {
	file_write tmp/temp.csv [deindent {
		a,b,c
		1,2,3
		a1,"b 2","c 3"
		"A 1",B2,"C 3"
		short,"a bit longer","very
		long indeed"
		4,5,6
	}]\n
	file_write tmp/exepected.tsv [deindent {
		a	b	c
		1	2	3
		a1	b 2	c 3
		A 1	B2	C 3
		short	a bit longer	very\nlong indeed
		4	5	6
	}]\n
	cg csv2tsv tmp/temp.csv tmp/result.tsv
	exec diff tmp/result.tsv tmp/exepected.tsv
} {}

test csv2tsv {csv2tsv long} {
	catch {close $o1} ; catch {close $o2}
	set o1 [open tmp/temp.csv w]
	set o2 [open tmp/expected.tsv w]
	for {set i 0} {$i < 10000} {incr i} {
		puts $o1 [deindent {
			a,b,c
			1,2,3
			a1,"b 2","c 3"
			"A 1",B2,"C 3"
			short,"a bit	longer","very
			long indeed"
			4,5,6
		}]
		puts $o2 [deindent {
			a	b	c
			1	2	3
			a1	b 2	c 3
			A 1	B2	C 3
			short	a bit\tlonger	very\nlong indeed
			4	5	6
		}]
	}
	close $o1 ; close $o2
	time {cg csv2tsv tmp/temp.csv tmp/result.tsv}
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test tsv2csv {tsv2csv} {
	file_write tmp/test.tsv [deindent {
		a	b	c
		1	2	3
		a1	b,2	c,3
		A,1	B2	C,3
		short	a bit longer	very\nlong indeed
		4	5	6
	}]\n
	file_write tmp/expected.csv [deindent {
		a,b,c
		1,2,3
		a1,"b,2","c,3"
		"A,1",B2,"C,3"
		short,a bit longer,very\nlong indeed
		4,5,6
	}]\n
	cg tsv2csv tmp/test.tsv tmp/result.csv
	exec diff tmp/result.csv tmp/expected.csv
} {}

test tsv2csv {tsv2csv long} {
	catch {close $o1} ; catch {close $o2}
	set o1 [open tmp/test.tsv w]
	set o2 [open tmp/expected.csv w]
	for {set i 0} {$i < 10000} {incr i} {
		puts $o1 [deindent {
			a	b	c
			1	2	3
			a1	b,2	c,3
			A,1	B2	C,3
			short	a bit longer	very\nlong indeed
			4	5	6
		}]
		puts $o2 [deindent {
			a,b,c
			1,2,3
			a1,"b,2","c,3"
			"A,1",B2,"C,3"
			short,a bit longer,very\nlong indeed
			4,5,6
		}]
	}
	close $o1 ; close $o2
	cg tsv2csv tmp/test.tsv tmp/result.csv
	exec diff tmp/result.csv tmp/expected.csv
} {}

test fasta2tsv {fasta2tsv} {
	file_write tmp/test.fasta [deindent {
		>name1
		AAG
		AAAA
		>name 2
		ACAAAAAAA
	}]\n
	cg fasta2tsv tmp/test.fasta
} {id	sequence
name1	AAGAAAA
name 2	ACAAAAAAA}

test fasta2tsv {fasta2tsv to outfile} {
	file_write tmp/test.fasta [deindent {
		>name1
		AAG
		AAAA
		>name 2
		ACAAAAAAA
	}]\n
	cg fasta2tsv tmp/test.fasta tmp/out.tsv
	file_read tmp/out.tsv
} {id	sequence
name1	AAGAAAA
name 2	ACAAAAAAA
}

test fasta2tsv {fasta2tsv redirect to outfile} {
	file_write tmp/test.fasta [deindent {
		>name1
		AAG
		AAAA
		>name 2
		ACAAAAAAA
	}]\n
	cg fasta2tsv tmp/test.fasta > tmp/out.tsv
	file_read tmp/out.tsv
} {id	sequence
name1	AAGAAAA
name 2	ACAAAAAAA
}

test fasta2tsv {fasta2tsv multiple infiles} {
	file_write tmp/test.fasta [deindent {
		>name1
		AAG
		AAAA
		>name 2
		ACAAAAAAA
	}]\n
	cg fasta2tsv tmp/test.fasta tmp/test.fasta tmp/out.tsv
	file_read tmp/out.tsv
} {id	sequence
name1	AAGAAAA
name 2	ACAAAAAAA
name1	AAGAAAA
name 2	ACAAAAAAA
}

test fasta2tsv {fasta2tsv multiple infiles to stdout} {
	file_write tmp/test.fasta [deindent {
		>name1
		AAG
		AAAA
		>name 2
		ACAAAAAAA
	}]\n
	cg fasta2tsv tmp/test.fasta tmp/test.fasta -
} {id	sequence
name1	AAGAAAA
name 2	ACAAAAAAA
name1	AAGAAAA
name 2	ACAAAAAAA}

test tsv2fasta {tsv2fasta} {
	file_write tmp/test.tsv [deindent {
		id	sequence
		name1	AAGAAAA
		name 2	ACAAAAAAA
	}]\n
	cg tsv2fasta tmp/test.tsv
} {>name1
AAGAAAA
>name 2
ACAAAAAAA}

test tsv2fasta {tsv2fasta to outfile} {
	file_write tmp/test.tsv [deindent {
		id	sequence
		name1	AAGAAAA
		name 2	ACAAAAAAA
	}]\n
	cg tsv2fasta tmp/test.tsv tmp/out.fasta
	file_read tmp/out.fasta
} {>name1
AAGAAAA
>name 2
ACAAAAAAA
}

test tsv2fasta {tsv2fasta redirect to outfile} {
	file_write tmp/test.tsv [deindent {
		id	sequence
		name1	AAGAAAA
		name 2	ACAAAAAAA
	}]\n
	cg tsv2fasta tmp/test.tsv > tmp/out.fasta
	file_read tmp/out.fasta
} {>name1
AAGAAAA
>name 2
ACAAAAAAA
}

test tsv2fasta {tsv2fasta stdin} {
	file_write tmp/test.tsv [deindent {
		id	sequence
		name1	AAGAAAA
		name 2	ACAAAAAAA
	}]\n
	exec cg tsv2fasta < tmp/test.tsv
} {>name1
AAGAAAA
>name 2
ACAAAAAAA}

test sam2tsv {sam2tsv} {
	exec samtools view --no-PG -h data/bwa.bam > tmp/bwa.sam
	cg sam2tsv tmp/bwa.sam tmp/bwa.tsv
	exec diff tmp/bwa.tsv data/bwa.tsv
} {}

test sam2tsv {sam2tsv split extra fields} {
	exec samtools view --no-PG -h data/bwa.bam > tmp/bwa.sam
	cg sam2tsv -fields {NM MD AS XS RG} tmp/bwa.sam tmp/bwa.tsv
	exec diff tmp/bwa.tsv data/bwa_splitextra.tsv
} {}

test sam2tsv {sam2tsv handles bam} {
	cg sam2tsv data/bwa.bam tmp/bwa.tsv
	exec diff tmp/bwa.tsv data/bwa.tsv
} {}

test bam2tsv {bam2tsv} {
	cg bam2tsv data/bwa.bam tmp/bwa.tsv
	exec diff tmp/bwa.tsv data/bwa.tsv
} {}

test sam2tsv {bam2tsv cram} {
	cg bam2cram -refseq $::refseqdir/hg19 data/bwa.bam tmp/bwa.cram
	cg sam2tsv -fields {AS XS MD NM RG} tmp/bwa.cram tmp/bwa.tsv
	exec diff tmp/bwa.tsv data/bwa-fields.tsv
} {}

test sam2tsv {sam2tsv cases} {
	file_write tmp/test.sam [deindent {
		@HD	VN:1.3	SO:coordinate
		@SQ	SN:chr1	LN:249250621
		@SQ	SN:chr2	LN:243199373
		@RG	ID:NA19240m	SM:NA19240m	PL:illumina	PU:NA19240m	LB:solexa-123
		@PG	ID:bwa	PN:bwa	VN:0.7.15-r1140	CL:bwa mem -t 2 -M -R @RG\tID:NA19240m\tSM:NA19240m\tPL:illumina\tPU:NA19240m\tLB:solexa-123 /home/peter/dev/genomecomb/tests/genomecomb.testdata/refseqtest/hg19/genome_hg19.ifas.bwa/genome_hg19.ifas /home/peter/dev/genomecomb/tests/tmp/seq_1.fastq /home/peter/dev/genomecomb/tests/tmp/seq_2.fastq
		r001	99	chr1	7	30	8M2I4M1D3M	=	37	39	TTAGATAAAGGATACTG	*
		r002	0	chr1	9	30	3S6M1P1I4M	*	0	0	AAAAGATAAGGATA	*
		r003	0	chr1	9	30	5S6M	*	0	0	GCCTAAGCTAA	*	SA:Z:chr1,29,-,6H5M,17,0;
		r004	0	chr1	16	30	6M14N5M	*	0	0	ATAGCTTCAGC	*
		r003	2064	chr1	29	17	6H5M	*	0	0	TAGGC	*	SA:Z:chr1,9,+,5S6M,30,1;
		r001	147	chr1	37	30	9M	=	7	-39	CAGCGGCAT	*	NM:i:1
		r005	99	chr1	100	30	4H2S2M2I2M4S2H	*	0	12	GATAAAGGATAG	*
	}]\n
	file_write tmp/expected.tsv [deindent {
		#filetype tsv/samfile
		#fileversion	0.99
		#fields	table
		#fields	field	number	type	description
		#fields	chromosome	1	String	Chromosome/Contig/Reference name
		#fields	begin	1	Integer	Begin of feature (0 based - half open)
		#fields	end	1	Integer	End of feature (0 based - half open)
		#fields	strand	1	String	strand (+ or -)
		#fields	qname	1	String	Query template NAME
		#fields	qstart	1	String	Start of the alignment in the query (0 based - half open)
		#fields	qend	1	String	End of the alignment in the query (0 based - half open)
		#fields	mapquality	1	String	MAPping Quality
		#fields	ref2	1	String	Chromosome/Contig/Reference name of the mate/next read
		#fields	begin2	1	String	Begin of the mate/next read (0 based - half open)
		#fields	strand2	1	String	String	strand (+ or -) of the mate/next read
		#fields	tlen	1	String	observed Template LENgth
		#fields	pair	1	Bool	template having multiple segments in sequencing
		#fields	properpair	1	Bool	each segment properly aligned according to the aligner
		#fields	unmapped	1	Bool	segment unmapped
		#fields	mateunmapped	1	Bool	next segment in the template unmapped
		#fields	read	1	Integer	read: 1 for first segment in the template, 2 for last segment in the template
		#fields	secondary	1	Bool	secondary alignment
		#fields	qcfail	1	Bool	not passing filters, such as platform/vendor quality controls
		#fields	duplicate	1	Bool	PCR or optical duplicate
		#fields	supplementary	1	Bool	supplementary alignment
		#fields	cigar	1	String	CIGAR string
		#fields	seqlen	1	String	sequence length
		#fields	seq	1	String	segment sequence.
		#fields	quality	1	String	ASCII of base quality plus 33
		#fields	other	L	String	
		#@HD	VN:1.3	SO:coordinate
		#@SQ	SN:chr1	LN:249250621
		#@SQ	SN:chr2	LN:243199373
		#@RG	ID:NA19240m	SM:NA19240m	PL:illumina	PU:NA19240m	LB:solexa-123
		#@PG	ID:bwa	PN:bwa	VN:0.7.15-r1140	CL:bwa mem -t 2 -M -R @RG\tID:NA19240m\tSM:NA19240m\tPL:illumina\tPU:NA19240m\tLB:solexa-123 /home/peter/dev/genomecomb/tests/genomecomb.testdata/refseqtest/hg19/genome_hg19.ifas.bwa/genome_hg19.ifas /home/peter/dev/genomecomb/tests/tmp/seq_1.fastq /home/peter/dev/genomecomb/tests/tmp/seq_2.fastq
		chromosome	begin	end	strand	qname	qstart	qend	mapquality	ref2	begin2	strand2	tlen	pair	properpair	unmapped	mateunmapped	read	secondary	qcfail	duplicate	supplementary	cigar	seqlen	seq	quality	other
		chr1	6	22	+	r001	0	17	30	=	36	-	39	1	1	0	0	1	0	0	0	0	8M2I4M1D3M	17	TTAGATAAAGGATACTG	*	RG\tID:NA19240m\tSM:NA19240m\tPL:illumina\tPU:NA19240m\tLB:solexa-123 /home/peter/dev/genomecomb/tests/genomecomb.testdata/refseqtest/hg19/genome_hg19.ifas.bwa/genome_hg19.ifas /home/peter/dev/genomecomb/tests/tmp/seq_1.fastq /home/peter/dev/genomecomb/tests/tmp/seq_2.fastq
		chr1	8	18	+	r002	3	14	30	*	-1	+	0	0	0	0	0	0	0	0	0	0	3S6M1P1I4M	14	AAAAGATAAGGATA	*	
		chr1	8	14	+	r003	5	11	30	*	-1	+	0	0	0	0	0	0	0	0	0	0	5S6M	11	GCCTAAGCTAA	*	SA:Z:chr1,29,-,6H5M,17,0;
		chr1	15	40	+	r004	0	11	30	*	-1	+	0	0	0	0	0	0	0	0	0	0	6M14N5M	11	ATAGCTTCAGC	*	
		chr1	28	33	-	r003	0	5	17	*	-1	+	0	0	0	0	0	0	0	0	0	1	6H5M	5	TAGGC	*	SA:Z:chr1,9,+,5S6M,30,1;
		chr1	36	45	-	r001	0	9	30	=	6	+	-39	1	1	0	0	2	0	0	0	0	9M	9	CAGCGGCAT	*	NM:i:1
		chr1	99	103	+	r005	2	8	30	*	-1	-	12	1	1	0	0	1	0	0	0	0	4H2S2M2I2M4S2H	12	GATAAAGGATAG	*	
	}]\n
	cg sam2tsv tmp/test.sam tmp/test.tsv
	exec diff tmp/test.tsv tmp/expected.tsv
} {}

test tsv2sam {tsv2sam} {
	exec samtools view --no-PG -h data/bwa.bam > tmp/expected.sam
	cg sam2tsv tmp/expected.sam tmp/bwa.tsv
	cg tsv2sam tmp/bwa.tsv tmp/bwa.sam
	exec diff tmp/bwa.sam tmp/expected.sam
} {}

test tsv2bam {tsv2bam} {
	exec samtools view --no-PG -h data/bwa.bam > tmp/expected.sam
	cg sam2tsv tmp/expected.sam tmp/bwa.tsv
	cg tsv2bam tmp/bwa.tsv tmp/bwa.bam
	exec samtools view --no-PG --no-PG -h tmp/bwa.bam > tmp/bwa.sam
	exec diff tmp/bwa.sam tmp/expected.sam
} {}

test tsv2sam {tsv2sam cases} {
	file_write tmp/test.tsv [deindent {
		#@HD	VN:1.3	SO:coordinate
		#@SQ	SN:chr1	LN:249250621
		#@SQ	SN:chr2	LN:243199373
		#@RG	ID:NA19240m	SM:NA19240m	PL:illumina	PU:NA19240m	LB:solexa-123
		#@PG	ID:bwa	PN:bwa	VN:0.7.15-r1140	CL:bwa mem -t 2 -M -R @RG\tID:NA19240m\tSM:NA19240m\tPL:illumina\tPU:NA19240m\tLB:solexa-123 /home/peter/dev/genomecomb/tests/genomecomb.testdata/refseqtest/hg19/genome_hg19.ifas.bwa/genome_hg19.ifas /home/peter/dev/genomecomb/tests/tmp/seq_1.fastq /home/peter/dev/genomecomb/tests/tmp/seq_2.fastq
		chromosome	begin	end	strand	qname	qstart	qend	mapquality	ref2	begin2	strand2	tlen	pair	properpair	unmapped	mateunmapped	read	secondary	qcfail	duplicate	supplementary	cigar	seqlen	seq	quality	other
		chr1	6	22	+	r001	0	17	30	=	36	-	39	1	1	0	0	1	0	0	0	0	8M2I4M1D3M	17	TTAGATAAAGGATACTG	*	
		chr1	8	18	+	r002	3	14	30	*	-1	+	0	0	0	0	0	0	0	0	0	0	3S6M1P1I4M	14	AAAAGATAAGGATA	*	
		chr1	8	14	+	r003	5	11	30	*	-1	+	0	0	0	0	0	0	0	0	0	0	5S6M	11	GCCTAAGCTAA	*	SA:Z:chr1,29,-,6H5M,17,0;
		chr1	15	40	+	r004	0	11	30	*	-1	+	0	0	0	0	0	0	0	0	0	0	6M14N5M	11	ATAGCTTCAGC	*	
		chr1	28	33	-	r003	0	5	17	*	-1	+	0	0	0	0	0	0	0	0	0	1	6H5M	5	TAGGC	*	SA:Z:chr1,9,+,5S6M,30,1;
		chr1	36	45	-	r001	0	9	30	=	6	+	-39	1	1	0	0	2	0	0	0	0	9M	9	CAGCGGCAT	*	NM:i:1
		chr1	99	103	+	r005	2	8	30	*	-1	-	12	1	1	0	0	1	0	0	0	0	4H2S2M2I2M4S2H	12	GATAAAGGATAG	*	
	}]\n
	file_write tmp/expected.sam [deindent {
		@HD	VN:1.3	SO:coordinate
		@SQ	SN:chr1	LN:249250621
		@SQ	SN:chr2	LN:243199373
		@RG	ID:NA19240m	SM:NA19240m	PL:illumina	PU:NA19240m	LB:solexa-123
		@PG	ID:bwa	PN:bwa	VN:0.7.15-r1140	CL:bwa mem -t 2 -M -R @RG\tID:NA19240m\tSM:NA19240m\tPL:illumina\tPU:NA19240m\tLB:solexa-123 /home/peter/dev/genomecomb/tests/genomecomb.testdata/refseqtest/hg19/genome_hg19.ifas.bwa/genome_hg19.ifas /home/peter/dev/genomecomb/tests/tmp/seq_1.fastq /home/peter/dev/genomecomb/tests/tmp/seq_2.fastq
		r001	99	chr1	7	30	8M2I4M1D3M	=	37	39	TTAGATAAAGGATACTG	*
		r002	0	chr1	9	30	3S6M1P1I4M	*	0	0	AAAAGATAAGGATA	*
		r003	0	chr1	9	30	5S6M	*	0	0	GCCTAAGCTAA	*	SA:Z:chr1,29,-,6H5M,17,0;
		r004	0	chr1	16	30	6M14N5M	*	0	0	ATAGCTTCAGC	*
		r003	2064	chr1	29	17	6H5M	*	0	0	TAGGC	*	SA:Z:chr1,9,+,5S6M,30,1;
		r001	147	chr1	37	30	9M	=	7	-39	CAGCGGCAT	*	NM:i:1
		r005	99	chr1	100	30	4H2S2M2I2M4S2H	*	0	12	GATAAAGGATAG	*
	}]\n
	cg tsv2sam tmp/test.tsv tmp/test.sam
	exec diff  tmp/test.sam tmp/expected.sam
	cg tsv2sam tmp/test.tsv tmp/test.bam
	exec samtools view --no-PG -h --no-PG tmp/test.bam > tmp/test2.sam
	exec diff tmp/test2.sam tmp/expected.sam
} {}

test fasta2cramref {cg_fasta2cramref} {
	set o [open tmp/src.fas w]
	puts $o ">r1"
	puts $o "GCCAGCGCAA\nGCCAGCGCAA\nTGGTGGGAGTGAGATGGTG"
	puts $o ">r2"
	puts $o "CAGACATGGACCCTGGCTCCTCCAGTTCTCTGGATCCCTCACTGG"
	close $o
	cg fasta2cramref tmp/src.fas tmp/forcram
	list [bsort [glob tmp/forcram/*]] \
		[file_read tmp/forcram/288894134ca536ab42862c6bc10abe50] \
		[file_read tmp/forcram/9ae0d1e948e7f87aa96f9e1bc6cfe86b] \
		[file_read tmp/forcram/mapping.tsv]
} {{tmp/forcram/9ae0d1e948e7f87aa96f9e1bc6cfe86b tmp/forcram/288894134ca536ab42862c6bc10abe50 tmp/forcram/mapping.tsv} GCCAGCGCAAGCCAGCGCAATGGTGGGAGTGAGATGGTG CAGACATGGACCCTGGCTCCTCCAGTTCTCTGGATCCCTCACTGG {reffile	chromosome	md5	size
src.fas	r1	288894134ca536ab42862c6bc10abe50	39
src.fas	r2	9ae0d1e948e7f87aa96f9e1bc6cfe86b	45
}}

test gff2tsv {basic} {
	file_write tmp/test.gff [string trim [deindent {
		##gff-version 3
		chr1	source1	type1	1	100	.	+	.	ID=test1
		chr1	source1	type1	1000	1100	.	+	.	ID=test1;Name=1
		chr2	source1	type2	11874	14409	.	+	.	ID=test2;Name=test2;description=descr
	}]]\n
	file_write tmp/expected.tsv [string trim [deindent {
		# -- tsv converted from gff, original comments follow --
		##gff-version 3
		# ----
		chromosome	begin	end	type	strand	transcript	gene	score	source	phase	comments	ID	Name	description
		chr1	0	100	type1	+	.	.	.	source1	.		test1		
		chr1	999	1100	type1	+	.	.	.	source1	.		test1	1	
		chr2	11873	14409	type2	+	.	.	.	source1	.		test2	test2	descr
	}]]\n
	cg gff2tsv tmp/test.gff tmp/test.tsv
	cg checktsv tmp/test.tsv
	cg tsvdiff tmp/test.tsv tmp/expected.tsv
} {}

test gtf2tsv {basic} {
	file_write tmp/test.gtf [string trim [deindent {
		chr3	source2	transcript	1000	2000	.	-	.	gene_id "gene1"; transcript_id "transcript1";
		chr3	source2	exon	1000	2000	.	-	5	gene_id "gene1"; transcript_id "transcript1"; exon_number "0";
		chr4	source2	transcript	1000	4000	.	+	.	gene_id "gene2"; transcript_id "transcript2-1";
		chr4	source2	exon	1000	1500	.	+	60	gene_id "gene2"; transcript_id "transcript2-1"; exon_number "0";
		chr4	source2	exon	3000	4000	.	+	60	gene_id "gene2"; transcript_id "transcript2-1"; exon_number "1";
		chr4	source2	gene	1000	5000	.	+	.	gene_id "gene2"; 
		chr4	source2	transcript	1000	5000	.	+	.	gene_id "gene2"; transcript_id "transcript2-2";
		chr4	source2	exon	1000	1500	.	+	60	gene_id "gene2"; transcript_id "transcript2-2"; exon_number "0";
		chr4	source2	exon	3000	4000	.	+	60	gene_id "gene2"; transcript_id "transcript2-2"; exon_number "1";
		chr4	source2	exon	4500	5000	.	+	60	gene_id "gene2"; transcript_id "transcript2-2"; exon_number "2";
	}]]\n
	file_write tmp/expected.tsv [string trim [deindent {
		#filetype	tsv/transcriptsfile
		#fileversion	0.99
		#fields	table
		#fields	field	number	type	description
		#fields	chromosome	1	String	Chromosome name
		#fields	begin	1	Integer	Transcription start position
		#fields	end	1	Integer	Transcription end position
		#fields	name	1	String	Name of transcript (usually transcript_id from GTF)
		#fields	gene	1	String	Alternate name / name of gene (e.g. gene_id from GTF)
		#fields	strand	1	String	+ or - for strand
		#fields	cdsStart	1	Integer	Coding region start
		#fields	cdsEnd	1	Integer	Coding region end
		#fields	exonCount	1	Integer	Number of exons
		#fields	exonStarts	E	Integer	Exon start positions
		#fields	exonEnds	E	Integer	Exon end positions
		#fields	source	1	String	Source of data
		# -- tsv converted from gtf, original comments follow --
		# ----
		chromosome	begin	end	type	transcript	gene	strand	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	source	gene_id	transcript_id	exon_number
		chr3	999	2000	transcript	transcript1	gene1	-			1	999	2000	source2	gene1	transcript1	0
		chr4	999	4000	transcript	transcript2-1	gene2	+			2	999,2999	1500,4000	source2	gene2	transcript2-1	0,1
		chr4	999	5000	transcript	transcript2-2	gene2	+			3	999,2999,4499	1500,4000,5000	source2	gene2	transcript2-2	0,1,2
	}]]\n
	cg gtf2tsv tmp/test.gtf tmp/test.tsv
	cg checktsv tmp/test.tsv
	cg tsvdiff tmp/test.tsv tmp/expected.tsv
} {}

test gtf2tsv {-separate 1} {
	file_write tmp/test.gtf [string trim [deindent {
		chr3	source2	transcript	1000	2000	.	-	.	gene_id "gene1"; transcript_id "transcript1";
		chr3	source2	exon	1000	2000	.	-	5	gene_id "gene1"; transcript_id "transcript1"; exon_number "0";
		chr4	source2	transcript	1000	4000	.	+	.	gene_id "gene2"; transcript_id "transcript2-1";
		chr4	source2	exon	1000	1500	.	+	60	gene_id "gene2"; transcript_id "transcript2-1"; exon_number "0";
		chr4	source2	exon	3000	4000	.	+	60	gene_id "gene2"; transcript_id "transcript2-1"; exon_number "1";
		chr4	source2	transcript	1000	5000	.	+	.	gene_id "gene2"; transcript_id "transcript2-2";
		chr4	source2	exon	1000	1500	.	+	60	gene_id "gene2"; transcript_id "transcript2-2"; exon_number "0";
		chr4	source2	exon	3000	4000	.	+	60	gene_id "gene2"; transcript_id "transcript2-2"; exon_number "1";
		chr4	source2	exon	4500	5000	.	+	60	gene_id "gene2"; transcript_id "transcript2-2"; exon_number "2";
	}]]\n
	file_write tmp/expected.tsv [string trim [deindent {
		#filetype	tsv/transcriptsfile
		#fileversion	0.99
		#fields	table
		#fields	field	number	type	description
		#fields	chromosome	1	String	Chromosome name
		#fields	begin	1	Integer	Transcription start position
		#fields	end	1	Integer	Transcription end position
		#fields	type	1	String	type of element (transcript/exon)
		#fields	transcript	1	String	Name of transcript (usually transcript_id from GTF)
		#fields	gene	1	String	Alternate name / name of gene (e.g. gene_id from GTF)
		#fields	strand	1	String	+ or - for strand
		#fields	source	1	String	Source of data
		#fields	comments	1	String	extra info on element
		# -- tsv converted from gtf, original comments follow --
		# ----
		chromosome	begin	end	type	transcript	gene	strand	source	comments	gene_id	transcript_id	exon_number
		chr3	999	2000	transcript	transcript1	gene1	-	source2		gene1	transcript1	
		chr3	999	2000	exon	transcript1	gene1	-	source2		gene1	transcript1	0
		chr4	999	4000	transcript	transcript2-1	gene2	+	source2		gene2	transcript2-1	
		chr4	999	1500	exon	transcript2-1	gene2	+	source2		gene2	transcript2-1	0
		chr4	2999	4000	exon	transcript2-1	gene2	+	source2		gene2	transcript2-1	1
		chr4	999	5000	transcript	transcript2-2	gene2	+	source2		gene2	transcript2-2	
		chr4	999	1500	exon	transcript2-2	gene2	+	source2		gene2	transcript2-2	0
		chr4	2999	4000	exon	transcript2-2	gene2	+	source2		gene2	transcript2-2	1
		chr4	4499	5000	exon	transcript2-2	gene2	+	source2		gene2	transcript2-2	2
	}]]\n
	cg gtf2tsv -separate 1 -sort 0 tmp/test.gtf | cg select -rf {score phase} > tmp/test.tsv
	cg tsvdiff tmp/test.tsv tmp/expected.tsv
} {}

test gtf2tsv {-separate 1 no attributes} {
	file_write tmp/test.gtf [string trim [deindent {
		chr3	source2	transcript	1000	2000	.	-	.	
	}]]\n
	file_write tmp/expected.tsv [deindent {
		chromosome	begin	end	type	strand	transcript	gene	score	phase	source	comments
		chr3	999	2000	transcript	-	.	.	.	.	source2	
	}]\n
	cg gtf2tsv -separate 1 tmp/test.gtf tmp/test.tsv
	cg tsvdiff tmp/test.tsv tmp/expected.tsv
} {}

test gtf2tsv {ENSG00000187583} {
	cg gtf2tsv data/ENSG00000187583.gtf tmp/test.tsv
	exec diff tmp/test.tsv data/ENSG00000187583.tsv
} {}

test gtf2tsv {gencode} {
	cg gtf2tsv data/gene_gencode_test.gtf tmp/test.tsv
	cg tsv2gtf -addgene 0 tmp/test.tsv tmp/test.gtf
#	set o [open tmp/expected.gtf w]
#	set f [open data/gene_gencode_test.gtf]
#	while {[gets $f line] != -1} {
#		if {[regexp "\tgene\t\|UTR\|stop_codon\|start_codon\|\#\#" $line]} continue
#		puts $o $line
#	}
#	close $o
#	close $f
	exec diff tmp/test.gtf data/gene_gencode_test-converted.gtf
} {}

test tsv2gtf {back and forth} {
	file_write tmp/expected.tsv [deindent {
		chromosome	begin	end	transcript	gene	source	strand	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds
		chr9	27535639	27573494	ENST00000673600.1	C9orf72	intgene	-	27548381	27567120	12	27535639,27548114,27548556,27550649,27556560,27558490,27560226,27561584,27562380,27565530,27566676,27573430	27536286,27548418,27548666,27550707,27556796,27558607,27560299,27561649,27562476,27565590,27567164,27573494
	}]\n
	cg tsv2gtf tmp/expected.tsv tmp/test.gtf
	cg gtf2tsv -stack 1 tmp/test.gtf tmp/result.tsv
	cg select -overwrite 1 -f {chromosome begin end transcript gene source strand cdsStart cdsEnd exonCount exonStarts exonEnds} -rc 1 tmp/result.tsv tmp/ncresult.tsv
	exec diff tmp/ncresult.tsv tmp/expected.tsv
} {}


test tsv2gtf {basic} {
	file_write tmp/test.tsv [string trim [deindent {
		#filetype	tsv/transcriptsfile
		#fileversion	0.99
		#fields	table
		#fields	field	number	type	description
		#fields	chromosome	1	String	Chromosome name
		#fields	begin	1	Integer	Transcription start position
		#fields	end	1	Integer	Transcription end position
		#fields	transcript	1	String	Name of transcript (usually transcript_id from GTF)
		#fields	gene	1	String	Alternate name / name of gene (e.g. gene_id from GTF)
		#fields	strand	1	String	+ or - for strand
		#fields	cdsStart	1	Integer	Coding region start
		#fields	cdsEnd	1	Integer	Coding region end
		#fields	exonCount	1	Integer	Number of exons
		#fields	exonStarts	E	Integer	Exon start positions
		#fields	exonEnds	E	Integer	Exon end positions
		#fields	source	1	String	Source of data
		# -- tsv converted from gtf, original comments follow --
		# ----
		chromosome	begin	end	transcript	gene	strand	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	source	gene_id	transcript_id	exon_number
		chr3	999	2000	transcript1	gene1	-			1	999	2000	source2	gene1	transcript1	0
		chr4	999	4000	transcript2-1	gene2	+			2	999,2999	1500,4000	source2	gene2	transcript2-1	1
		chr4	999	5000	transcript2-2	gene2	+			3	999,2999,4499	1500,4000,5000	source2	gene2	transcript2-2	2
	}]]\n
	file_write tmp/expected.gtf [string trim [deindent {
		chr3	source2	transcript	1000	2000	.	-	.	gene_id "gene1"; transcript_id "transcript1"; gene_name "gene1";
		chr3	source2	exon	1000	2000	.	-	.	gene_id "gene1"; transcript_id "transcript1"; gene_name "gene1";
		chr4	source2	transcript	1000	4000	.	+	.	gene_id "gene2"; transcript_id "transcript2-1"; gene_name "gene2";
		chr4	source2	exon	1000	1500	.	+	.	gene_id "gene2"; transcript_id "transcript2-1"; gene_name "gene2";
		chr4	source2	exon	3000	4000	.	+	.	gene_id "gene2"; transcript_id "transcript2-1"; gene_name "gene2";
		chr4	source2	transcript	1000	5000	.	+	.	gene_id "gene2"; transcript_id "transcript2-2"; gene_name "gene2";
		chr4	source2	exon	1000	1500	.	+	.	gene_id "gene2"; transcript_id "transcript2-2"; gene_name "gene2";
		chr4	source2	exon	3000	4000	.	+	.	gene_id "gene2"; transcript_id "transcript2-2"; gene_name "gene2";
		chr4	source2	exon	4500	5000	.	+	.	gene_id "gene2"; transcript_id "transcript2-2"; gene_name "gene2";
	}]]\n
	cg tsv2gtf tmp/test.tsv tmp/test.gtf
	exec diff tmp/test.gtf tmp/expected.gtf
} {}

test tsv2gtf {basic with CDS} {
	file_write tmp/test.tsv [string trim [deindent {
		#filetype	tsv/transcriptsfile
		#fileversion	0.99
		#fields	table
		#fields	field	number	type	description
		#fields	chromosome	1	String	Chromosome name
		#fields	begin	1	Integer	Transcription start position
		#fields	end	1	Integer	Transcription end position
		#fields	transcript	1	String	Name of transcript (usually transcript_id from GTF)
		#fields	gene	1	String	Alternate name / name of gene (e.g. gene_id from GTF)
		#fields	strand	1	String	+ or - for strand
		#fields	cdsStart	1	Integer	Coding region start
		#fields	cdsEnd	1	Integer	Coding region end
		#fields	exonCount	1	Integer	Number of exons
		#fields	exonStarts	E	Integer	Exon start positions
		#fields	exonEnds	E	Integer	Exon end positions
		#fields	source	1	String	Source of data
		# -- tsv converted from gtf, original comments follow --
		# ----
		chromosome	begin	end	transcript	gene	strand	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	source	gene_id	transcript_id	exon_number
		chr1	923922	943380	testt	test testt	+	924431	943380	12	923922,925921,930154,931038,935771,939039,939274,941143,942135,942409,942558,943252,	924948,926013,930336,931089,935896,939129,939412,941306,942251,942488,943058,943380,	HAVANA	test testt	12
		chr4	999	4000	transcript2-1	gene2	-	1030	2030	2	999,2999	1500,4000	source2	gene2	transcript2-1	1
		chr4	999	5000	transcript2-2	gene2	+	1030	4800	3	999,2999,4499	1500,4000,5000	source2	gene2	transcript2-2	2
	}]]\n
	file_write tmp/expected.gtf [string trim [deindent {
		chr1	HAVANA	transcript	923923	943380	.	+	.	gene_id "test testt"; transcript_id "testt"; gene_name "test testt";
		chr1	HAVANA	exon	923923	924948	.	+	.	gene_id "test testt"; transcript_id "testt"; gene_name "test testt";
		chr1	HAVANA	CDS	924432	924948	.	+	0	gene_id "test testt"; transcript_id "testt"; gene_name "test testt";
		chr1	HAVANA	exon	925922	926013	.	+	.	gene_id "test testt"; transcript_id "testt"; gene_name "test testt";
		chr1	HAVANA	CDS	925922	926013	.	+	2	gene_id "test testt"; transcript_id "testt"; gene_name "test testt";
		chr1	HAVANA	exon	930155	930336	.	+	.	gene_id "test testt"; transcript_id "testt"; gene_name "test testt";
		chr1	HAVANA	CDS	930155	930336	.	+	0	gene_id "test testt"; transcript_id "testt"; gene_name "test testt";
		chr1	HAVANA	exon	931039	931089	.	+	.	gene_id "test testt"; transcript_id "testt"; gene_name "test testt";
		chr1	HAVANA	CDS	931039	931089	.	+	1	gene_id "test testt"; transcript_id "testt"; gene_name "test testt";
		chr1	HAVANA	exon	935772	935896	.	+	.	gene_id "test testt"; transcript_id "testt"; gene_name "test testt";
		chr1	HAVANA	CDS	935772	935896	.	+	1	gene_id "test testt"; transcript_id "testt"; gene_name "test testt";
		chr1	HAVANA	exon	939040	939129	.	+	.	gene_id "test testt"; transcript_id "testt"; gene_name "test testt";
		chr1	HAVANA	CDS	939040	939129	.	+	2	gene_id "test testt"; transcript_id "testt"; gene_name "test testt";
		chr1	HAVANA	exon	939275	939412	.	+	.	gene_id "test testt"; transcript_id "testt"; gene_name "test testt";
		chr1	HAVANA	CDS	939275	939412	.	+	2	gene_id "test testt"; transcript_id "testt"; gene_name "test testt";
		chr1	HAVANA	exon	941144	941306	.	+	.	gene_id "test testt"; transcript_id "testt"; gene_name "test testt";
		chr1	HAVANA	CDS	941144	941306	.	+	2	gene_id "test testt"; transcript_id "testt"; gene_name "test testt";
		chr1	HAVANA	exon	942136	942251	.	+	.	gene_id "test testt"; transcript_id "testt"; gene_name "test testt";
		chr1	HAVANA	CDS	942136	942251	.	+	1	gene_id "test testt"; transcript_id "testt"; gene_name "test testt";
		chr1	HAVANA	exon	942410	942488	.	+	.	gene_id "test testt"; transcript_id "testt"; gene_name "test testt";
		chr1	HAVANA	CDS	942410	942488	.	+	2	gene_id "test testt"; transcript_id "testt"; gene_name "test testt";
		chr1	HAVANA	exon	942559	943058	.	+	.	gene_id "test testt"; transcript_id "testt"; gene_name "test testt";
		chr1	HAVANA	CDS	942559	943058	.	+	1	gene_id "test testt"; transcript_id "testt"; gene_name "test testt";
		chr1	HAVANA	exon	943253	943380	.	+	.	gene_id "test testt"; transcript_id "testt"; gene_name "test testt";
		chr1	HAVANA	CDS	943253	943377	.	+	2	gene_id "test testt"; transcript_id "testt"; gene_name "test testt";
		chr4	source2	transcript	1000	4000	.	-	.	gene_id "gene2"; transcript_id "transcript2-1"; gene_name "gene2";
		chr4	source2	exon	1000	1500	.	-	.	gene_id "gene2"; transcript_id "transcript2-1"; gene_name "gene2";
		chr4	source2	CDS	1034	1500	.	-	1	gene_id "gene2"; transcript_id "transcript2-1"; gene_name "gene2";
		chr4	source2	exon	3000	4000	.	-	.	gene_id "gene2"; transcript_id "transcript2-1"; gene_name "gene2";
		chr4	source2	CDS	3000	4000	.	-	0	gene_id "gene2"; transcript_id "transcript2-1"; gene_name "gene2";
		chr4	source2	transcript	1000	5000	.	+	.	gene_id "gene2"; transcript_id "transcript2-2"; gene_name "gene2";
		chr4	source2	exon	1000	1500	.	+	.	gene_id "gene2"; transcript_id "transcript2-2"; gene_name "gene2";
		chr4	source2	CDS	1031	1500	.	+	0	gene_id "gene2"; transcript_id "transcript2-2"; gene_name "gene2";
		chr4	source2	exon	3000	4000	.	+	.	gene_id "gene2"; transcript_id "transcript2-2"; gene_name "gene2";
		chr4	source2	CDS	3000	4000	.	+	1	gene_id "gene2"; transcript_id "transcript2-2"; gene_name "gene2";
		chr4	source2	exon	4500	5000	.	+	.	gene_id "gene2"; transcript_id "transcript2-2"; gene_name "gene2";
		chr4	source2	CDS	4500	4797	.	+	2	gene_id "gene2"; transcript_id "transcript2-2"; gene_name "gene2";
	}]]\n
	cg tsv2gtf -stack 1 tmp/test.tsv tmp/test.gtf
	exec diff tmp/test.gtf tmp/expected.gtf
} {}

test tsv2gtf {basic with CDS bug fix} {
	file_write tmp/test.tsv [string trim [deindent {
		#filetype	tsv/transcriptsfile
		#fileversion	0.99
		#fields	table
		#fields	field	number	type	description
		#fields	chromosome	1	String	Chromosome name
		#fields	begin	1	Integer	Transcription start position
		#fields	end	1	Integer	Transcription end position
		#fields	transcript	1	String	Name of transcript (usually transcript_id from GTF)
		#fields	gene	1	String	Alternate name / name of gene (e.g. gene_id from GTF)
		#fields	strand	1	String	+ or - for strand
		#fields	cdsStart	1	Integer	Coding region start
		#fields	cdsEnd	1	Integer	Coding region end
		#fields	exonCount	1	Integer	Number of exons
		#fields	exonStarts	E	Integer	Exon start positions
		#fields	exonEnds	E	Integer	Exon end positions
		#fields	source	1	String	Source of data
		# -- tsv converted from gtf, original comments follow --
		# ----
		chromosome	begin	end	transcript	gene	strand	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	source	gene_id	transcript_id	exon_number
		chr9	34551431	34589706	test2	CNTFR	-	34552070	34568981	9	34551431,34552160,34552673,34556254,34557525,34557866,34564598,34568896,34589554,	34552071,34552329,34552854,34556418,34557692,34557984,34564832,34568981,34589706,	HAVANA	CNTFR	test2	9
	}]]\n
	file_write tmp/expected.gtf [string trim [deindent {
		chr9	HAVANA	transcript	34551432	34589706	.	-	.	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	exon	34551432	34552071	.	-	.	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	exon	34552161	34552329	.	-	.	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	CDS	34552164	34552329	.	-	2	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	exon	34552674	34552854	.	-	.	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	CDS	34552674	34552854	.	-	0	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	exon	34556255	34556418	.	-	.	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	CDS	34556255	34556418	.	-	2	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	exon	34557526	34557692	.	-	.	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	CDS	34557526	34557692	.	-	1	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	exon	34557867	34557984	.	-	.	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	CDS	34557867	34557984	.	-	2	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	exon	34564599	34564832	.	-	.	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	CDS	34564599	34564832	.	-	2	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	exon	34568897	34568981	.	-	.	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	CDS	34568897	34568981	.	-	0	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	exon	34589555	34589706	.	-	.	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
	}]]\n
	cg tsv2gtf -stack 1 tmp/test.tsv tmp/test.gtf
	exec diff tmp/test.gtf tmp/expected.gtf
} {}

test tsv2gtf {basic with CDS bug fix check around} {
	file_write tmp/test.tsv [string trim [deindent {
		chromosome	begin	end	transcript	gene	strand	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	source	gene_id	transcript_id	exon_number
		chr9	34551431	34552854	test2	CNTFR	-	34552068	34552850	3	34551431,34552160,34552673,	34552071,34552329,34552854,	HAVANA	CNTFR	test2	3
	}]]\n
	file_write tmp/expected.gtf [string trim [deindent {
		chr9	HAVANA	transcript	34551432	34552854	.	-	.	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	exon	34551432	34552071	.	-	.	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	exon	34552161	34552329	.	-	.	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	CDS	34552162	34552329	.	-	0	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	exon	34552674	34552854	.	-	.	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	CDS	34552674	34552850	.	-	0	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
	}]]\n
	cg tsv2gtf -stack 1 tmp/test.tsv tmp/test.gtf
	exec diff tmp/test.gtf tmp/expected.gtf
	#
	file_write tmp/test.tsv [string trim [deindent {
		chromosome	begin	end	transcript	gene	strand	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	source	gene_id	transcript_id	exon_number
		chr9	34551431	34552854	test2	CNTFR	-	34552067	34552850	3	34551431,34552160,34552673,	34552071,34552329,34552854,	HAVANA	CNTFR	test2	3
	}]]\n
	file_write tmp/expected.gtf [string trim [deindent {
		chr9	HAVANA	transcript	34551432	34552854	.	-	.	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	exon	34551432	34552071	.	-	.	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	exon	34552161	34552329	.	-	.	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	CDS	34552161	34552329	.	-	0	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	exon	34552674	34552854	.	-	.	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	CDS	34552674	34552850	.	-	0	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
	}]]\n
	cg tsv2gtf -stack 1 tmp/test.tsv tmp/test.gtf
	exec diff tmp/test.gtf tmp/expected.gtf
	#
	file_write tmp/test.tsv [string trim [deindent {
		chromosome	begin	end	transcript	gene	strand	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	source	gene_id	transcript_id	exon_number
		chr9	34551431	34552854	test2	CNTFR	-	34552066	34552850	3	34551431,34552160,34552673,	34552071,34552329,34552854,	HAVANA	CNTFR	test2	3
	}]]\n
	file_write tmp/expected.gtf [string trim [deindent {
		chr9	HAVANA	transcript	34551432	34552854	.	-	.	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	exon	34551432	34552071	.	-	.	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	CDS	34552070	34552071	.	-	2	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	exon	34552161	34552329	.	-	.	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	CDS	34552161	34552329	.	-	0	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	exon	34552674	34552854	.	-	.	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
		chr9	HAVANA	CDS	34552674	34552850	.	-	0	gene_id "CNTFR"; transcript_id "test2"; gene_name "CNTFR";
	}]]\n
	cg tsv2gtf -stack 1 tmp/test.tsv tmp/test.gtf
	exec diff tmp/test.gtf tmp/expected.gtf
	#
	file_write tmp/test.tsv [string trim [deindent {
		chromosome	begin	end	transcript	gene	strand	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	source	gene_id	transcript_id	exon_number
		chr1	1000	2000	test1	testg	+	1020	1600	3	1000,1500,1800,	1100,1600,2000,	test	testg	test1	3
		chr1	1000	2000	test2	testg	+	1020	1800	3	1000,1500,1800,	1100,1600,2000,	test	testg	test2	3
		chr1	1000	2000	test2	testg	+	1020	1801	3	1000,1500,1800,	1100,1600,2000,	test	testg	test2	3
		chr1	1000	2000	test2	testg	+	1020	1802	3	1000,1500,1800,	1100,1600,2000,	test	testg	test2	3
		chr1	1000	2000	test2	testg	+	1020	1803	3	1000,1500,1800,	1100,1600,2000,	test	testg	test2	3
		chr1	1000	2000	test2	testg	+	1020	1804	3	1000,1500,1800,	1100,1600,2000,	test	testg	test2	3
	}]]\n
	file_write tmp/expected.gtf [string trim [deindent {
		chr1	test	transcript	1001	2000	.	+	.	gene_id "testg"; transcript_id "test1"; gene_name "testg";
		chr1	test	exon	1001	1100	.	+	.	gene_id "testg"; transcript_id "test1"; gene_name "testg";
		chr1	test	CDS	1021	1100	.	+	0	gene_id "testg"; transcript_id "test1"; gene_name "testg";
		chr1	test	exon	1501	1600	.	+	.	gene_id "testg"; transcript_id "test1"; gene_name "testg";
		chr1	test	CDS	1501	1597	.	+	1	gene_id "testg"; transcript_id "test1"; gene_name "testg";
		chr1	test	exon	1801	2000	.	+	.	gene_id "testg"; transcript_id "test1"; gene_name "testg";
		chr1	test	transcript	1001	2000	.	+	.	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	exon	1001	1100	.	+	.	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	CDS	1021	1100	.	+	0	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	exon	1501	1600	.	+	.	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	CDS	1501	1597	.	+	1	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	exon	1801	2000	.	+	.	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	transcript	1001	2000	.	+	.	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	exon	1001	1100	.	+	.	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	CDS	1021	1100	.	+	0	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	exon	1501	1600	.	+	.	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	CDS	1501	1598	.	+	1	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	exon	1801	2000	.	+	.	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	transcript	1001	2000	.	+	.	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	exon	1001	1100	.	+	.	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	CDS	1021	1100	.	+	0	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	exon	1501	1600	.	+	.	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	CDS	1501	1599	.	+	1	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	exon	1801	2000	.	+	.	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	transcript	1001	2000	.	+	.	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	exon	1001	1100	.	+	.	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	CDS	1021	1100	.	+	0	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	exon	1501	1600	.	+	.	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	CDS	1501	1600	.	+	1	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	exon	1801	2000	.	+	.	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	transcript	1001	2000	.	+	.	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	exon	1001	1100	.	+	.	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	CDS	1021	1100	.	+	0	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	exon	1501	1600	.	+	.	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	CDS	1501	1600	.	+	1	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	exon	1801	2000	.	+	.	gene_id "testg"; transcript_id "test2"; gene_name "testg";
		chr1	test	CDS	1801	1801	.	+	0	gene_id "testg"; transcript_id "test2"; gene_name "testg";
	}]]\n
	cg tsv2gtf -stack 1 tmp/test.tsv tmp/test.gtf
	exec diff tmp/test.gtf tmp/expected.gtf
} {}

test tsv2gtf {basic with CDS bug fix empty CDS} {
	file_write tmp/test.tsv [string trim [deindent {
		chrom	start	end	strand	geneid	transcript	score	bin	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	name2	cdsStartStat	cdsEndStat	exonFrames	ROW
		chr1	3857603	3883741	+	DFFB	ENST00000341385.4	0	614	3857603	3872484	4	3857603,3858717,3872483,3883506,	3857717,3858851,3872572,3883741,	DFFB	cmpl	cmpl	0,0,2,-1,	328
	}]]\n
	file_write tmp/expected.gtf [string trim [deindent {
		chr1	genepred2tsv	transcript	3857604	3883741	0	+	.	gene_id "DFFB"; transcript_id "ENST00000341385.4"; gene_name "DFFB";
		chr1	genepred2tsv	exon	3857604	3857717	0	+	.	gene_id "DFFB"; transcript_id "ENST00000341385.4"; gene_name "DFFB";
		chr1	genepred2tsv	CDS	3857604	3857717	0	+	0	gene_id "DFFB"; transcript_id "ENST00000341385.4"; gene_name "DFFB";
		chr1	genepred2tsv	exon	3858718	3858851	0	+	.	gene_id "DFFB"; transcript_id "ENST00000341385.4"; gene_name "DFFB";
		chr1	genepred2tsv	CDS	3858718	3858849	0	+	0	gene_id "DFFB"; transcript_id "ENST00000341385.4"; gene_name "DFFB";
		chr1	genepred2tsv	exon	3872484	3872572	0	+	.	gene_id "DFFB"; transcript_id "ENST00000341385.4"; gene_name "DFFB";
		chr1	genepred2tsv	exon	3883507	3883741	0	+	.	gene_id "DFFB"; transcript_id "ENST00000341385.4"; gene_name "DFFB";
	}]]\n
	cg tsv2gtf -stack 1 tmp/test.tsv tmp/test.gtf
	exec diff tmp/test.gtf tmp/expected.gtf
} {}

test tsv2gtf {bug fix} {
	file_write tmp/test.tsv [string trim [deindent {
		chromosome	begin	end	transcript	gene	strand	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	source	gene_type	havana_transcript	transcript_name	exon_number	gene_id	tag	transcript_support_level	gene_name	level	transcript_type	transcript_id	exon_id	ont	hgnc_id	havana_gene	protein_id	ccdsid
		chr4	573879	576295	ENST00000610212.3	ENSG00000273238.3	-	574904	576062	1	573879,	576295,	HAVANA	protein_coding	OTTHUMT00000492917.2	TMEM271-201	1	ENSG00000273238.3	CAGE_supported_TSS,basic,Ensembl_canonical,MANE_Select,appris_principal_1	NA	TMEM271	2	protein_coding	ENST00000610212.3	ENSE00003703066.3		HGNC:53639	OTTHUMG00000192218.2	ENSP00000493161.1	
	}]]\n
	file_write tmp/expected.gtf [string trim [deindent {
		chr4	HAVANA	transcript	573880	576295	.	-	.	gene_id "ENSG00000273238.3"; transcript_id "ENST00000610212.3"; gene_name "TMEM271"; gene_type "protein_coding"; transcript_type "protein_coding"; transcript_name "TMEM271-201"; level "2"; transcript_support_level "NA"; hgnc_id "HGNC:53639"; tag "CAGE_supported_TSS,basic,Ensembl_canonical,MANE_Select,appris_principal_1"; havana_gene "OTTHUMG00000192218.2"; havana_transcript "OTTHUMT00000492917.2";
		chr4	HAVANA	exon	573880	576295	.	-	.	gene_id "ENSG00000273238.3"; transcript_id "ENST00000610212.3"; gene_name "TMEM271"; gene_type "protein_coding"; transcript_type "protein_coding"; transcript_name "TMEM271-201"; level "2"; transcript_support_level "NA"; hgnc_id "HGNC:53639"; tag "CAGE_supported_TSS,basic,Ensembl_canonical,MANE_Select,appris_principal_1"; havana_gene "OTTHUMG00000192218.2"; havana_transcript "OTTHUMT00000492917.2";
		chr4	HAVANA	CDS	574908	576062	.	-	0	gene_id "ENSG00000273238.3"; transcript_id "ENST00000610212.3"; gene_name "TMEM271"; gene_type "protein_coding"; transcript_type "protein_coding"; transcript_name "TMEM271-201"; level "2"; transcript_support_level "NA"; hgnc_id "HGNC:53639"; tag "CAGE_supported_TSS,basic,Ensembl_canonical,MANE_Select,appris_principal_1"; havana_gene "OTTHUMG00000192218.2"; havana_transcript "OTTHUMT00000492917.2";
	}]]\n
	cg tsv2gtf -stack 1 tmp/test.tsv tmp/test.gtf
	exec diff tmp/test.gtf tmp/expected.gtf
} {}


test convert_pipe {convert_pipe various} {
	set errors {}
	foreach {test expected} {
		{convert_pipe a.tsv.zst dg.bam -endpipe 1} {zstd-mt -T 1 -k -q -d -c a.tsv.zst | cg tsv2sam > dg.bam}
		{convert_pipe a.bam dg.tsv.gz -endpipe 1} {cg sam2tsv a.bam | bgzip -l 6 -c > dg.tsv.gz}
		{convert_pipe a.vcf dg.vcf.gz -endpipe 1} {bgzip -l 6 -c a.vcf > dg.vcf.gz}
		{convert_pipe a.vcf dg.vcf.gz} {bgzip -l 6 -c a.vcf}
		{convert_pipe a.vcf dg.bed} {cg vcf2tsv a.vcf | cg tsv2bed}
		{convert_pipe -.vcf -.bed} {cg vcf2tsv | cg tsv2bed}
		{convert_pipe a.vcf.gz dg.bed.zst} {zcat a.vcf.gz | cg vcf2tsv | cg tsv2bed | zstd-mt -k -q -8 -b 0.5 -T 1 -c}
		{convert_pipe a.tsv b.tsv} {}
		{convert_pipe a.tsv b.tsv -endpipe 1} {cp a.tsv b.tsv}
		{convert_pipe a.tsv b.tsv -endpipe 1 -cpcmd {ln -s}} {ln -s a.tsv b.tsv}
		{convert_pipe -.tsv b.tsv -endpipe 1} {> b.tsv}
		{convert_pipe a.tsv -.tsv -endpipe 1} {cat a.tsv}
	} {
		set result [{*}$test]
		if {$result ne $expected} {
			lappend errors "error testing \"$test\": gave \"$result\" instead of \"$expected\""
		}
	}
	join $errors \n
} {}

test convert {gb2fasta} {
	file delete tmp/results.fas
	cg gb2fasta -clean_names 0 data/plasmids.gb tmp/results.fas
	exec diff tmp/results.fas data/plasmids.fas
} {}

test convert {gb2fasta redirect} {
	file delete tmp/results.fas
	cg gb2fasta -clean_names 0 data/plasmids.gb > tmp/results.fas
	exec diff tmp/results.fas data/plasmids.fas
} {}

test 10x {tsv210x and 10x2tsv} {
	set tsvfile data/sc_gene_counts-test.tsv
	set 10xdir tmp/test.10x
	set genefields {geneid gene gene_id} ; set cellbarcodefield {} ; set countfield {}
	cg tsv210x data/sc_gene_counts-test.tsv tmp/test.10x
	cg 10x2tsv tmp/test.10x tmp/sc_gene_counts-result.tsv
	exec diff tmp/sc_gene_counts-result.tsv data/sc_gene_counts-test.tsv
} {}

test bam2cram {bam2cram} {
	exec samtools view -b -h data/bwa.sam > tmp/bwa.bam
	set refseq $::refseqdir/hg19
	cg bam2cram -refseq $refseq tmp/bwa.bam tmp/bwa.cram
	cg sam2tsv -fields {NM MD AS XS RG MQ MC ms} data/bwa.sam | cg select -rc 1 > tmp/test.tsv
	cg sam2tsv -fields {NM MD AS XS RG MQ MC ms} tmp/bwa.cram | cg select -rc 1 > tmp/bwa.tsv
	cg tsvdiff tmp/bwa.tsv tmp/test.tsv
} {}

test bam2cram {bams2crams} {
	exec samtools view -b -h data/bwa.sam > tmp/bwa.bam
	exec samtools view -b -h data/bwa.sam > tmp/bwa2.bam
	set refseq $::refseqdir/hg19
	cg bams2crams -refseq $refseq tmp/bwa.bam tmp/bwa2.bam
	cg sam2tsv -fields {NM MD AS XS RG MQ MC ms} data/bwa.sam | cg select -rc 1 > tmp/test.tsv
	cg sam2tsv -fields {NM MD AS XS RG MQ MC ms} tmp/bwa.cram | cg select -rc 1 > tmp/bwa.tsv
	cg sam2tsv -fields {NM MD AS XS RG MQ MC ms} tmp/bwa2.cram | cg select -rc 1 > tmp/bwa2.tsv
	cg tsvdiff tmp/bwa.tsv tmp/test.tsv
	cg tsvdiff tmp/bwa2.tsv tmp/test.tsv
} {}

test bam2cram {bams2crams -links rename} {
	exec samtools view -b -h data/bwa.sam > tmp/bwa.bam
	exec samtools index tmp/bwa.bam
	mklink tmp/bwa.bam tmp/bwa2.bam
	mklink tmp/bwa.bam.bai tmp/bwa2.bam.bai
	set refseq $::refseqdir/hg19
	cg bams2crams -links rename -refseq $refseq tmp/bwa.bam tmp/bwa2.bam
	cg sam2tsv -fields {NM MD AS XS RG MQ MC ms} data/bwa.sam | cg select -rc 1 > tmp/test.tsv
	cg sam2tsv -fields {NM MD AS XS RG MQ MC ms} tmp/bwa.cram | cg select -rc 1 > tmp/bwa.tsv
	cg sam2tsv -fields {NM MD AS XS RG MQ MC ms} tmp/bwa2.cram | cg select -rc 1 > tmp/bwa2.tsv
	cg tsvdiff tmp/bwa.tsv tmp/test.tsv
	cg tsvdiff tmp/bwa2.tsv tmp/test.tsv
} {}

testsummarize
