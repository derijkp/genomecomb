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

test_cleantmp

test vcf2tsv {vcf2tsv} {
	exec cg vcf2tsv data/test1000glow.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test1000glow.vcf2tsv
} {}

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
	cg select -q {$name2 == "HES4"} data/gene_test.tsv tmp/genetmp.tsv
	cg gene2reg tmp/genetmp.tsv tmp/result.tsv
	exec diff tmp/result.tsv data/expected-gene2reg.tsv
} {}

test correctvariants {basic} {
	test_cleantmp
	cg select -f {chromosome begin end type ref alt sequenced=$sequenced-sample1 zyg=$zyg-sample1 alleleSeq1=$alleleSeq1-sample1 alleleSeq2=$alleleSeq2-sample1} data/var_h19referrors.tsv tmp/vars.tsv
	exec cg correctvariants -f 1 -s 1 tmp/vars.tsv tmp/temp.tsv /complgen/refseq/hg19
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
	exec cg correctvariants -f 1 -s 1 tmp/vars.tsv tmp/temp.tsv /complgen/refseq/hg19
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
	exec cg correctvariants -f 1 -s 1 tmp/vars.tsv tmp/temp.tsv /complgen/refseq/hg19
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
	exec cg correctvariants -f 1 -s 1 data/var_h19referrors.tsv tmp/temp.tsv /complgen/refseq/hg19
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
	exec cg correctvariants -f 1 -s 0 tmp/temp.tsv tmp/result.tsv /complgen/refseq/hg19
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
	exec cg correctvariants -c 1 data/updatavartest.tsv tmp/temp.tsv /complgen/refseq/hg18
	exec diff tmp/temp.tsv data/expected-updatavartest.tsv
} {}

test correctvariants {-f} {
	test_cleantmp
	exec cg select -f {chromosome begin end type ref alt} data/updatavartest.tsv tmp/temp.tsv
	exec cg correctvariants -c 1 tmp/temp.tsv tmp/temp2.tsv /complgen/refseq/hg18
	exec diff tmp/temp2.tsv data/expected-updatavartest2.tsv
} {}

test renamesamples {basic} {
	test_cleantmp
	test_genomecombdir
	file copy tmp/test/compar/compar-test.tsv tmp/test/compar/annot_compar-test.tsv
	cg razip tmp/test/compar/compar-test.tsv
	cg lz4 tmp/test/samples/sample1/var-sample1.tsv
	cg renamesamples tmp/test sample1 new1 sample2 new2 sample3 new3
	lsort -dict [glob tmp/test/* tmp/test/*/* tmp/test/*/*/*]
} {tmp/test/compar tmp/test/compar/annot_compar-test.tsv tmp/test/compar/compar-test.tsv.index tmp/test/compar/compar-test.tsv.index/multicompar tmp/test/compar/compar-test.tsv.rz tmp/test/samples tmp/test/samples/new1 tmp/test/samples/new1/sreg-new1.tsv tmp/test/samples/new1/var-new1.tsv.lz4 tmp/test/samples/new2 tmp/test/samples/new2/sreg-new2.tsv tmp/test/samples/new2/var-new2.tsv tmp/test/samples/new3 tmp/test/samples/new3/prevar-new3.tsv tmp/test/samples/new3/sreg-new3.tsv tmp/test/samples/new3/var-new3.tsv}

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

testsummarize
