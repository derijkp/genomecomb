#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test cg2tsv {cg2tsv} {
	exec cg cg2tsv -split 0 data/var-cgtest.tsv data/gene-cgtest.tsv tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-cgtest.tsv
} {}

test cg2tsv {cg2tsv without genefile} {
	exec cg cg2tsv -split 0 data/var-cgtest.tsv tmp/temp.tsv
	exec cg select -f {locus chromosome begin end type reference alt zyg alleleSeq1 alleleSeq2 totalScore1 totalScore2 xRef} data/expected-cgtest.tsv tmp/expected.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test cg2tsv {cg2tsv (var2annot)} {
	exec cg var2annot -split 0 data/var-cgtest.tsv data/gene-cgtest.tsv tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-cgtest.tsv
} {}

test cg2tsv {cg2tsv -split 1} {
	exec cg cg2tsv -split 1 data/var-cgtest.tsv data/gene-cgtest.tsv tmp/temp.tsv
	catch {exec diff tmp/temp.tsv data/expected-cgtest.tsv} result
	set result
} {2c2
< #split	1
---
> #split	0
21,22c21
< 15950476	chr21	9438355	9438355	ins		CCCCAACGCCGCGGCTTTTTGT	c	CCCCAACGCCGCGGCTTTTTGT	CCCCAACGCCGCGGGTTTTTCT	39	262	VQLOW															
< 15950476	chr21	9438355	9438355	ins		CCCCAACGCCGCGGGTTTTTCT	c	CCCCAACGCCGCGGCTTTTTGT	CCCCAACGCCGCGGGTTTTTCT	39	262	VQLOW															
---
> 15950476	chr21	9438355	9438355	ins		CCCCAACGCCGCGGCTTTTTGT,CCCCAACGCCGCGGGTTTTTCT	c	CCCCAACGCCGCGGCTTTTTGT	CCCCAACGCCGCGGGTTTTTCT	39	262	VQLOW															
26,27c25
< 15979836	chr21	10701846	10701847	snp	C	G	c	G	T	433	1017	VQHIGH															
< 15979836	chr21	10701846	10701847	snp	C	T	c	G	T	433	1017	VQHIGH															
---
> 15979836	chr21	10701846	10701847	snp	C	G,T	c	G	T	433	1017	VQHIGH															
34,35c32
< 16185771	chr21	43169356	43169357	snp	C	G	c	G	T	511	534	VQHIGH	54101	NM_020639.2	NP_065690.2	RIPK4	-	CDS	3	Y	SYNONYMOUS	677	209	A	A	A	PFAM:PF00069:Pkinase
< 16185771	chr21	43169356	43169357	snp	C	T	c	G	T	511	534	VQHIGH	54101	NM_020639.2	NP_065690.2	RIPK4	-	CDS	3	Y	SYNONYMOUS	677	209	A	A	A	PFAM:PF00069:Pkinase
---
> 16185771	chr21	43169356	43169357	snp	C	G,T	c	G	T	511	534	VQHIGH	54101	NM_020639.2	NP_065690.2	RIPK4	-	CDS	3	Y	SYNONYMOUS	677	209	A	A	A	PFAM:PF00069:Pkinase
child process exited abnormally}

test cg2tsv {cg2tsv -ref} {
	exec cg cg2tsv -ref test -split 0 data/var-cgtest.tsv data/gene-cgtest.tsv tmp/temp.tsv
	catch {exec diff tmp/temp.tsv data/expected-cgtest.tsv} result
	set result
} {3c3
< #ref	test
---
> #ref	hg19
child process exited abnormally}

test cg2tsv {cg2tsv ome ins} {
	write_tab tmp/testvar.tsv {
		locus ploidy allele chromosome begin end varType reference alleleSeq totalScore hapLink xRef
		15949538 2 1 chr21 9412358 9412358 ref {} {} 178 172 VQHIGH
		15949538 2 2 chr21 9412358 9412358 ins {} A 157 172 VQHIGH
	}
	cg cg2tsv -split 0 tmp/testvar.tsv tmp/temp.tsv
	write_tab tmp/expected.tsv {
		#genomecomb 0.10
		#split 0
		#
		locus chromosome begin end type reference alt zyg alleleSeq1 alleleSeq2 totalScore1 totalScore2 xRef
		15949538 chr21 9412358 9412358 ins {} A t A {} 157 178 VQHIGH
	}
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test cg2tsv {cg2tsv various ins, snp and del} {
	write_tab tmp/testvar.tsv {
		locus ploidy allele chromosome begin end varType reference alleleSeq totalScore hapLink xRef
		15949514 2 1 chr21 9411326 9411327 ref C C 398 373 VQHIGH
		15949514 2 2 chr21 9411326 9411327 snp C G 162 208 VQHIGH dbsnp.131:rs75025155
		15949986 2 1 chr21 9425146 9425147 snp G T 434 414 VQHIGH 2722593 dbsnp.131:rs79408562
		15949986 2 2 chr21 9425146 9425147 snp G T 141 85 VQHIGH 2722594 dbsnp.131:rs79408562
		15949538 2 1 chr21 9412358 9412358 ref {} {} 178 172 VQHIGH
		15949538 2 2 chr21 9412358 9412358 ins {} A 157 172 VQHIGH
		15949540 2 1 chr21 9412441 9412444 ref ATA ATA 160 158 VQHIGH
		15949540 2 2 chr21 9412441 9412444 del ATA {} 87 117 VQHIGH
		15950476 2 1 chr21 9438355 9438355 ins {} CCCCAACGCCGCGGCTTTTTGT 39 10 VQLOW
		15950476 2 2 chr21 9438355 9438355 ins {} CCCCAACGCCGCGGGTTTTTCT 262 234 VQHIGH
		15950530 2 1 chr21 9438963 9438963 ins {} TTTTTGCCGCCGCGGCTTTTTGCCCCCGCCGCCGCGGG 89 78 VQHIGH
		15950530 2 2 chr21 9438963 9438963 ref {} {} 187 177 VQHIGH
		18824384 2 1 chr21 9553322 9553327 del GATGA {} 83
		18824384 2 2 chr21 9553322 9553327 ref GATGA GATGA 108
	}
	cg cg2tsv -split 0 tmp/testvar.tsv tmp/temp.tsv
	write_tab tmp/expected.tsv {
		#genomecomb 0.10
		#split 0
		#
		locus chromosome begin end type reference alt zyg alleleSeq1 alleleSeq2 totalScore1 totalScore2 xRef
		15949514 chr21 9411326 9411327 snp C G t G C 162 398 VQHIGH
		15949538 chr21 9412358 9412358 ins {} A t A {} 157 178 VQHIGH
		15949540 chr21 9412441 9412444 del ATA {} t {} ATA 87 160 VQHIGH
		15949986 chr21 9425146 9425147 snp G T m T T 434 141 VQHIGH
		15950476 chr21 9438355 9438355 ins {} CCCCAACGCCGCGGCTTTTTGT,CCCCAACGCCGCGGGTTTTTCT c CCCCAACGCCGCGGCTTTTTGT CCCCAACGCCGCGGGTTTTTCT 39 262 VQLOW
		15950530 chr21 9438963 9438963 ins {} TTTTTGCCGCCGCGGCTTTTTGCCCCCGCCGCCGCGGG t TTTTTGCCGCCGCGGCTTTTTGCCCCCGCCGCCGCGGG {} 89 187 VQHIGH
		18824384 chr21 9553322 9553327 del GATGA {} t {} GATGA 83 108 {}
	}
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test cg2tsv {cg2tsv various ins, snp and del -split 1} {
	write_tab tmp/testvar.tsv {
		locus ploidy allele chromosome begin end varType reference alleleSeq totalScore hapLink xRef
		15949514 2 1 chr21 9411326 9411327 ref C C 398 373 VQHIGH
		15949514 2 2 chr21 9411326 9411327 snp C G 162 208 VQHIGH dbsnp.131:rs75025155
		15949986 2 1 chr21 9425146 9425147 snp G T 434 414 VQHIGH 2722593 dbsnp.131:rs79408562
		15949986 2 2 chr21 9425146 9425147 snp G T 141 85 VQHIGH 2722594 dbsnp.131:rs79408562
		15949538 2 1 chr21 9412358 9412358 ref {} {} 178 172 VQHIGH
		15949538 2 2 chr21 9412358 9412358 ins {} A 157 172 VQHIGH
		15949540 2 1 chr21 9412441 9412444 ref ATA ATA 160 158 VQHIGH
		15949540 2 2 chr21 9412441 9412444 del ATA {} 87 117 VQHIGH
		15950476 2 1 chr21 9438355 9438355 ins {} CCCCAACGCCGCGGCTTTTTGT 39 10 VQLOW
		15950476 2 2 chr21 9438355 9438355 ins {} CCCCAACGCCGCGGGTTTTTCT 262 234 VQHIGH
		15950530 2 1 chr21 9438963 9438963 ins {} TTTTTGCCGCCGCGGCTTTTTGCCCCCGCCGCCGCGGG 89 78 VQHIGH
		15950530 2 2 chr21 9438963 9438963 ref {} {} 187 177 VQHIGH
		18824384 2 1 chr21 9553322 9553327 del GATGA {} 83
		18824384 2 2 chr21 9553322 9553327 ref GATGA GATGA 108
	}
	cg cg2tsv -split 1 tmp/testvar.tsv tmp/temp.tsv
	write_tab tmp/expected.tsv {
		#genomecomb 0.10
		#split 1
		#
		locus chromosome begin end type reference alt zyg alleleSeq1 alleleSeq2 totalScore1 totalScore2 xRef
		15949514 chr21 9411326 9411327 snp C G t G C 162 398 VQHIGH
		15949538 chr21 9412358 9412358 ins {} A t A {} 157 178 VQHIGH
		15949540 chr21 9412441 9412444 del ATA {} t {} ATA 87 160 VQHIGH
		15949986 chr21 9425146 9425147 snp G T m T T 434 141 VQHIGH
		15950476 chr21 9438355 9438355 ins {} CCCCAACGCCGCGGCTTTTTGT c CCCCAACGCCGCGGCTTTTTGT CCCCAACGCCGCGGGTTTTTCT 39 262 VQLOW
		15950476 chr21 9438355 9438355 ins {} CCCCAACGCCGCGGGTTTTTCT c CCCCAACGCCGCGGCTTTTTGT CCCCAACGCCGCGGGTTTTTCT 39 262 VQLOW
		15950530 chr21 9438963 9438963 ins {} TTTTTGCCGCCGCGGCTTTTTGCCCCCGCCGCCGCGGG t TTTTTGCCGCCGCGGCTTTTTGCCCCCGCCGCCGCGGG {} 89 187 VQHIGH
		18824384 chr21 9553322 9553327 del GATGA {} t {} GATGA 83 108 {}
	}
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test vcf2tsv {vcf2tsv} {
	exec cg vcf2tsv data/test.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test.vcf2tsv
} {}

test vcf2tsv {vcf2tsv ins and del} {
	exec cg vcf2tsv data/test2.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test2.vcf2tsv
} {}

test vcf2tsv {vcf2tsv 1000glow} {
	exec cg vcf2tsv data/test1000glow.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test1000glow.vcf2tsv
} {}

test vcf2tsv {vcf2tsv split} {
	exec cg vcf2tsv -s 1 data/test.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-tests.vcf2tsv
} {}

test vcf2tsv {vcf2tsv ins and del split} {
	exec cg vcf2tsv -s 1 data/test2.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test2s.vcf2tsv
} {}

test bed2sft {bed2sft} {
	exec cg bed2sft data/sample.bed
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

file delete tmp/temp.tsv

test sft2gff {sft2gff} {
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

test collapsealleles {collapsealleles2} {
	write_tab tmp/test.tsv {
		chromosome begin end type ref alt freq-sample
	 	chr1 4200 4200 ins {} A 0.5
	 	chr1 4200 4200 ins {} A 0.8
	 	chr1 5000 5001 snp G A 0.5
	}
	exec cg collapsealleles < tmp/test.tsv > tmp/temp.tsv
	write_tab tmp/expected.tsv {
		chromosome begin end type ref alt freq-sample
	 	chr1 4200 4200 ins {} A 0.5,0.8
	 	chr1 5000 5001 snp G A 0.5
	}
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

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
	exec cg wide tmp/long.tsv tmp/wide.tsv 2>/dev/null
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
	exec cg wide -s {sample samplea} -f id tmp/long.tsv tmp/wide.tsv 2>/dev/null
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
	exec cg wide -s {sample samplea} -f id tmp/long.tsv tmp/wide.tsv 2>/dev/null
	write_tab tmp/expected.tsv {
		id data-s1-sa1 data-s1-sa2 data-s2-sa1
	 	1 i1d11 i1d12 i1d21
	 	2 i2d11 ? i2d21 
	} {#type	wide tsv
#samplefields	sample samplea}
	exec diff tmp/wide.tsv tmp/expected.tsv
} {}

test gene2reg {gene2reg basic} {
	cg select -q {$name2 == "HES4"} data/gene_test.tsv tmp/genetmp.tsv
	cg gene2reg tmp/genetmp.tsv tmp/result.tsv
	exec diff tmp/result.tsv data/expected-gene2reg.tsv
} {}

testsummarize
