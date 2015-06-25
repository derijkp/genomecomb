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
} {3c3
< #split	1
---
> #split	0
22,23c22
< 15950476	chr21	9438355	9438355	ins		CCCCAACGCCGCGGCTTTTTGT	c	CCCCAACGCCGCGGCTTTTTGT	CCCCAACGCCGCGGGTTTTTCT	39	262	VQLOW															
< 15950476	chr21	9438355	9438355	ins		CCCCAACGCCGCGGGTTTTTCT	c	CCCCAACGCCGCGGCTTTTTGT	CCCCAACGCCGCGGGTTTTTCT	39	262	VQLOW															
---
> 15950476	chr21	9438355	9438355	ins		CCCCAACGCCGCGGCTTTTTGT,CCCCAACGCCGCGGGTTTTTCT	c	CCCCAACGCCGCGGCTTTTTGT	CCCCAACGCCGCGGGTTTTTCT	39	262	VQLOW															
27,28c26
< 15979836	chr21	10701846	10701847	snp	C	G	c	G	T	433	1017	VQHIGH															
< 15979836	chr21	10701846	10701847	snp	C	T	c	G	T	433	1017	VQHIGH															
---
> 15979836	chr21	10701846	10701847	snp	C	G,T	c	G	T	433	1017	VQHIGH															
35,36c33
< 16185771	chr21	43169356	43169357	snp	C	G	c	G	T	511	534	VQHIGH	54101	NM_020639.2	NP_065690.2	RIPK4	-	CDS	3	Y	SYNONYMOUS	677	209	A	A	A	PFAM:PF00069:Pkinase
< 16185771	chr21	43169356	43169357	snp	C	T	c	G	T	511	534	VQHIGH	54101	NM_020639.2	NP_065690.2	RIPK4	-	CDS	3	Y	SYNONYMOUS	677	209	A	A	A	PFAM:PF00069:Pkinase
---
> 16185771	chr21	43169356	43169357	snp	C	G,T	c	G	T	511	534	VQHIGH	54101	NM_020639.2	NP_065690.2	RIPK4	-	CDS	3	Y	SYNONYMOUS	677	209	A	A	A	PFAM:PF00069:Pkinase
child process exited abnormally}

test cg2tsv {cg2tsv -ref} {
	exec cg cg2tsv -ref test -split 0 data/var-cgtest.tsv data/gene-cgtest.tsv tmp/temp.tsv
	catch {exec diff tmp/temp.tsv data/expected-cgtest.tsv} result
	set result
} {4c4
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
		#filetype	tsv/varfile
		#fileversion	0.10.0
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
		#filetype	tsv/varfile
		#fileversion	0.10.0
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
		#filetype	tsv/varfile
		#fileversion	0.10.0
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

test vcf2tsv {vcf2tsv vars_mirna.vcf} {
	exec cg vcf2tsv -s 1 data/vars_mirna.vcf tmp/temp.tsv
	cg select -rc 1 -rf {name quality filter totalcoverage	allelecount	totalallelecount} tmp/temp.tsv tmp/temp2.tsv
	cg select -rc 1 -rf {name} data/vars_mirna.tsv tmp/expected.tsv
	exec diff tmp/temp2.tsv tmp/expected.tsv
} {22,23c22,23
< chr1	1102505	1102508	del	NNN	
< chr1	1102520	1102537	del	NNNNNNNNNNNNNNNNN	
---
> chr1	1102505	1102508	del	3	
> chr1	1102520	1102537	del	17	
child process exited abnormally} error

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

test_cleantmp

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

test liftover {basic no correctvariants} {
	test_cleantmp
	cg select -rf {test} data/expected-var_lift-hg18tohg19.tsv tmp/expected.tsv
	exec cg liftover -correctvariants 0 data/var_lift.tsv tmp/temp.tsv /complgen/refseq/liftover/hg18ToHg19.over.chain
	cg select -rf {test} tmp/temp.tsv tmp/temp2.tsv
	exec diff tmp/temp.tsv.unmapped data/expected-var_lift-hg18tohg19.tsv.unmapped
	exec diff tmp/temp2.tsv tmp/expected.tsv
} {12,13c12,13
< 1	2492275	2492276	snp	C	A	B2	r	o	19	G	G	v	t	18	A	C	1	2482141	2482142	C
< 1	2492275	2492276	snp	C	G	B1	v	m	18	G	G	r	o	18	A	C	1	2482141	2482142	C
---
> 1	2492275	2492276	snp	G	C	B1	v	m	18	C	C	r	o	18	T	G	1	2482141	2482142	C
> 1	2492275	2492276	snp	G	T	B2	r	o	19	C	C	v	t	18	T	G	1	2482141	2482142	C
} error regexp

test liftover {basic liftover with correctvariants} {
	test_cleantmp
	exec cg liftover data/var_lift.tsv tmp/temp.tsv /complgen/refseq/liftover/hg18ToHg19.over.chain
	exec diff tmp/temp.tsv data/expected-var_lift-hg18tohg19.tsv
	exec diff tmp/temp.tsv.unmapped data/expected-var_lift-hg18tohg19.tsv.unmapped
} {}

test liftregion {basic} {
	test_cleantmp
	exec cg liftregion data/reg_lift.tsv tmp/temp.tsv /complgen/refseq/liftover/hg18ToHg19.over.chain
	exec diff tmp/temp.tsv data/expected-reg_lift-hg18tohg19.tsv
	exec diff tmp/temp.tsv.unmapped data/expected-reg_lift-hg18tohg19.tsv.unmapped
} {}

test liftregion {changed refseq} {
	test_cleantmp
	# 1: alt becomes ref
	# 2,3: two alleles (one becomes ref)
	# 4: complement with refchange to complement!
	# 5: complement with refchange to other base
	write_tab tmp/temp.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	score-sample1	alleleSeq1-sample1	alleleSeq2-sample1	sequenced-sample2	zyg-sample2	score-sample2	alleleSeq1-sample2	alleleSeq2-sample2
		1	1609629	1609630	snp	T	C	v	m	1	C	C	r	r	2	T	T
		1	2474628	2474629	snp	G	C	v	t	2	G	C	r	o	4	T	T
		1	2474628	2474629	snp	G	T	r	o	3	G	C	v	m	6	T	T
		1	2489102	2489103	snp	A	T	v	m	4	T	T	r	r	5	A	A
		2	89594172	89594173	snp	T	A	v	t	5	T	A	r	r	8	T	T
		2	110743803	110743804	snp	A	T	v	t	6	A	T	r	r	10	A	A
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	score-sample1	alleleSeq1-sample1	alleleSeq2-sample1	sequenced-sample2	zyg-sample2	score-sample2	alleleSeq1-sample2	alleleSeq2-sample2	hg18_chromosome	hg18_begin	hg18_end	hg18_ref
		1	1619766	1619767	snp	C	T	r	r	1	C	C	v	m	2	T	T	1	1609629	1609630	T
		1	2484768	2484769	snp	T	C	v	c	2	G	C	r	r	4	T	T	1	2474628	2474629	G
		1	2484768	2484769	snp	T	G	v	c	2	G	C	r	r	4	T	T	1	2474628	2474629	G
		1	2499242	2499243	snp	C	A	r	o	4	T	T	v	m	5	A	A	1	2489102	2489103	A
		1	2499242	2499243	snp	C	T	v	m	4	T	T	r	o	5	A	A	1	2489102	2489103	A
		2	89563941	89563942	snp	T	A	v	t	5	A	T	v	m	8	A	A	2	89594172	89594173	T
		2	110577444	110577445	snp	C	A	v	c	6	T	A	r	o	10	T	T	2	110743803	110743804	A
		2	110577444	110577445	snp	C	T	v	c	6	T	A	v	m	10	T	T	2	110743803	110743804	A
	}
	file delete tmp/lifted.tsv
	exec cg liftover tmp/temp.tsv tmp/lifted.tsv /complgen/refseq/liftover/hg18ToHg19.over.chain
	cg select -rc 1 tmp/lifted.tsv tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test liftregion {changed refseq, unsplit} {
	test_cleantmp
	# 1: alt becomes ref
	# 2,3: two alleles (one becomes ref)
	# 4: complement with refchange to complement!
	# 5: complement with refchange to other base
	write_tab tmp/temp.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	score-sample1	alleleSeq1-sample1	alleleSeq2-sample1	sequenced-sample2	zyg-sample2	score-sample2	alleleSeq1-sample2	alleleSeq2-sample2
		1	1609629	1609630	snp	T	C	v	m	1	C	C	r	r	2	T	T
		1	2474628	2474629	snp	G	C,T	v	t	2	G	C	v	m	4	T	T
		1	2489102	2489103	snp	A	T	v	m	4	T	T	r	r	5	A	A
		2	89594172	89594173	snp	T	A	v	t	5	T	A	r	r	8	T	T
		2	110743803	110743804	snp	A	T	v	t	6	A	T	r	r	10	A	A
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	score-sample1	alleleSeq1-sample1	alleleSeq2-sample1	sequenced-sample2	zyg-sample2	score-sample2	alleleSeq1-sample2	alleleSeq2-sample2	hg18_chromosome	hg18_begin	hg18_end hg18_ref
		1	1619766	1619767	snp	C	T	r	r	1	C	C	v	m	2	T	T	1	1609629	1609630	T
		1	2484768	2484769	snp	T	C,G	v	c	2	G	C	r	r	4	T	T	1	2474628	2474629	G
		1	2499242	2499243	snp	C	A,T	v	m	4	T	T	v	m	5	A	A	1	2489102	2489103	A
		2	89563941	89563942	snp	T	A	v	t	5	A	T	v	m	8	A	A	2	89594172	89594173	T
		2	110577444	110577445	snp	C	A,T	v	c	6	T	A	v	m	10	T	T	2	110743803	110743804	A
	}
	file delete tmp/lifted.tsv
	exec cg liftover -split 0 tmp/temp.tsv tmp/lifted.tsv /complgen/refseq/liftover/hg18ToHg19.over.chain
	cg select -rc 1 tmp/lifted.tsv tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test liftover {variants with different ref ending up in same spot -s 1} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		chromosome	begin	end	type	ref	alt
		2       109820576       109820577       snp     G       T
		2	110858379	110858380	snp	A	A
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	hg18_chromosome	hg18_begin	hg18_end	hg18_ref
		2       110463287       110463288       snp     G       T	2	109820576       109820577	G
	}
	exec cg liftover tmp/temp.tsv tmp/lifted.tsv /complgen/refseq/liftover/hg18ToHg19.over.chain
	exec cg correctvariants -f 1 -s 1 tmp/lifted.tsv tmp/result.tsv.temp /complgen/refseq/hg19
	cg select -rc 1 tmp/result.tsv.temp tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test liftover {variants with different alt ending up in same spot -s 0} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		chromosome	begin	end	type	ref	alt
		2       109820576       109820577       snp     G       C
		2	110858379	110858380	snp	A	T
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	hg18_chromosome	hg18_begin	hg18_end	hg18_ref
		2       110463287       110463288       snp     G       A,C,T	2	110858379	110858380	A
	}
	exec cg liftover tmp/temp.tsv tmp/lifted.tsv /complgen/refseq/liftover/hg18ToHg19.over.chain
	file delete tmp/result.tsv.temp
	exec cg correctvariants -f 1 -s 0 tmp/lifted.tsv tmp/result.tsv.temp /complgen/refseq/hg19
	cg select -rc 1 tmp/result.tsv.temp tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test liftover {changed refseq, regionsfile} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		chromosome	begin	end	type	ref	alt
		1	1610084	1610085	snp	A	G
	}
	write_tab tmp/reg.tsv {
		chromosome	begin	end
		1	609620	2474628
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	hg18_chromosome	hg18_begin	hg18_end	hg18_ref
		1	1619766	1619767	snp	C	T	1	1609629	1609630	T
		1	1620224	1620225	snp	G	A	1	1610084	1610085	A
		1	1620884	1620885	snp	T	C	1	1610744	1610745	C
	}
	file delete tmp/lifted.tsv
	exec cg liftover -regionfile tmp/reg.tsv tmp/temp.tsv tmp/lifted.tsv /complgen/refseq/liftover/hg18ToHg19.over.chain
	cg select -rc 1 tmp/lifted.tsv tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test liftover {changed refseq, regionsfile with multiple samples} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	alleleSeq1-sample1	alleleSeq2-sample1	sequenced-sample2	zyg-sample2	alleleSeq1-sample2	alleleSeq2-sample2
		1	1610084	1610085	snp	A	G	v	m	G	G	r	r	A	A
	}
	write_tab tmp/reg.tsv {
		chromosome	begin	end	sreg-sample1 sreg-sample2
		1	609620	2474628	1	0
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	alleleSeq1-sample1	alleleSeq2-sample1	sequenced-sample2	zyg-sample2	alleleSeq1-sample2	alleleSeq2-sample2	hg18_chromosome	hg18_begin	hg18_end	hg18_ref
		1	1619766	1619767	snp	C	T	v	m	T	T	u	u	?	?	1	1609629	1609630	T
		1	1620224	1620225	snp	G	A	r	r	G	G	v	m	A	A	1	1610084	1610085	A
		1	1620884	1620885	snp	T	C	v	m	C	C	u	u	?	?	1	1610744	1610745	C
	}
	file delete tmp/lifted.tsv
	exec cg liftover -regionfile tmp/reg.tsv tmp/temp.tsv tmp/lifted.tsv /complgen/refseq/liftover/hg18ToHg19.over.chain
	cg select -rc 1 tmp/lifted.tsv tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

testsummarize
