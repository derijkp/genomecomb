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
		#genomecomb 0.9
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
		#genomecomb 0.9
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
		#genomecomb 0.9
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

testsummarize
