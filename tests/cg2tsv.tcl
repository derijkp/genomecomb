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

test cg2tsv {cg2tsv various ins with -split 1 sort order} {
	write_tab tmp/testvar.tsv {
		locus ploidy allele chromosome begin end varType reference alleleSeq totalScore hapLink xRef
		15949538 2 2 chr21 9412358 9412358 ins {} AA 157 172 VQHIGH
		15949538 2 2 chr21 9412358 9412358 ins {} A 178 172 VQHIGH
	}
	cg cg2tsv -split 1 tmp/testvar.tsv tmp/temp.tsv
	write_tab tmp/expected.tsv {
		#filetype	tsv/varfile
		#fileversion	0.10.0
		#split 1
		#
		locus chromosome begin end type reference alt zyg alleleSeq1 alleleSeq2 totalScore1 totalScore2 xRef
		15949538 chr21 9412358 9412358 ins {} A c AA A 157 178 VQHIGH
		15949538 chr21 9412358 9412358 ins {} AA c AA A 157 178 VQHIGH
	}
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

testsummarize

