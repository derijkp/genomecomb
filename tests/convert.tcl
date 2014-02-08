#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test select {cg2tsv} {
	exec cg cg2tsv -split 0 data/var-cgtest.tsv data/gene-cgtest.tsv tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-cgtest.tsv
} {}

test select {cg2tsv without genefile} {
	exec cg cg2tsv -split 0 data/var-cgtest.tsv tmp/temp.tsv
	exec cg select -f {locus chromosome begin end type reference alt zyg alleleSeq1 alleleSeq2 totalScore1 totalScore2 xRef} data/expected-cgtest.tsv tmp/expected.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test select {cg2tsv (var2annot)} {
	exec cg var2annot -split 0 data/var-cgtest.tsv data/gene-cgtest.tsv tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-cgtest.tsv
} {}

test select {cg2tsv -split 1} {
	exec cg cg2tsv -split 1 data/var-cgtest.tsv data/gene-cgtest.tsv tmp/temp.tsv
	catch {exec diff tmp/temp.tsv data/expected-cgtest.tsv} result
	set result
} {2c2
< #split	1
---
> #split	0
21,22c21
< 15979836	chr21	10701846	10701847	snp	C	G	c	G	T	433	1017	VQHIGH															
< 15979836	chr21	10701846	10701847	snp	C	T	c	G	T	433	1017	VQHIGH															
---
> 15979836	chr21	10701846	10701847	snp	C	G,T	c	G	T	433	1017	VQHIGH															
29,30c28
< 16185771	chr21	43169356	43169357	snp	C	G	c	G	T	511	534	VQHIGH	54101	NM_020639.2	NP_065690.2	RIPK4	-	CDS	3	Y	SYNONYMOUS	677	209	A	A	A	PFAM:PF00069:Pkinase
< 16185771	chr21	43169356	43169357	snp	C	T	c	G	T	511	534	VQHIGH	54101	NM_020639.2	NP_065690.2	RIPK4	-	CDS	3	Y	SYNONYMOUS	677	209	A	A	A	PFAM:PF00069:Pkinase
---
> 16185771	chr21	43169356	43169357	snp	C	G,T	c	G	T	511	534	VQHIGH	54101	NM_020639.2	NP_065690.2	RIPK4	-	CDS	3	Y	SYNONYMOUS	677	209	A	A	A	PFAM:PF00069:Pkinase
child process exited abnormally}

test select {cg2tsv -ref} {
	exec cg cg2tsv -ref test -split 0 data/var-cgtest.tsv data/gene-cgtest.tsv tmp/temp.tsv
	catch {exec diff tmp/temp.tsv data/expected-cgtest.tsv} result
	set result
} {3c3
< #ref	test
---
> #ref	hg19
child process exited abnormally}

test select {vcf2tsv} {
	exec cg vcf2tsv data/test.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test.vcf2tsv
} {}

test select {vcf2tsv ins and del} {
	exec cg vcf2tsv data/test2.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test2.vcf2tsv
} {}

test select {vcf2tsv 1000glow} {
	exec cg vcf2tsv data/test1000glow.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test1000glow.vcf2tsv
} {}

test select {vcf2tsv split} {
	exec cg vcf2tsv -s 1 data/test.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-tests.vcf2tsv
} {}

test select {vcf2tsv ins and del split} {
	exec cg vcf2tsv -s 1 data/test2.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test2s.vcf2tsv
} {}

test select {bed2sft} {
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

test select {sft2gff} {
	exec cg vcf2tsv data/test1000glow.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test1000glow.vcf2tsv
} {}

testsummarize
