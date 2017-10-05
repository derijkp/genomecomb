#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

proc write_testdata {} {
	write_tab tmp/var-sample1.tsv {
		chromosome	begin	end	type	ref	alt	sequenced	zyg alleleSeq1	alleleSeq2	name	score
		1	10	11	snp	G	A	r	o	G	C	1A10	1
		1	10	11	snp	G	C	v	t	G	C	1C10	2
		1	10	11	snp	G	G	r	o	G	C	1G10	3
		1	11	12	snp	T	A	v	m	A	A	1A11	4
	}
	cg multicompar -split 1 tmp/multicompar1.tsv tmp/var-sample1.tsv
	write_tab tmp/var-sample2.tsv {
		chromosome	begin	end	type	ref	alt	sequenced	zyg alleleSeq1	alleleSeq2	name	score
		1	10	11	snp	G	C	v	m	C	C	2C10	5
		1	20	21	snp	G	T	v	m	T	T	2T20	6
	}
	cg multicompar -split 1 tmp/multicompar2.tsv tmp/var-sample2.tsv
	write_tab tmp/var-sample3.tsv {
		chromosome	begin	end	type	ref	alt	sequenced	zyg alleleSeq1	alleleSeq2	name score
		1	10	11	snp	G	C	v	t	G	C	3C10	7
		1	20	21	snp	G	T	v	m	T	T	3T20	8
	}
	cg multicompar -split 1 tmp/multicompar3.tsv tmp/var-sample3.tsv
	write_tab tmp/sreg-sample1.tsv {
		chromosome	begin	end
		1	0	20
	}
	file copy tmp/sreg-sample1.tsv tmp/sreg-sample2.tsv
	cg multicompar -split 1 tmp/multicompar12.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
}

test select "2 samples" {
	write_testdata
	exec cg multiselect -q {$name matches "*10"} tmp/var-sample1.tsv tmp/var-sample2.tsv > tmp/result.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1 alleleSeq1-sample1	alleleSeq2-sample1	name-sample1	score-sample1	sequenced-sample2	zyg-sample2 alleleSeq1-sample2	alleleSeq2-sample2	name-sample2	score-sample2
		1	10	11	snp	G	A	r	o	G	C	1A10	1	?	?	?	?	?	?
		1	10	11	snp	G	C	v	t	G	C	1C10	2	v	m	C	C	2C10	5
		1	10	11	snp	G	G	r	o	G	C	1G10	3	?	?	?	?	?	?
	}
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test select "2 samples -combine files" {
	write_testdata
	exec cg multiselect -q {$name matches "*10"} -combine files -o tmp/result_ tmp/var-sample1.tsv tmp/var-sample2.tsv
	cg select -q {$name matches "*10"} tmp/var-sample1.tsv tmp/expected1.tsv
	cg select -q {$name matches "*10"} tmp/var-sample2.tsv tmp/expected2.tsv
	exec diff tmp/result_var-sample1.tsv tmp/expected1.tsv
	exec diff tmp/result_var-sample2.tsv tmp/expected2.tsv
} {}

test select "2 multicompars" {
	write_testdata
	exec cg multiselect -split 1 -q {scount($name matches "*10")} tmp/multicompar1.tsv tmp/multicompar2.tsv > tmp/result.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1 alleleSeq1-sample1	alleleSeq2-sample1	name-sample1	score-sample1	sequenced-sample2	zyg-sample2 alleleSeq1-sample2	alleleSeq2-sample2	name-sample2	score-sample2
		1	10	11	snp	G	A	r	o	G	C	1A10	1	?	?	?	?	?	?
		1	10	11	snp	G	C	v	t	G	C	1C10	2	v	m	C	C	2C10	5
		1	10	11	snp	G	G	r	o	G	C	1G10	3	?	?	?	?	?	?
	}
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test select "-g" {
	write_testdata
	exec cg multiselect -combine cat -g {sample {} alt} tmp/multicompar1.tsv tmp/multicompar2.tsv
} {sample	alt	count
sample1	A	2
sample1	C	1
sample1	G	1
sample2	C	1
sample2	T	1}

test select "-samples" {
	write_testdata
	exec cg multiselect -samples {sample1 sample3} -o tmp/result.tsv tmp/multicompar12.tsv tmp/multicompar3.tsv
	cg multicompar -split 1 tmp/expected.tsv tmp/var-sample1.tsv tmp/multicompar3.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

testsummarize
