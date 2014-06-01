#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test select {group simple} {
	cg select -g type data/expected_near-vars1-reg_annot.sft
} {type	count
del	1
ins	1
snp	12}

test select {group query} {
	cg select -q {$coverage-sample1 > 2} -g type data/expected_near-vars1-reg_annot.sft
} {type	count
del	1
ins	1
snp	8}

test select {group 2 fields} {
	cg select -g {type {} alt {}} data/expected-annotate-vars_annottest-gene_test.tsv
} {type	alt	count
del		4
ins	GG	1
snp	A	11
snp	A,T	1
snp	C	6
snp	C,T	1
snp	G	5
snp	T	15
snp	T,C,A	1
sub	GGG	1}

test select {group 2 fields, one calculated} {
	cg select -g {type {} {altx="${alt}x"} {}} data/expected-annotate-vars_annottest-gene_test.tsv
} {type	altx	count
del	x	4
ins	GGx	1
snp	A,Tx	1
snp	Ax	11
snp	C,Tx	1
snp	Cx	6
snp	Gx	5
snp	T,C,Ax	1
snp	Tx	15
sub	GGGx	1}

test select {group calc col, multiple cols} {
	cg select -f {{coding=if($test_impact regexp "CDS","coding","noncoding")} *} \
	-g test_gene \
	-gc {type {snp del} coding {} count,min(begin),max(begin),max(end),avg(begin)} \
	data/expected-annotate-vars_annottest-gene_test.tsv tmp/result.tsv
	exec diff data/expected-selectgroup.tsv tmp/result.tsv
} {}

test select {group calc col in groupcol, multiple cols} {
	cg select -g test_gene \
	-gc {type {snp del} {coding=if($test_impact regexp "CDS","coding","noncoding")} {} count,min(begin),max(begin),max(end),avg(begin)} \
	data/expected-annotate-vars_annottest-gene_test.tsv tmp/result.tsv
	exec diff data/expected-selectgroup.tsv tmp/result.tsv
} {}

test select {group calc col used in group, multiple cols} {
	cg select -f {{coding=if($test_impact regexp "CDS","coding","noncoding")} *} \
	-g {coding} \
	-gc {type {snp del} count,min(begin),max(begin),max(end),avg(begin)} \
	data/expected-annotate-vars_annottest-gene_test.tsv tmp/result.tsv
	exec diff data/expected-selectgroupcoding.tsv tmp/result.tsv
} {}

test select {group calc col defined in group, multiple cols} {
	cg select -g {{coding=if($test_impact regexp "CDS","coding","noncoding")}} \
	-gc {type {snp del} count,min(begin),max(begin),max(end),avg(begin)} \
	data/expected-annotate-vars_annottest-gene_test.tsv tmp/result.tsv
	exec diff data/expected-selectgroupcoding.tsv tmp/result.tsv
} {}

test select {group calc col in func} {
	cg select -f {{size=$end - $begin}} -g type -gc {count,distinct(size)} data/expected_near-vars1-reg_annot.sft
} {type	count	distinct_size
del	1	10
ins	1	10
snp	12	1,901}

test select {group calc col in func extra space in def} {
	cg select -f {{size =$end - $begin}} -g type -gc {count,distinct(size)} data/expected_near-vars1-reg_annot.sft
} {type	count	distinct_size
del	1	10
ins	1	10
snp	12	1,901}

test select {group distinct} {
	cg select -g type -gc {sample {} count,distinct(coverage)} data/expected_near-vars1-reg_annot.sft
} {type	sample1-count	sample1-distinct_coverage	sample2-count	sample2-distinct_coverage
del	1	32	1	41
ins	1	32	1	41
snp	12	1,47,54	12	0,35,52}

test select {group ucount} {
	cg select -g type -gc {sample {} count,ucount(coverage)} data/expected_near-vars1-reg_annot.sft
} {type	sample1-count	sample1-ucount_coverage	sample2-count	sample2-ucount_coverage
del	1	1	1	1
ins	1	1	1	1
snp	12	3	12	3}

test select {group list} {
	cg select -g type -gc {sample {} list(coverage)} data/expected_near-vars1-reg_annot.sft
} {type	sample1-list_coverage	sample2-list_coverage
del	32	41
ins	32	41
snp	1,1,47,54,1,1,47,54,54,54,54,54	0,0,35,52,0,0,35,52,52,52,52,52}

test select {group sample sum} {
	cg select -g chromosome -gc {sample {} sequenced v sum(coverage)} data/expected_near-vars1-reg_annot.sft
} {chromosome	sample1-v-sum_coverage	sample2-v-sum_coverage
chr1	135	35
chr2	351	139}

test select {group sample percent} {
	cg select -g chromosome -gc {sample {} sequenced v percent} data/expected_near-vars1-reg_annot.sft
} {chromosome	sample1-v-percent	sample2-v-percent
chr1	35.71	37.50
chr2	64.29	62.50}

test select {group sample gpercent} {
	cg select -g chromosome -gc {sample {} sequenced v gpercent} data/expected_near-vars1-reg_annot.sft
} {chromosome	sample1-v-gpercent	sample2-v-gpercent
chr1	62.50	37.50
chr2	64.29	35.71}

test select {group sample and query} {
	cg select -q {$coverage-sample1 > 2} -g chromosome -gc {sample {} sequenced v count} data/expected_near-vars1-reg_annot.sft
} {chromosome	sample1-v-count	sample2-v-count
chr1	3	1
chr2	7	3}

test select {group sample with non-sample var} {
	cg select -g chromosome -gc {sample {} sequenced v type {} count} data/expected_near-vars1-reg_annot.sft
} {chromosome	sample1-v-del-count	sample1-v-ins-count	sample1-v-snp-count	sample2-v-snp-count
chr1	1	0	4	3
chr2	0	1	8	5}

test select {group where -g has sample} {
	cg select -g coverage -gc {sample {} count} data/expected_near-vars1-reg_annot.sft
} {coverage	sample1-count	sample2-count
0	0	4
1	4	0
32	2	0
35	0	2
41	0	2
47	2	0
52	0	6
54	6	0}

test select {group where -g has sample and non-sample field} {
	cg select -g {coverage {} type {}} -gc {sample {} count} data/expected_near-vars1-reg_annot.sft
} {coverage	type	sample1-count	sample2-count
0	snp	0	4
1	snp	4	0
32	del	1	0
32	ins	1	0
35	snp	0	2
41	del	0	1
41	ins	0	1
47	snp	2	0
52	snp	0	6
54	snp	6	0}

test select {group query use calculated column in query} {
	cg select -f {cov=$coverage-sample1} -q {$cov > 2} -g type data/expected_near-vars1-reg_annot.sft
} {type	count
del	1
ins	1
snp	8}

test select {group query field with and without sample preference} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	val	val-s1	val-s2
		i1	1	0	1
		i2	1	1	0
		i3	2	1	1
		i4	0	0	0
		i5	1	0	1
	}
	cg select -q {$val == 1} -g {sample {} val 1} -gc count tmp/temp.tsv
} {sample	val	count
s1	1	1
s2	1	2}

test select {group with wildcard calc col} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2
		1	A	B
		2	B	B
	}
	exec cg select -f {{typex-*="${type-*}x"}} -g {sample {} typex {}} -gc {count} tmp/temp.tsv
} {sample	typex	count
sample1	Ax	1
sample1	Bx	1
sample2	Bx	2}

test select {group with wildcard calc col and sampledata} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2
		1	A	B
		2	B	B
	}
	write_tab tmp/temp.sampledata.tsv {
		id	gender
		sample1	m
		sample2	f
	}
	exec cg select -f {{typex-*="${type-*}${gender-*}"}} -g {sample {} typex {}} -gc {count} tmp/temp.tsv
} {sample	typex	count
sample1	Am	1
sample1	Bm	1
sample2	Bf	2}

test select {sampledata in group} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2	type-sample3
		1	A	B	A
		2	B	B	B
	}
	write_tab tmp/temp.sampledata.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select -g {type {} gender {}} -gc {sample {} count} tmp/temp.tsv
} {type	gender	sample1-count	sample2-count	sample3-count
A	f	0	0	1
A	m	1	0	0
B	f	0	2	1
B	m	1	0	0}

test select {sampledata in group, filter} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2	type-sample3
		1	A	B	A
		2	B	B	B
	}
	write_tab tmp/temp.sampledata.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select -g {type {} gender f} -gc {sample {} count} tmp/temp.tsv
} {type	gender	sample2-count	sample3-count
A	f	0	1
B	f	2	1}

test select {sampledata in gc} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2	type-sample3
		1	A	B	A
		2	B	B	B
	}
	write_tab tmp/temp.sampledata.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select -g {type {}} -gc {sample {} gender f count} tmp/temp.tsv
} {type	sample2-f-count	sample3-f-count
A	0	1
B	2	1}

test select {sampledata in gc, sample in g} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2	type-sample3
		1	A	B	A
		2	B	B	B
	}
	write_tab tmp/temp.sampledata.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select -g {sample {}} -gc {type {} gender {} count} tmp/temp.tsv
} {sample	A-f-count	A-m-count	B-f-count	B-m-count
sample1	0	1	0	1
sample2	0	0	2	0
sample3	1	0	1	0}

test select {sampledata in agregate, sample in g} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2	type-sample3
		1	A	B	A
		2	B	B	B
	}
	write_tab tmp/temp.sampledata.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select -g {sample {}} -gc {type {} list(gender)} tmp/temp.tsv
} {sample	A-list_gender	B-list_gender
sample1	m	m
sample2		f,f
sample3	f	f}

test select {sampledata in agregate, sample in g} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2	type-sample3
		1	A	B	A
		2	B	B	B
	}
	write_tab tmp/temp.sampledata.tsv {
		id	num
		sample1	1
		sample2	2
		sample3	3
	}
	exec cg select -g {sample {}} -gc {type {} list(num),sum(num),avg(num)} tmp/temp.tsv
} {sample	A-list_num	A-sum_num	A-avg_num	B-list_num	B-sum_num	B-avg_num
sample1	1	1	1.0	1	1	1.0
sample2				2,2	4	2.0
sample3	3	3	3.0	3	3	3.0}

test select {hidden sample} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2	type-sample3
		1	A	B	A
		2	B	B	B
	}
	exec cg select -g {type {}} -gc {count} tmp/temp.tsv
} {type	count
A	2
B	4}

test select {hidden sample, sampledata in gc} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2	type-sample3
		1	A	B	A
		2	B	B	B
	}
	write_tab tmp/temp.sampledata.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select -g {type {}} -gc {gender {} count} tmp/temp.tsv
} {type	f-count	m-count
A	1	1
B	3	1}

test select {hidden sample in sampledata} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2	type-sample3
		1	A	B	A
		2	B	B	B
	}
	write_tab tmp/temp.sampledata.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select -g {gender {}} -gc {count} tmp/temp.tsv
} {gender	count
f	4
m	2}

test select {hidden sample in sampledata, list} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2	type-sample3
		1	A	B	A
		2	B	B	B
	}
	write_tab tmp/temp.sampledata.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select -g {gender {}} -gc {list(type)} tmp/temp.tsv
} {gender	list_type
f	B,A,B,B
m	A,B}

test select {hidden sample in wildcard calc cols} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2	type-sample3
		1	A	B	A
		2	B	B	B
	}
	exec cg select -f {{typex-*="x$type-*"}} -g {typex {}} -gc {count} tmp/temp.tsv
} {typex	count
xA	2
xB	4}

test select {hidden sample in wildcard calc cols in aggregate} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2	type-sample3
		1	A	B	A
		2	B	B	B
	}
	exec cg select -f {{typex-*="x$type-*"}} -g {sample {}} -gc {list(typex)} tmp/temp.tsv
} {sample	list_typex
sample1	xA,xB
sample2	xB,xB
sample3	xA,xB}

testsummarize
