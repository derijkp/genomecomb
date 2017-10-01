#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test select_group {group simple} {
	cg select -g type data/expected_near-vars1-reg_annot.sft
} {type	count
del	1
ins	1
snp	12}

test select_group {group query} {
	cg select -q {$coverage-sample1 > 2} -g type data/expected_near-vars1-reg_annot.sft
} {type	count
del	1
ins	1
snp	8}

test select_group {group filter} {
	cg select -g {type snp} data/expected_near-vars1-reg_annot.sft
} {type	count
snp	12}

test select_group {group filter sample} {
	write_tab tmp/temp.tsv {
		sample	val
		s1	1
		s2	1
		s1	2
	}
	cg select -g {sample s1} tmp/temp.tsv
} {sample	count
s1	2}

test select_group {group query -sr} {
	cg select -q {$coverage-sample1 > 2} -g type -sr count data/expected_near-vars1-reg_annot.sft
} {type	count
snp	8
del	1
ins	1}

test select_group {group 2 fields} {
	cg select -g {type {} alt {}} data/expected-annotate-vars_annottest-gene_test.tsv
} {type	alt	count
del		4
ins	GG	1
snp	A	11
snp	A,T	1
snp	C	6
snp	C,T	1
snp	G	6
snp	T	15
snp	T,C,A	1
sub	GGG	1}

test select_group {group 2 fields, one calculated} {
	cg select -g {type {} {altx="${alt}x"} {}} data/expected-annotate-vars_annottest-gene_test.tsv
} {type	altx	count
del	x	4
ins	GGx	1
snp	A,Tx	1
snp	Ax	11
snp	C,Tx	1
snp	Cx	6
snp	Gx	6
snp	T,C,Ax	1
snp	Tx	15
sub	GGGx	1}

test select_group {group calc col, multiple cols} {
	cg select -f {{coding=if($test_impact regexp "CDS","coding","noncoding")} *} \
	-g test_gene \
	-gc {type {snp del} coding {} count,min(begin),max(begin),max(end),avg(begin)} \
	data/expected-annotate-vars_annottest-gene_test.tsv tmp/result.tsv
	exec diff data/expected-selectgroup.tsv tmp/result.tsv
} {}

test select_group {group calc col in groupcol, multiple cols} {
	cg select -g test_gene \
	-gc {type {snp del} {coding=if($test_impact regexp "CDS","coding","noncoding")} {} count,min(begin),max(begin),max(end),avg(begin)} \
	data/expected-annotate-vars_annottest-gene_test.tsv tmp/result.tsv
	exec diff data/expected-selectgroup.tsv tmp/result.tsv
} {}

test select_group {group calc col used in group, multiple cols} {
	cg select -f {{coding=if($test_impact regexp "CDS","coding","noncoding")} *} \
	-g {coding} \
	-gc {type {snp del} count,min(begin),max(begin),max(end),avg(begin)} \
	data/expected-annotate-vars_annottest-gene_test.tsv tmp/result.tsv
	exec diff data/expected-selectgroupcoding.tsv tmp/result.tsv
} {}

test select_group {group calc col defined in group, multiple cols} {
	cg select -g {{coding=if($test_impact regexp "CDS","coding","noncoding")}} \
	-gc {type {snp del} count,min(begin),max(begin),max(end),avg(begin)} \
	data/expected-annotate-vars_annottest-gene_test.tsv tmp/result.tsv
	exec diff data/expected-selectgroupcoding.tsv tmp/result.tsv
} {}

test select_group {group calc col in func} {
	cg select -f {{size=$end - $begin}} -g type -gc {count,distinct(size)} data/expected_near-vars1-reg_annot.sft
} {type	count	distinct_size
del	1	10
ins	1	0
snp	12	1}

test select_group {group calc col in func extra space in def} {
	cg select -f {{size =$end - $begin}} -g type -gc {count,distinct(size)} data/expected_near-vars1-reg_annot.sft
} {type	count	distinct_size
del	1	10
ins	1	0
snp	12	1}

test select_group {group distinct} {
	cg select -g type -gc {sample {} count,distinct(coverage)} data/expected_near-vars1-reg_annot.sft
} {type	count-sample1	distinct_coverage-sample1	count-sample2	distinct_coverage-sample2
del	1	32	1	41
ins	1	32	1	41
snp	12	1,47,54	12	0,35,52}

test select_group {group distinct sample with list} {
	cg select -g type -gc {sample {sample2} count,distinct(coverage)} data/expected_near-vars1-reg_annot.sft
} {type	count-sample2	distinct_coverage-sample2
del	1	41
ins	1	41
snp	12	0,35,52}

test select_group {group ucount} {
	cg select -g type -gc {sample {} count,ucount(coverage)} data/expected_near-vars1-reg_annot.sft
} {type	count-sample1	ucount_coverage-sample1	count-sample2	ucount_coverage-sample2
del	1	1	1	1
ins	1	1	1	1
snp	12	3	12	3}

test select_group {group list} {
	cg select -g type -gc {sample {} list(coverage)} data/expected_near-vars1-reg_annot.sft
} {type	list_coverage-sample1	list_coverage-sample2
del	32	41
ins	32	41
snp	1,1,47,54,1,1,47,54,54,54,54,54	0,0,35,52,0,0,35,52,52,52,52,52}

test select_group {group sample sum} {
	cg select -g chromosome -gc {sequenced v sample {} sum(coverage)} data/expected_near-vars1-reg_annot.sft
} {chromosome	sum_coverage-v-sample1	sum_coverage-v-sample2
chr1	135	35
chr2	351	139}

test select_group {group sample percent} {
	cg select -g chromosome -gc {sample {} sequenced v percent} data/expected_near-vars1-reg_annot.sft
} {chromosome	percent-sample1-v	percent-sample2-v
chr1	35.71	37.50
chr2	64.29	62.50}

test select_group {group sample gpercent} {
	cg select -g chromosome -gc {sample {} sequenced v gpercent} data/expected_near-vars1-reg_annot.sft
} {chromosome	gpercent-sample1-v	gpercent-sample2-v
chr1	62.50	37.50
chr2	64.29	35.71}

test select_group {group sample and query} {
	cg select -q {$coverage-sample1 > 2} -g chromosome -gc {sample {} sequenced v count} data/expected_near-vars1-reg_annot.sft
} {chromosome	count-sample1-v	count-sample2-v
chr1	3	1
chr2	7	3}

test select_group {group sample with non-sample var} {
	cg select -g chromosome -gc {sample {} sequenced v type {} count} data/expected_near-vars1-reg_annot.sft
} {chromosome	count-sample1-v-del	count-sample1-v-ins	count-sample1-v-snp	count-sample2-v-snp
chr1	1	0	4	3
chr2	0	1	8	5}

test select_group {group where -g has sample} {
	cg select -g coverage -gc {sample {} count} data/expected_near-vars1-reg_annot.sft
} {coverage	count-sample1	count-sample2
0	0	4
1	4	0
32	2	0
35	0	2
41	0	2
47	2	0
52	0	6
54	6	0}

test select_group {group where -g has sample and non-sample field} {
	cg select -g {coverage {} type {}} -gc {sample {} count} data/expected_near-vars1-reg_annot.sft
} {coverage	type	count-sample1	count-sample2
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

test select_group {group query use calculated column in query} {
	cg select -f {cov=$coverage-sample1} -q {$cov > 2} -g type data/expected_near-vars1-reg_annot.sft
} {type	count
del	1
ins	1
snp	8}

test select_group {group query field with and without sample preference} {
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

test select_group {group query missing field in one sample} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	val-s1	test-s2
		i1	0	1
		i2	1	0
		i3	1	1
		i4	0	0
		i5	0	1
	}
	cg select -g {sample {} val 1} -gc count tmp/temp.tsv
} {sample	val	count
s1	1	2}

test select_group {group with wildcard calc col} {
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

test select_group {group with sample in -gc} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2
		1	A	B
		2	B	B
	}
	exec cg select -g {type} -gc {sample {} count} tmp/temp.tsv
} {type	count-sample1	count-sample2
A	1	0
B	1	2}

test select_group {group with sample in -gc where not all samples have the -g field} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2 other-sample3
		1	A	B	a
		2	B	B	b
	}
	exec cg select -g {type} -gc {sample {} count} tmp/temp.tsv
} {type	count-sample1	count-sample2
A	1	0
B	1	2}

test select_group {group with sample in -gc where not all samples have the -g field and limit sample} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2 other-sample3
		1	A	B	a
		2	B	B	b
	}
	exec cg select -g {type} -gc {sample {sample2} count} tmp/temp.tsv
} {type	count-sample2
B	2}

test select_group {group with wildcard calc col and sampleinfo} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2
		1	A	B
		2	B	B
	}
	write_tab tmp/temp.sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	f
	}
	exec cg select -f {{typex-*="${type-*}${gender-*}"}} -g {sample {} typex {}} -gc {count} tmp/temp.tsv
} {sample	typex	count
sample1	Am	1
sample1	Bm	1
sample2	Bf	2}

test select_group {group with wildcard calc col and sampleinfo, ignoring sample prefix} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-gatk-crsbwa-sample1	type-gatk-crsbwa-sample2
		1	A	B
		2	B	B
	}
	write_tab tmp/temp.sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	f
	}
	exec cg select -f {{typex-*="${type-*}${gender-*}"}} -g {sample {} typex {}} -gc {count} tmp/temp.tsv
} {sample	typex	count
gatk-crsbwa-sample1	Am	1
gatk-crsbwa-sample1	Bm	1
gatk-crsbwa-sample2	Bf	2}

test select_group {sampleinfo in group} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2	type-sample3
		1	A	B	A
		2	B	B	B
	}
	write_tab tmp/temp.sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select -g {type {} gender {}} -gc {sample {} count} tmp/temp.tsv
} {type	gender	count-sample1	count-sample2	count-sample3
A	f	0	0	1
A	m	1	0	0
B	f	0	2	1
B	m	1	0	0}

test select_group {sampleinfo in group, filter} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2	type-sample3
		1	A	B	A
		2	B	B	B
	}
	write_tab tmp/temp.sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select -g {type {} gender f} -gc {sample {} count} tmp/temp.tsv
} {type	gender	count-sample1	count-sample2	count-sample3
A	f	0	0	1
B	f	0	2	1}

test select_group {sampleinfo in gc} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2	type-sample3
		1	A	B	A
		2	B	B	B
	}
	write_tab tmp/temp.sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select -g {type {}} -gc {sample {} gender f count} tmp/temp.tsv
} {type	count-sample1-f	count-sample2-f	count-sample3-f
A	0	0	1
B	0	2	1}

test select_group {sampleinfo in gc, sample in g} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2	type-sample3
		1	A	B	A
		2	B	B	B
	}
	write_tab tmp/temp.sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select -g {sample {}} -gc {type {} gender {} count} tmp/temp.tsv
} {sample	count-A-f	count-A-m	count-B-f	count-B-m
sample1	0	1	0	1
sample2	0	0	2	0
sample3	1	0	1	0}

test select_group {sampleinfo in agregate, sample in g} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2	type-sample3
		1	A	B	A
		2	B	B	B
	}
	write_tab tmp/temp.sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select -g {sample {}} -gc {type {} list(gender)} tmp/temp.tsv
} {sample	list_gender-A	list_gender-B
sample1	m	m
sample2		f,f
sample3	f	f}

test select_group {sampleinfo in agregate, sample in g} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2	type-sample3
		1	A	B	A
		2	B	B	B
	}
	write_tab tmp/temp.sampleinfo.tsv {
		id	num
		sample1	1
		sample2	2
		sample3	3
	}
	exec cg select -g {sample {}} -gc {type {} list(num),sum(num),avg(num)} tmp/temp.tsv
} {sample	list_num-A	sum_num-A	avg_num-A	list_num-B	sum_num-B	avg_num-B
sample1	1	1	1.0	1	1	1.0
sample2				2,2	4	2.0
sample3	3	3	3.0	3	3	3.0}

test select_group {hidden sample} {
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

test select_group {hidden sample, sampleinfo in gc} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2	type-sample3
		1	A	B	A
		2	B	B	B
	}
	write_tab tmp/temp.sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select -g {type {}} -gc {gender {} count} tmp/temp.tsv
} {type	count-f	count-m
A	1	1
B	3	1}

test select_group {hidden sample in sampleinfo} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2	type-sample3
		1	A	B	A
		2	B	B	B
	}
	write_tab tmp/temp.sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select -g {gender {}} -gc {count} tmp/temp.tsv
} {gender	count
f	4
m	2}

test select_group {hidden sample in sampleinfo, list} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2	type-sample3
		1	A	B	A
		2	B	B	B
	}
	write_tab tmp/temp.sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select -g {gender {}} -gc {list(type)} tmp/temp.tsv
} {gender	list_type
f	B,A,B,B
m	A,B}

test select_group {hidden sample in wildcard calc cols} {
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

test select_group {hidden sample in wildcard calc cols in aggregate} {
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

test select_group {sampleinfo in code of query} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type	freq-sample1	freq-sample2	freq-sample3
		1	A	0.4	0.8	1.0
		2	B	0.8	0.9	0.3
		3	A	0.9	0.9	0.8
		4	B	0.8	0.9	0.7
	}
	write_tab tmp/temp.sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select -q {scount($freq > 0.5 and $gender eq "f") > 1} -g type tmp/temp.tsv
} {type	count
A	2
B	1}

test select_group {wildcard calc column used in gc} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type	val-sample1	val-sample2	val-sample3
		1	A	2	0	1
		2	B	5	0	2
		3	A	0	1	3
		4	B	1	9	4
	}
	exec cg select -f {{valb-*=if($val-* == 0,0,1)}} -g type -gc {valb {} count} tmp/temp.tsv
} {type	count-0	count-1
A	2	4
B	1	5}

test select_group "long format -g with sampleinfo" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		sample id	freq
		sample1	1	0.4
		sample2	1	0.8
		sample3	1	1.0
	}
	write_tab tmp/temp.tsv.sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select -q {$gender eq "f"} -g all tmp/temp.tsv
} {all	count
all	2}

test select_group "long format -g with sampleinfo in calc" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		sample id	freq
		sample1	1	0.4
		sample2	1	0.8
		sample3	1	1.0
	}
	write_tab tmp/temp.tsv.sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select -f {{t=concat($sample,$gender)}} -g t tmp/temp.tsv
} {t	count
sample1m	1
sample2f	1
sample3f	1}

test select_group "long format -g with sampleinfo" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		sample id	freq
		sample1	1	0.4
		sample2	1	0.8
		sample3	1	1.0
	}
	write_tab tmp/temp.tsv.sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select -g gender tmp/temp.tsv
} {gender	count
f	2
m	1}

test select_group "long format -gc with sampleinfo" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		sample id
		sample1	1
		sample2	1
		sample3	1
	}
	write_tab tmp/temp.tsv.sampleinfo.tsv {
		id	gender	freq
		sample1	m	0.4
		sample2	f	0.8
		sample3	f	1.0
	}
	exec cg select -g all -gc {gender {} count} tmp/temp.tsv
} {all	count-f	count-m
all	2	1}

test select_group "long format -g and -gc aggregate with sampleinfo" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		sample id
		sample1	1
		sample2	1
		sample3	1
	}
	write_tab tmp/temp.tsv.sampleinfo.tsv {
		id	gender	freq
		sample1	m	0.4
		sample2	f	0.8
		sample3	f	1.0
	}
	exec cg select -g gender -gc {list(freq)} tmp/temp.tsv
} {gender	list_freq
f	0.8,1.0
m	0.4}

test select_group "long format -g and -gc with sampleinfo" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		sample exp
		sample1	1
		sample2	1
		sample3	2
		sample4	2
	}
	write_tab tmp/temp.tsv.sampleinfo.tsv {
		id	gender	freq
		sample1	m	0.4
		sample2	f	0.8
		sample3	f	1.0
		sample4	m	0.5
	}
	exec cg select -g gender -gc {exp {} avg(freq)} tmp/temp.tsv
} {gender	avg_freq-1	avg_freq-2
f	0.8	1.0
m	0.4	0.5}

test select_group "long format -g and -gc and filter with sampleinfo" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		sample exp
		sample1	1
		sample2	1
		sample3	2
		sample4	2
	}
	write_tab tmp/temp.tsv.sampleinfo.tsv {
		id	gender	freq
		sample1	m	0.4
		sample2	f	0.8
		sample3	f	1.0
		sample4	m	0.5
	}
	exec cg select -g {gender f} -gc {exp {} avg(freq)} tmp/temp.tsv
} {gender	avg_freq-1	avg_freq-2
f	0.8	1.0}

test select_group "long format -g with sample with sampleinfo" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		sample exp
		sample1	1
		sample2	1
		sample3	2
		sample4	2
	}
	write_tab tmp/temp.tsv.sampleinfo.tsv {
		id	gender	freq
		sample1	m	0.4
		sample2	f	0.8
		sample3	f	1.0
		sample4	m	0.5
	}
	exec cg select -g {sample {} freq {0.*}} tmp/temp.tsv
} {sample	freq	count
sample1	0.4	1
sample2	0.8	1
sample4	0.5	1}

test select_group "long format -gc with sample with sampleinfo" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		sample exp
		sample1	1
		sample2	1
		sample3	2
		sample4	2
	}
	write_tab tmp/temp.tsv.sampleinfo.tsv {
		id	gender	freq
		sample1	m	0.4
		sample2	f	0.8
		sample3	f	1.0
		sample4	m	0.5
	}
	exec cg select -g gender -gc {sample {} list(freq)} tmp/temp.tsv
} {gender	list_freq-sample1	list_freq-sample2	list_freq-sample3	list_freq-sample4
f		0.8	1.0	
m	0.4			0.5}

test select_group "median" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		name v-sample1 v-sample2
		n1	1	10
		n2	2	20
		n3	3	30
		n4	4	40
		n4	5	50
	}
	exec cg select -g sample -gc {q1(v),median(v),q3(v)} tmp/temp.tsv
} {sample	q1_v	median_v	q3_v
sample1	1.5	3	4.5
sample2	15.0	30	45.0}

test select_group "decompose lists" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		name value
		n1	1
		n1,n2	2
		n2	3
	}
	exec cg select -g name tmp/temp.tsv
} {name	count
n1	1
n1,n2	1
n2	1}

test select_group "decompose lists" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		name value
		n1	1
		n1,n2	2
		n2	3
	}
	exec cg select -g +name tmp/temp.tsv
} {name	count
n1	2
n2	2}

test select_group "decompose lists with duplicates (+)" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		name value
		n1	1
		n1,n2,n2	2
		n2	3
	}
	exec cg select -g +name tmp/temp.tsv
} {name	count
n1	2
n2	3}

test select_group "decompose lists with duplicates (-)" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		name value
		n1	1
		n1,n2,n2	2
		n2	3
	}
	exec cg select -g -name tmp/temp.tsv
} {name	count
n1	2
n2	2}

test select_group "decompose lists filter" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		name value
		n1	1
		n1,n2	2
		n2	3
	}
	exec cg select -g {+name n1} tmp/temp.tsv
} {name	count
n1	2}

test select_group "decompose lists filter with duplicates (+)" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		name value
		n1	1
		n1,n2,n2	2
		n2	3
	}
	exec cg select -g {+name n2} tmp/temp.tsv
} {name	count
n2	3}

test select_group "decompose lists with duplicates (-)" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		name value
		n1	1
		n1,n2,n2	2
		n2	3
	}
	exec cg select -g {-name n2} tmp/temp.tsv
} {name	count
n2	2}

test select_group "decompose lists multiple" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		name test	value
		n1	t1,t2	1
		n1,n2	t1,t2	2
		n2	t2	3
	}
	exec cg select -g {+name * +test} tmp/temp.tsv
} {name	test	count
n1	t1	2
n1	t2	2
n2	t1	1
n2	t2	2}

test select_group "decompose lists groupcol" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		name value
		n1	1
		n1,n2	2
		n2	3
	}
	exec cg select -stack 1 -g all -gc {+name * count} tmp/temp.tsv
} {all	count-n1	count-n2
all	2	2}

test select_group "decompose lists groupcol duplicates +" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		name value
		n1	1
		n1,n2,n2	2
		n2	3
	}
	exec cg select -stack 1 -g all -gc {+name * count} tmp/temp.tsv
} {all	count-n1	count-n2
all	2	3}

test select_group "decompose lists groupcol duplicates -" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		name value
		n1	1
		n1,n2,n2	2
		n2	3
	}
	exec cg select -stack 1 -g all -gc {-name * count} tmp/temp.tsv
} {all	count-n1	count-n2
all	2	2}

test select_group "decompose lists group and groupcol" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		name test	value
		n1	t1,t2	1
		n1,n2	t1,t2	2
		n2	t2	3
	}
	exec cg select -g +name -gc {+test * count} tmp/temp.tsv
} {name	count-t1	count-t2
n1	2	2
n2	1	2}

test select_group "decompose lists group and groupcol duplicates +" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		name test	value
		n1	t1,t2	1
		n1,n1,n2	t1,t2,t2	2
		n2	t2	3
	}
	exec cg select -g +name -gc {+test * count} tmp/temp.tsv
} {name	count-t1	count-t2
n1	3	5
n2	1	3}

test select_group "decompose lists group and groupcol duplicates -" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		name test	value
		n1	t1,t2	1
		n1,n1,n2	t1,t2,t2	2
		n2	t2	3
	}
	exec cg select -g -name -gc {-test * count} tmp/temp.tsv
} {name	count-t1	count-t2
n1	2	2
n2	1	2}

test select_group "decompose lists multiple sample" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		name test-sample1	test-sample2	value
		n1	t1,t2	t1	1
		n1,n2	t1,t2	t1,t2	2
		n2	t2	t2	3
	}
	exec cg select -g {sample * +name * +test} tmp/temp.tsv
} {sample	name	test	count
sample1	n1	t1	2
sample1	n1	t2	2
sample1	n2	t1	1
sample1	n2	t2	2
sample2	n1	t1	2
sample2	n1	t2	1
sample2	n2	t1	1
sample2	n2	t2	2}

test select_group "decompose lists multiple sample duplicates -" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		name test-sample1	test-sample2	value
		n1	t1,t2	t1	1
		n1,n1,n2	t1,t2,t2	t1,t2	2
		n2	t2	t2	3
	}
	exec cg select -g {sample * -name * -test} tmp/temp.tsv
} {sample	name	test	count
sample1	n1	t1	2
sample1	n1	t2	2
sample1	n2	t1	1
sample1	n2	t2	2
sample2	n1	t1	2
sample2	n1	t2	1
sample2	n2	t1	1
sample2	n2	t2	2}

test select_group "decompose lists group and groupcol and sample" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		name test-sample1	test-sample2	value
		n1	t1,t2	t1	1
		n1,n2	t1,t2	t1,t2	2
		n2	t2	t2	3
	}
	exec cg select -g {sample * +name} -gc {+test * count} tmp/temp.tsv
} {sample	name	count-t1	count-t2
sample1	n1	2	2
sample1	n2	1	2
sample2	n1	2	1
sample2	n2	1	2}

test select_group {decompose lists sampleinfo} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	seq-sample1	seq-sample2
		1	v	r
		2	v	v
	}
	write_tab tmp/temp.sampleinfo.tsv {
		id	test
		sample1	t1,t2
		sample2	t2
	}
	exec cg select -g {seq v +test} tmp/temp.tsv
} {seq	test	count
v	t1	2
v	t2	3}

test select_group {decompose lists sampleinfo} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	seq-sample1	seq-sample2
		1	v	r
		2	v	v
	}
	write_tab tmp/temp.sampleinfo.tsv {
		id	test
		sample1	t1,t2
		sample2	t2
	}
	exec cg select -g {sample *} -gc {seq v +test * count} tmp/temp.tsv
} {sample	count-v-t1	count-v-t2
sample1	2	2
sample2	0	1}

test select_group {decompose lists sampleinfo duplicates +} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	seq-sample1	seq-sample2
		1	v	r
		2	v	v
	}
	write_tab tmp/temp.sampleinfo.tsv {
		id	test
		sample1	t1,t2,t2
		sample2	t2
	}
	exec cg select -g {sample *} -gc {seq v +test * count} tmp/temp.tsv
} {sample	count-v-t1	count-v-t2
sample1	2	4
sample2	0	1}

test select_group {decompose lists sampleinfo duplicates -} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	seq-sample1	seq-sample2
		1	v	r
		2	v	v
	}
	write_tab tmp/temp.sampleinfo.tsv {
		id	test
		sample1	t1,t2,t2
		sample2	t2
	}
	exec cg select -g {sample *} -gc {seq v -test * count} tmp/temp.tsv
} {sample	count-v-t1	count-v-t2
sample1	2	2
sample2	0	1}

test select_group {group simple empty} {
	write_tab tmp/temp.tsv {
		chromosome	begin		end	type
	}
	cg select -g all tmp/temp.tsv
} {all	count
all	0}

test select_group {group simple empty} {
	write_tab tmp/temp.tsv {
		chromosome	begin		end	type
	}
	cg select -g all -gc {type {ins del} count} tmp/temp.tsv
} {all	count-del	count-ins
all	0	0}

test select_group {group with sample in -gc, field in -g missing for some samples} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2	other-sample3
		1	A	B	X
		2	B	B	Y
	}
	exec cg select -g {type} -gc {sample {} count} tmp/temp.tsv
} {type	count-sample1	count-sample2
A	1	0
B	1	2}

test select_group {group with sample in -gc, field in -gc missing for some samples} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2	other-sample3
		1	A	B	X
		2	B	B	Y
	}
	exec cg select -g all -gc {type {} sample {} count} tmp/temp.tsv
} {all	count-A-sample1	count-B-sample1	count-B-sample2
all	1	1	2}

test select_group {group with sample in -gc, field in -gc missing for some samples, filter} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	type-sample1	type-sample2	other-sample3
		1	A	B	X
		2	B	B	Y
	}
	exec cg select -g all -gc {type {A B} sample {} count} tmp/temp.tsv
} {all	count-A-sample1	count-A-sample2	count-B-sample1	count-B-sample2
all	1	0	1	2}

testsummarize
