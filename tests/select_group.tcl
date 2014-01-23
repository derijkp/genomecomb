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

testsummarize
