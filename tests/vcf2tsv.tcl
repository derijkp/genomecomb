#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

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

test vcf2tsv {vcf2tsv ins and del split} {
	exec cg vcf2tsv -typelist . -s 1 data/test2.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test2s.vcf2tsv
} {27,28c27,28
< 20	1110695	1110696	snp	A	G	rs6040355	67	PASS	G	T	c	1	1,2	21	6	23,27	T	G	c	1	2,1	2	0	18,2	T	T	o	0	2;2	35	4		2	10	0.333,0.667	T	1	
< 20	1110695	1110696	snp	A	T	rs6040355	67	PASS	G	T	c	1	2,1	21	6	23,27	T	G	c	1	1,2	2	0	18,2	T	T	m	0	1;1	35	4		2	10	0.333,0.667	T	1	
---
> 20	1110695	1110696	snp	A	G	rs6040355	67	PASS	G	T	c	1	1,2	21	6	23,27	T	G	c	1	2,1	2	0	18,2	T	T	o	0	2;2	35	4		2	10	0.333	T	1	
> 20	1110695	1110696	snp	A	T	rs6040355	67	PASS	G	T	c	1	2,1	21	6	23,27	T	G	c	1	1,2	2	0	18,2	T	T	m	0	1;1	35	4		2	10	0.667	T	1	
child process exited abnormally} error

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

test vcf2tsv {cg_vcf2tsv info handling} {
	exec head -8 data/test2.vcf > tmp/temp.vcf
	exec tail -n+11 data/test2.vcf >> tmp/temp.vcf
	set o [open tmp/msgs.txt w]
	cg vcf2tsv -s 1 tmp/temp.vcf tmp/temp.tsv 2>@ $o
	close $o
	file_read tmp/msgs.txt
} {line 16: info field DB not described in header, skipping
line 18: info field AA not described in header, skipping
line 19: info field AA not described in header, skipping
line 19: info field DB not described in header, skipping
line 19: info field AA not described in header, skipping
line 19: info field DB not described in header, skipping
line 20: info field AA not described in header, skipping
line 21: info field AA not described in header, skipping
line 21: info field AA not described in header, skipping
line 22: info field AA not described in header, skipping
line 23: info field AA not described in header, skipping
}

testsummarize
