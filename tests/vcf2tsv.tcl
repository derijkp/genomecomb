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

testsummarize
