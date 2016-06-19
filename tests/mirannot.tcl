#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test mir_annot {basic mir annotation} {
	test_cleantmp
	exec cg annotate -dbdir /complgen/refseq/hg19 data/vars_mirna.tsv tmp/annot_test.tsv data/mir_small.tsv
	exec diff tmp/annot_test.tsv data/expected-vars_mirna_annot.tsv
} {} 

test mir_annot {-upstream} {
	test_cleantmp
	exec cg annotate --upstreamsize 100 -dbdir /complgen/refseq/hg19 data/vars_mirna.tsv tmp/annot_test.tsv data/mir_small.tsv
	exec diff tmp/annot_test.tsv data/expected-vars_mirna_annot.tsv
} {2c2
< chr1	567603	567604	snp	G	C	pre		
---
> chr1	567603	567604	snp	G	C	pre	mir1+:upstream(a-101e);mir1-a:downstream(a+e101);mir1-b:downstream(a+e101)	mir1+;mir1-;mir1-
16c16
< chr1	567893	567894	snp	A	C	post		
---
> chr1	567893	567894	snp	A	C	post	mir1+:downstream(a+e101);mir1-a:upstream(a-101e);mir1-b:upstream(a-101e)	mir1+;mir1-;mir1-
*} error match

set ::env(PATH) $keeppath

testsummarize
