#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test mir_annot {basic mir annotation} {
	test_cleantmp
	exec cg annotate -dbdir /complgen/refseq/hg19 data/vars_mirna.tsv tmp/annot_test.tsv data/mir_small.tsv
	exec diff tmp/annot_test.tsv data/expected-vars_mirna_annot.tsv
} {} 

set ::env(PATH) $keeppath

testsummarize
