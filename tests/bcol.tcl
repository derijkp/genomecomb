#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

file delete -force temp.sft
file delete -force temp.sft.index

test bcol_index {basic} {
	file copy -force data/expected-annotate-vars_annottest-gene_test.tsv temp.sft
	exec cg size temp.sft
} {46} 

test bcol_index {basic} {
	file copy -force data/expected-annotate-vars_annottest-gene_test.tsv temp.sft
	exec cg index temp.sft
	exec cg bcol get temp.sft.index/lines.bcol 0 2
} {76 103 228} 

file delete -force temp.sft
file delete -force temp.sft.index

set ::env(PATH) $keeppath

testsummarize
