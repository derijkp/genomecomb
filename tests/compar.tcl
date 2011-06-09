#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test multicompar {basic} {
	file delete -force temp.sft
	cg multicompar temp.sft data/var_annot.sft data/var_annot2.sft
	exec diff temp.sft data/expected-multicompar-var_annotvar_annot2.sft
} {} 

file delete -force temp.sft

set ::env(PATH) $keeppath

testsummarize
