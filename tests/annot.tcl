#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test reg_annot {basic} {
	exec cg annotate data/vars1.sft data/temp.sft data/reg_annot.sft 2> /dev/null
	exec diff data/temp.sft data/expected-vars1-reg_annot.sft
} {} 

test reg_annot {near} {
	exec cg annotate -near 1000 data/vars1.sft data/temp.sft data/reg_annot.sft 2> /dev/null
	exec diff data/temp.sft data/expected_near-vars1-reg_annot.sft
} {} 

test reg_annot {near indels} {
	exec cg select -q {$type == "del" || $type == "ins"} data/vars1.sft data/indels1.sft 2> /dev/null
	exec cg annotate -near 50 data/vars1.sft data/temp.sft data/indels1.sft 2> /dev/null
	exec diff data/temp.sft data/expected-vars1-indels.sft
} {} 

test var_annot {basic} {
	exec cg annotate data/vars1.sft data/temp.sft data/var_annot.sft 2> /dev/null
	exec diff data/temp.sft data/expected-vars1-var_annot.sft
} {} 

set ::env(PATH) $keeppath

testsummarize
