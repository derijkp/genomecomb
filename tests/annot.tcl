#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test reg_annot {basic} {
	exec cg annotate data/vars1.sft tmp/temp.sft data/reg_annot.sft 2> /dev/null
	exec cg select -rf {list} tmp/temp.sft tmp/temp2.sft
	exec diff tmp/temp2.sft data/expected-vars1-reg_annot.sft
} {} 

test reg_annot {basic, multiple fields} {
	file mkdir tmp
	cg select -f {chromosome begin end type ref alt} data/vars1.sft tmp/vars.sft
	file copy -force data/reg_annot.sft tmp/reg_annot.sft
	file_write tmp/reg_annot.sft.opt "fields\t{type begin end}\n"
	exec cg annotate tmp/vars.sft tmp/temp.sft tmp/reg_annot.sft 2> /dev/null
	exec diff tmp/temp.sft data/expected-vars1-reg_annot-multi.sft
} {} 

test reg_annot {near} {
	exec cg annotate -near 1000 data/vars1.sft tmp/temp.sft data/reg_annot.sft 2> /dev/null
	exec cg select -rf {list} tmp/temp.sft tmp/temp2.sft
	exec diff tmp/temp2.sft data/expected_near-vars1-reg_annot.sft
} {} 

test reg_annot {near indels} {
	exec cg select -q {$type == "del" || $type == "ins"} data/vars1.sft data/indels1.sft 2> /dev/null
	exec cg annotate -near 50 data/vars1.sft tmp/temp.sft data/indels1.sft 2> /dev/null
	exec cg select -rf {list} tmp/temp.sft tmp/temp2.sft
	exec diff tmp/temp2.sft data/expected-vars1-indels.sft
} {} 

test var_annot {basic} {
	exec cg annotate data/vars1.sft tmp/temp.sft data/var_annot.sft 2> /dev/null
	exec cg select -rf {list} tmp/temp.sft tmp/temp2.sft
	exec diff tmp/temp2.sft data/expected-vars1-var_annot.sft
} {} 

test var_annot {different types on same pos} {
	exec cg annotate data/vars2.tsv tmp/temp.tsv data/var_annot3.tsv 2> /dev/null
	exec diff tmp/temp.tsv data/expected-vars2-var_annot3.tsv
} {} 

test var_annot {gene} {
	exec cg annotate -dbdir /complgen/refseq/hg18 data/vars_annottest.sft tmp/temp.sft data/gene_test.tsv 2> /dev/null
	exec diff tmp/temp.sft data/expected-annotate-vars_annottest-gene_test.tsv
} {} 

test reg_annot {basic, extra comments} {
	file_write tmp/temp2.sft "# a comment\n"
	exec cat data/vars1.sft >> tmp/temp2.sft
	exec cg annotate tmp/temp2.sft tmp/temp.sft data/reg_annot.sft 2> /dev/null
	exec cg select -rf {list} tmp/temp.sft tmp/temp2.sft
	exec diff tmp/temp2.sft data/expected-vars1-reg_annot.sft
} {1d0
< # a comment	
child process exited abnormally} error

test var_annot {different types on same pos, extra comments} {
	file_write tmp/temp2.sft "# a comment\n"
	exec cat data/vars2.tsv >> tmp/temp2.sft
	exec cg annotate tmp/temp2.sft tmp/temp.tsv data/var_annot3.tsv 2> /dev/null
	exec diff tmp/temp.tsv data/expected-vars2-var_annot3.tsv
} {1d0
< # a comment	
child process exited abnormally} error

test var_annot {gene, extra comments} {
	file_write tmp/temp2.sft "# a comment\n# another comment\n"
	exec cat data/vars_annottest.sft >> tmp/temp2.sft
	exec cg annotate -dbdir /complgen/refseq/hg18 tmp/temp2.sft tmp/temp.sft data/gene_test.tsv 2> /dev/null
	exec diff tmp/temp.sft data/expected-annotate-vars_annottest-gene_test.tsv
} {1,2d0
< # a comment	
< # another comment	 
child process exited abnormally} error

file delete -force tmp/temp.sft
file delete -force tmp/temp2.sft

set ::env(PATH) $keeppath

testsummarize
