#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test exportplink {basic} {
	test_cleantmp
	exec cg exportplink data/vars3.sft tmp/temp
	exec diff tmp/temp.tfam.pre data/expected-vars3.tfam.pre
	exec diff tmp/temp.tped data/expected-vars3.tped
} {}

test exportplink {names with -} {
	catch {file delete {*}[glob tmp/temp*]}
	set header [cg select -h data/vars3.sft]
	set header [string_change $header {-sample1 -m1-sample1 -sample2 -m2-sample1}]
	cg select -nh $header data/vars3.sft tmp/tempsrc.tsv
	set c [file_read data/expected-vars3.tfam.pre]
	set c [string_change $c {sample1 m1-sample1 sample2 m2-sample1}]
	file_write tmp/expected-temp.tfam.pre $c
	exec cg exportplink tmp/tempsrc.tsv tmp/temp
	exec diff tmp/temp.tfam.pre tmp/expected-temp.tfam.pre
	exec diff tmp/temp.tped data/expected-vars3.tped
} {}

file delete -force tmp/temp.tsv tmp/temp.tsv.old

testsummarize
