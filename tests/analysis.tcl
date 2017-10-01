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

test exportplink {del} {
	test_cleantmp
	set c [file_read data/vars3.sft]
	append c "chr4\t4000\t4001\tdel\tG\t\t0.5\t\t\tv\t0.4\t\tG\tv\n"
	file_write tmp/vars.tsv $c
	exec cg exportplink tmp/vars.tsv tmp/temp
	exec diff tmp/temp.tped data/expected-vars3.tped
} {9d8
< 4	4-4000-4001-del-G-	0.004000	4000	-	-	-	G
child process exited abnormally} error

test exportplink {codegeno} {
	test_cleantmp
	exec cg exportplink -c 1 data/vars3.sft tmp/temp
	exec diff tmp/temp.tfam.pre data/expected-vars3.tfam.pre
	exec diff tmp/temp.tped data/expected-vars3.codedtped
} {}

test exportplink {samples} {
	test_cleantmp
	file_write tmp/h.tsv [join {f1 f2 f3 f4 f5 f6 f7 f8} \t]
	cg select -hf tmp/h.tsv  data/expected-vars3.tped tmp/temp.tsv
	cg select -f {f1 f2 f3 f4 f7 f8} -sh /dev/null tmp/temp.tsv tmp/expected.tsv
	exec cg exportplink -s sample2 data/vars3.sft tmp/temp
	exec diff tmp/temp.tped tmp/expected.tsv
	exec diff tmp/temp.tfam.pre data/expected-vars3.tfam.pre
} {0a1
> fam	sample1	0	0	0	-9
child process exited abnormally} error

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

test exportplink {query} {
	test_cleantmp
	exec cg exportplink -q {$chromosome == "chr1"} data/vars3.sft tmp/temp
	exec diff tmp/temp.tped data/expected-vars3.tped
} {6a7,8
> 2	2-4001-4002-snp-A-G	0.004001	4001	G	G	G	G
> 2	2-4099-5000-snp-C-T	0.004099	4099	T	T	T	T
child process exited abnormally} error

test exportplink {query error} {
	test_cleantmp
	exec cg exportplink -q {$bla == "chr1"} data/vars3.sft tmp/temp
	exec diff tmp/temp.tped data/expected-vars3.tped
} {error querying file: field bla not present in file (or sampleinfo)*} error match

file delete -force tmp/temp.tsv tmp/temp.tsv.old

testsummarize
