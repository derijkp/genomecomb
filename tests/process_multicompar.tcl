#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl
package require genomecomb

proc reorder {file dest} {
	set header [cg select -h $file]
	set samples [samples $header]
	set newheader {}
	foreach sample $samples {
		set temp [list_sub $header [list_find -glob $header *-$sample]]
		set pos1 [lsearch $temp sequenced-$sample]
		if {$pos1 != -1} {lappend newheader sequenced-$sample}
		set pos2 [lsearch $temp zyg-$sample]
		if {$pos2 != -1} {lappend newheader zyg-$sample}
		set pos3 [lsearch $temp alleleSeq1-$sample]
		if {$pos3 != -1} {lappend newheader alleleSeq1-$sample}
		set pos4 [lsearch $temp alleleSeq2-$sample]
		if {$pos4 != -1} {lappend newheader alleleSeq2-$sample}
		lappend newheader {*}[list_sub $temp -exclude [list $pos1 $pos2]]
	}
	set poss [tsv_basicfields $header 6 0]
	set poss [list_remove $poss -1]
	set header [list_concat [list_sub $header $poss] $newheader]
	cg select -f $header $file $dest
}

set testname _direct
set jobopts {}

foreach {testname jobopts} {
	_direct {}
	_d4 {-d 4}
} {

test process_multicompar$testname {process_multicompar} {
	test_cleantmp
	file mkdir tmp/samples/annot1
	file mkdir tmp/samples/annot2
	cg select -f {* zyg=zyg("")} data/var_annot.sft tmp/samples/annot1/var-annot1.tsv
	cg select -f {* zyg=zyg("")} data/var_annot2.sft tmp/samples/annot2/var-annot2.tsv
	file copy data/sreg-annot1.sft tmp/samples/annot1/sreg-annot1.tsv
	file copy data/sreg-annot2.sft tmp/samples/annot2/sreg-annot2.tsv
	cg process_multicompar {*}$::jobopts -dbdir /complgen/refseq/hg19_test -split 0 tmp
	reorder data/expected-multicompar_reannot-var_annotvar_annot2.sft tmp/expected.tsv
	exec diff tmp/compar/compar-tmp.tsv tmp/expected.tsv
	exec diff tmp/compar/sreg-tmp.tsv data/expected-sreg-multicompar.tsv
} {}

test process_multicompar$testname {process_multicompar missing sreg} {
	test_cleantmp
	file mkdir tmp/samples/annot1
	file mkdir tmp/samples/annot2
	cg select -f {* zyg=zyg("")} data/var_annot.sft tmp/samples/annot1/var-annot1.tsv
	cg select -f {* zyg=zyg("")} data/var_annot2.sft tmp/samples/annot2/var-annot2.tsv
	file copy data/sreg-annot1.sft tmp/samples/annot1/sreg-annot1.tsv
	cg process_multicompar {*}$::jobopts -dbdir /complgen/refseq/hg19_test -split 0 tmp
	reorder data/expected-multicompar_reannot-var_annotvar_annot2.sft tmp/expected.tsv
	exec diff tmp/compar/compar-tmp.tsv tmp/expected.tsv
} {4c4
< 1	4050	4060	snp	G	T	v	m	T	T	test3	0.3	?	?	?	?	?	?
---
> 1	4050	4060	snp	G	T	v	m	T	T	test3	0.3	u	u	?	?	?	?
8c8
< 2	4001	4002	snp	A	G,C	v	c	G	C	test7,test8	0.7,0.8	?	?	?	?	?	?
---
> 2	4001	4002	snp	A	G,C	v	c	G	C	test7,test8	0.7,0.8	r	r	A	A	?	?
child process exited abnormally} error

test process_multicompar$testname {process_multicompar missing sreg} {
	test_cleantmp
	file mkdir tmp/samples/annot1
	file mkdir tmp/samples/annot2
	cg select -f {* zyg=zyg("")} data/var_annot.sft tmp/samples/annot1/var-annot1.tsv
	cg select -f {* zyg=zyg("")} data/var_annot2.sft tmp/samples/annot2/var-annot2.tsv
	file copy data/sreg-annot1.sft tmp/samples/annot1/sreg-annot1.tsv
	cg process_multicompar {*}$::jobopts -skipincomplete 0 -dbdir /complgen/refseq/hg19_test -split 0 tmp
	reorder data/expected-multicompar_reannot-var_annotvar_annot2.sft tmp/expected.tsv
	exec diff tmp/compar/compar-tmp.tsv tmp/expected.tsv
} {*no sorted region file (*/sreg-annot2.tsv) or varallfile (*/varall-annot2.tsv) found: not properly processed sample*} error match

test process_multicompar$testname {process_multicompar split 3 samples} {
	test_cleantmp
	file mkdir tmp/samples/sample1
	file mkdir tmp/samples/sample2
	file mkdir tmp/samples/sample3
	cg splitalleles data/var_annot.sft > tmp/samples/sample1/var-sample1.tsv
	cg splitalleles data/var_annot2.sft > tmp/samples/sample2/var-sample2.tsv
	cg splitalleles data/var_annot2seq.sft > tmp/samples/sample3/prevar-sample3.tsv
	cg select -f {sequenced *} tmp/samples/sample3/prevar-sample3.tsv tmp/samples/sample3/var-sample3.tsv
	file copy data/sreg-annot1.sft tmp/samples/sample1/sreg-sample1.tsv
	file copy data/sreg-annot2.sft tmp/samples/sample2/sreg-sample2.tsv
	file copy data/sreg-annot2.sft tmp/samples/sample3/sreg-sample3.tsv
	cg process_multicompar {*}$::jobopts -dbdir /complgen/refseq/hg19_test -split 1 tmp
	exec diff tmp/compar/compar-tmp.tsv data/expected-multicompar-split-reannot.sft
	exec diff tmp/compar/sreg-tmp.tsv data/expected-sreg-split-reannot.sft
} {} 

test process_multicompar$testname {process_multicompar varall} {
	test_cleantmp
	file mkdir tmp/samples/annot1
	file mkdir tmp/samples/annot2
	file copy data/var_annot.sft tmp/samples/annot1/var-annot1.tsv
	file copy data/var_annot2.sft tmp/samples/annot2/var-annot2.tsv
	file copy data/sreg-annot1.sft tmp/samples/annot1/sreg-annot1.tsv
	file copy data/sreg-annot2.sft tmp/samples/annot2/sreg-annot2.tsv
	catch {file delete tmp.tsv}
	# make tmp/samples/annot1/varall-annot1.tsv
	cg select -f {* sequenced="v"} data/var_annot.sft tmp/temp.tsv
	set f [open tmp/temp.tsv a]
	puts $f [join {chr2 4009 4010 snp C C C C teste 0.0 r} \t]
	puts $f [join {chr2 4010 4011 snp A G G C test7e 0.1 u} \t]
	close $f
	cg select -s - tmp/temp.tsv tmp/samples/annot1/varall-annot1.tsv
	mklink data/expected-multicompar_reannot_varall-var_annotvar_annot2.sft tmp/expected.tsv
	# make tmp/samples/annot2/varall-annot2.tsv
	cg select -f {* sequenced="v"} data/var_annot2.sft tmp/temp.tsv
	set f [open tmp/temp.tsv a]
	puts $f [join {chr1 4050 4060 snp G G G G test3e 0.3 r} \t]
	close $f
	cg select -s - tmp/temp.tsv tmp/samples/annot2/varall-annot2.tsv
	mklink data/expected-pmulticompar_reannot_varall-var_annotvar_annot2.tsv tmp/expected.tsv
	catch {file delete tmp/temp.tsv}
	cg process_multicompar {*}$::jobopts -dbdir /complgen/refseq/hg19_test -split 0 tmp
	exec diff tmp/compar/compar-tmp.tsv tmp/expected.tsv
	exec diff tmp/compar/sreg-tmp.tsv data/expected-sreg-2sample.sft
} {} 

test process_multicompar$testname {process_multicompar -varfiles} {
	test_cleantmp
	file mkdir tmp/samples/sample1
	file mkdir tmp/samples/sample2
	file mkdir tmp/samples/sample3
	cg splitalleles data/var_annot.sft > tmp/samples/sample1/var-sample1.tsv
	cg splitalleles data/var_annot2.sft > tmp/samples/sample2/var-sample2.tsv
	cg splitalleles data/var_annot2seq.sft > tmp/samples/sample3/prevar-sample3.tsv
	cg select -f {sequenced *} tmp/samples/sample3/prevar-sample3.tsv tmp/samples/sample3/var-sample3.tsv
	file copy data/sreg-annot1.sft tmp/samples/sample1/sreg-sample1.tsv
	file copy data/sreg-annot2.sft tmp/samples/sample2/sreg-sample2.tsv
	file copy data/sreg-annot2.sft tmp/samples/sample3/sreg-sample3.tsv
	cg process_multicompar {*}$::jobopts -dbdir /complgen/refseq/hg19_test -split 1 -varfiles {tmp/samples/sample1/var-sample1.tsv tmp/samples/sample2/var-sample2.tsv} tmp
	cg select -rf {*-sample3} data/expected-multicompar-split-reannot.sft tmp/expected.tsv.temp
	cg select -q {scount($sequenced eq "v") > 0} tmp/expected.tsv.temp tmp/expected.tsv
	exec diff tmp/compar/compar-tmp.tsv tmp/expected.tsv
} {} 

test process_multicompar$testname {process_multicompar -varfiles pattern} {
	test_cleantmp
	file mkdir tmp/samples/sample1
	file mkdir tmp/samples/sample2
	file mkdir tmp/samples/s3
	cg splitalleles data/var_annot.sft > tmp/samples/sample1/var-sample1.tsv
	cg splitalleles data/var_annot2.sft > tmp/samples/sample2/var-sample2.tsv
	cg splitalleles data/var_annot2seq.sft > tmp/samples/s3/prevar-s3.tsv
	cg select -f {sequenced *} tmp/samples/s3/prevar-s3.tsv tmp/samples/s3/var-s3.tsv
	file copy data/sreg-annot1.sft tmp/samples/sample1/sreg-sample1.tsv
	file copy data/sreg-annot2.sft tmp/samples/sample2/sreg-sample2.tsv
	file copy data/sreg-annot2.sft tmp/samples/s3/sreg-sample3.tsv
	cg process_multicompar {*}$::jobopts -dbdir /complgen/refseq/hg19_test -split 1 -varfiles {tmp/samples/sample*/var-sample*.tsv} tmp
	cg select -rf {*-sample3} data/expected-multicompar-split-reannot.sft tmp/expected.tsv.temp
	cg select -q {scount($sequenced eq "v") > 0} tmp/expected.tsv.temp tmp/expected.tsv
	exec diff tmp/compar/compar-tmp.tsv tmp/expected.tsv
} {} 

}

test_cleantmp

set ::env(PATH) $keeppath

testsummarize
