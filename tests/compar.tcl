#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl
package require genomecomb

test nat_compare {chr1 1} {
	nat_compare chr10 chr2
} 1

test nat_compare {chr1 1} {
	nat_compare chr1 1
} 50

test loc_compare {2 2} {
	loc_compare 2 2
} 0

test loc_compare {1 10} {
	loc_compare 1 10
} -1

test loc_compare {2 10} {
	loc_compare 2 10
} -1

test loc_compare {10 2} {
	loc_compare 10 2
} 1

test loc_compare {chr1 chr2} {
	loc_compare chr1 chr2
} -1

test loc_compare {chr1 1} {
	loc_compare chr1 1
} 0

test loc_compare {chr2 chr1} {
	loc_compare chr2 Chr1
} 1

test loc_compare {chr2 chr1} {
	loc_compare chr-2 1
} 1

test loc_compare {chr2 chr1} {
	loc_compare ab2 ab10
} -1

test loc_compare {{chr2 10} {chr2 2}} {
	loc_compare {chr2 10} {chr2 2}
} 1

test multicompar {basic} {
	test_cleantmp
	cg multicompar tmp/temp.sft data/var_annot.sft data/var_annot2.sft
	exec diff tmp/temp.sft data/expected-multicompar-var_annotvar_annot2.sft
} {} 

test multicompar {basic split} {
	test_cleantmp
	cg splitalleles data/var_annot.sft > tmp/var-sample1.tsv
	cg splitalleles data/var_annot2.sft > tmp/var-sample2.tsv
	cg splitalleles data/var_annot2seq.sft > tmp/var-sample3.tsv
	cg multicompar -split 1 tmp/temp.sft tmp/var-sample1.tsv tmp/var-sample2.tsv tmp/var-sample3.tsv
	exec diff tmp/temp.sft data/expected-multicompar-split.sft
} {} 

test multicompar {basic reannot} {
	test_cleantmp
	file copy data/var_annot.sft tmp/var-annot1.tsv
	file copy data/var_annot2.sft tmp/var-annot2.tsv
	file copy data/sreg-annot1.sft tmp/sreg-annot1.tsv
	file copy data/sreg-annot2.sft tmp/sreg-annot2.tsv
	cg multicompar tmp/temp.tsv tmp/var-annot1.tsv tmp/var-annot2.tsv
	cg multicompar_reannot tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-multicompar_reannot-var_annotvar_annot2.sft
} {} 

test multicompar {basic split reannot} {
	test_cleantmp
	cg splitalleles data/var_annot.sft > tmp/var-sample1.tsv
	cg splitalleles data/var_annot2.sft > tmp/var-sample2.tsv
	cg splitalleles data/var_annot2seq.sft > tmp/var-sample3.tsv
	file copy data/sreg-annot1.sft tmp/sreg-sample1.tsv
	file copy data/sreg-annot2.sft tmp/sreg-sample2.tsv
	file copy data/sreg-annot2.sft tmp/sreg-sample3.tsv
	cg multicompar -split 1 tmp/temp.sft tmp/var-sample1.tsv tmp/var-sample2.tsv tmp/var-sample3.tsv
	cg multicompar_reannot tmp/temp.sft
	exec diff tmp/temp.sft data/expected-multicompar-split-reannot.sft
} {} 

test multicompar {basic, sequenced already present} {
	test_cleantmp
	cg multicompar tmp/temp.sft data/var_annot.sft data/var_annot2seq.sft
	catch {exec diff tmp/temp.sft data/expected-multicompar-var_annotvar_annot2.sft} e
	set e
} {1c1
< chromosome	begin	end	type	ref	alt	sequenced-var_annot	alleleSeq1-var_annot	alleleSeq2-var_annot	name-var_annot	freq-var_annot	sequenced-var_annot2seq	alleleSeq1-var_annot2seq	alleleSeq2-var_annot2seq	name-var_annot2seq	freq-var_annot2seq
---
> chromosome	begin	end	type	ref	alt	sequenced-var_annot	alleleSeq1-var_annot	alleleSeq2-var_annot	name-var_annot	freq-var_annot	sequenced-var_annot2	alleleSeq1-var_annot2	alleleSeq2-var_annot2	name-var_annot2	freq-var_annot2
3c3
< 1	4001	4002	snp	A	C	v	A	C	test2	0.2	r	A	C	test2	0.2
---
> 1	4001	4002	snp	A	C	v	A	C	test2	0.2	v	A	C	test2	0.2
5c5
< 1	5000	5010	del	AGCGTGGCAA		v	AGCGTGGCAA		test4	0.4	r	AGCGTGGCAA		test4	0.4
---
> 1	5000	5010	del	AGCGTGGCAA		v	AGCGTGGCAA		test4	0.4	v	AGCGTGGCAA		test4	0.4
child process exited abnormally} 

test multicompar {noalt} {
	test_cleantmp
	cg multicompar tmp/temp.sft data/var_annotnoalt.sft data/var_annot2noalt.sft
	exec diff tmp/temp.sft data/expected-multicompar-var_annotvar_annot2noalt.sft
} {} 

test svmulticompar {basic} {
	test_cleantmp
	cg svmulticompar tmp/temp.sft data/cgsv1.tsv data/cgsv2.tsv data/cgsv3.tsv
	exec diff tmp/temp.sft data/expected-svmulticompar.tsv
} {} 

test makepvt {basic} {
	cg makepvt -sorted 0 data/expected-multireg-reg1-reg4.sft tmp/temp.sft
	file_read tmp/temp.sft
} {reg1	reg4	numbases	numlines
0	1	195	7
1	0	300	8
1	1	1160	9
}

test makepvt {fields} {
	cg makepvt -sorted 0 data/expected-multireg-reg1-reg4.sft tmp/temp.sft {chromosome reg1}
	file_read tmp/temp.sft
} {chromosome	reg1	numbases	numlines
1	0	10	2
1	1	20	4
2	0	170	3
2	1	130	6
3	1	200	2
M	0	5	1
M	1	10	1
X	0	10	1
X	1	100	1
Y	1	1000	3
} 

test multicompar {var and mapper naming convention} {
	test_cleantmp
	catch {file link tmp/var-varcaller1-mapper1-sample1.sft ../data/var_annot.sft}
	catch {file link tmp/var-varcaller2-mapper2-sample1.sft ../data/var_annot.sft}
	catch {file link tmp/var-varcaller1-mapper1-sample2.sft ../data/var_annot2.sft}
	catch {file link tmp/var-varcaller2-mapper2-sample2.sft ../data/var_annot2.sft}
	cg multicompar tmp/temp.sft tmp/var-varcaller1-mapper1-sample1.sft tmp/var-varcaller2-mapper2-sample1.sft tmp/var-varcaller1-mapper1-sample2.sft tmp/var-varcaller2-mapper2-sample2.sft
	set fields {chromosome	begin	end	type	ref	alt}
	foreach {from to} {
		var_annot varcaller1-mapper1-sample1 var_annot varcaller2-mapper2-sample1
		var_annot2 varcaller1-mapper1-sample2 var_annot2 varcaller2-mapper2-sample2
	} {
		foreach field {sequenced alleleSeq1 alleleSeq2 name freq} {
			lappend fields $field-$to=\$$field-$from
		}
	}
	cg select -f $fields data/expected-multicompar-var_annotvar_annot2.sft tmp/expected.sft
	exec diff tmp/temp.sft tmp/expected.sft
} {} 

test multicompar {basic reannot var-all} {
	test_cleantmp
	file copy data/var_annot.sft tmp/var-annot1.tsv
	file copy data/var_annot2.sft tmp/var-annot2.tsv
	file copy data/sreg-annot1.sft tmp/sreg-annot1.tsv
	file copy data/sreg-annot2.sft tmp/sreg-annot2.tsv
	catch {file delete temp.tsv}
	# make tmp/varall-annot1.tsv
	cg select -f {* sequenced="v"} data/var_annot.sft tmp/temp.tsv
	set f [open tmp/temp.tsv a]
	puts $f [join {chr2 4009 4010 snp C C C C teste 0.0 r} \t]
	puts $f [join {chr2 4010 4011 snp A G G C test7e 0.1 r} \t]
	close $f
	cg select -s - tmp/temp.tsv tmp/varall-annot1.tsv
	# make tmp/varall-annot2.tsv
	catch {file delete temp.tsv}
	cg select -f {* sequenced="v"} data/var_annot2.sft tmp/temp.tsv
	set f [open tmp/temp.tsv a]
	puts $f [join {chr1 4050 4060 snp N T T T test3e 0.3 r} \t]
	puts $f [join {chr2 4001 4002 snp A C G C test8e 0.2 r} \t]
	close $f
	cg select -s - tmp/temp.tsv tmp/varall-annot2.tsv
	file delete tmp/temp.tsv
	cg multicompar tmp/temp.tsv tmp/var-annot1.tsv tmp/var-annot2.tsv
	cg multicompar_reannot tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-multicompar_reannot_varall-var_annotvar_annot2.sft
} {} 

test multicompar {merge} {
	test_cleantmp
	file copy data/vars1.sft tmp/mcompar.tsv
	set c [file_read data/vars3.sft]
	regsub -all sample1 $c sample3 c
	regsub -all sample2 $c sample4 c
	file_write tmp/merge.tsv $c
	cg multicompar -listfields list tmp/mcompar.tsv tmp/merge.tsv
	exec diff tmp/mcompar.tsv data/expected-multicompar-merge.tsv
} {}

test multicompar {sort empty bug split} {
	test_cleantmp
	# this gave an incorrectly sorted file
	cg multicompar -split 1 tmp/temp.tsv data/var-compartest1.tsv data/var-compartest2.tsv
	cg checksort tmp/temp.tsv
} {}

test multicompar {error on split files without split option} {
	test_cleantmp
	# this gave an incorrectly sorted file
	cg multicompar tmp/temp.tsv data/var-compartest1.tsv data/var-compartest2.tsv
	cg checksort tmp/temp.tsv
} {*error in "*var-compartest2.tsv": file uses split alleles ("1 207806142 207806170 sub" occurs more than once and you are not running multicompar with the -split option)*} match error

test multicompar {error on badly sorted files} {
	test_cleantmp
	# this gave an incorrectly sorted file
	cg multicompar tmp/temp.tsv data/vars_sorterror1.sft data/vars_sorterror2.sft
} {*sorting error in "*vars_sorterror1.sft": "10 43198434 43198435 snp" comes before "3 52847042 52847060 del"*} match error

test_cleantmp

set ::env(PATH) $keeppath

testsummarize
