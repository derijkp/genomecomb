#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test multicompar {basic} {
	file delete -force temp.sft
	cg multicompar temp.sft data/var_annot.sft data/var_annot2.sft
	exec diff temp.sft data/expected-multicompar-var_annotvar_annot2.sft
} {} 

test multicompar {basic, sequenced already present} {
	file delete -force tmp/temp.sft
	cg multicompar tmp/temp.sft data/var_annot.sft data/var_annot2seq.sft
	catch {exec diff tmp/temp.sft data/expected-multicompar-var_annotvar_annot2.sft} e
	set e
} {1c1
< chromosome	begin	end	type	ref	alt	alleleSeq1-var_annot	alleleSeq2-var_annot	name-var_annot	freq-var_annot	sequenced-var_annot	alleleSeq1-var_annot2seq	alleleSeq2-var_annot2seq	name-var_annot2seq	freq-var_annot2seq	sequenced-var_annot2seq
---
> chromosome	begin	end	type	ref	alt	alleleSeq1-var_annot	alleleSeq2-var_annot	name-var_annot	freq-var_annot	sequenced-var_annot	alleleSeq1-var_annot2	alleleSeq2-var_annot2	name-var_annot2	freq-var_annot2	sequenced-var_annot2
3c3
< 1	4001	4002	snp	A	C	A	C	test2	0.2	v	A	C	test2	0.2	r
---
> 1	4001	4002	snp	A	C	A	C	test2	0.2	v	A	C	test2	0.2	v
5c5
< 1	5000	5010	del	AGCGTGGCAA		AGCGTGGCAA		test4	0.4	v	AGCGTGGCAA		test4	0.4	r
---
> 1	5000	5010	del	AGCGTGGCAA		AGCGTGGCAA		test4	0.4	v	AGCGTGGCAA		test4	0.4	v
child process exited abnormally} 

test multicompar {noalt} {
	file delete -force temp.sft
	cg multicompar temp.sft data/var_annotnoalt.sft data/var_annot2noalt.sft
	exec diff temp.sft data/expected-multicompar-var_annotvar_annot2noalt.sft
} {} 

test makepvt {basic} {
	cg makepvt -sorted 0 data/expected-multireg-reg1-reg4.sft temp.sft
	file_read temp.sft
} {reg1	reg4	numbases	numlines
0	1	195	7
1	0	300	8
1	1	1160	9
}

test makepvt {fields} {
	cg makepvt -sorted 0 data/expected-multireg-reg1-reg4.sft temp.sft {chromosome reg1}
	file_read temp.sft
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

file delete -force temp.sft

set ::env(PATH) $keeppath

testsummarize
