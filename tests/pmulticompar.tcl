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

test pmulticompar {basic} {
	test_cleantmp
	cg pmulticompar tmp/temp.sft data/var_annot.sft data/var_annot2.sft
	reorder data/expected-multicompar-var_annotvar_annot2.sft tmp/expected.tsv
	exec diff tmp/temp.sft tmp/expected.tsv
} {} 

test pmulticompar {basic with 3} {
	test_cleantmp
	cg pmulticompar tmp/temp.tsv data/var_annot.sft data/var_annot2.sft data/var_annot3.sft
	reorder data/expected-multicompar-var_annotvar_annot3.sft tmp/expected.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar {basic split} {
	test_cleantmp
	cg splitalleles data/var_annot.sft > tmp/var-sample1.tsv
	cg splitalleles data/var_annot2.sft > tmp/var-sample2.tsv
	cg splitalleles data/var_annot3.sft > tmp/var-sample3.tsv
	cg pmulticompar -split 1 tmp/temp.sft tmp/var-sample1.tsv tmp/var-sample2.tsv tmp/var-sample3.tsv
	exec diff tmp/temp.sft data/expected-multicompar-split.sft
} {} 

test pmulticompar {basic reannot} {
	test_cleantmp
	cg select -f {* zyg=zyg("")} data/var_annot.sft tmp/var-annot1.tsv
	cg select -f {* zyg=zyg("")} data/var_annot2.sft tmp/var-annot2.tsv
	file copy data/sreg-annot1.sft tmp/sreg-annot1.tsv
	file copy data/sreg-annot2.sft tmp/sreg-annot2.tsv
	cg pmulticompar tmp/temp.tsv tmp/var-annot1.tsv tmp/var-annot2.tsv
	reorder data/expected-multicompar_reannot-var_annotvar_annot2.sft tmp/expected.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar {basic split reannot} {
	test_cleantmp
	cg splitalleles data/var_annot.sft > tmp/var-sample1.tsv
	cg splitalleles data/var_annot2.sft > tmp/var-sample2.tsv
	cg splitalleles data/var_annot2seq.sft > tmp/prevar-sample3.tsv
	cg select -f {sequenced *} tmp/prevar-sample3.tsv tmp/var-sample3.tsv
	file copy data/sreg-annot1.sft tmp/sreg-sample1.tsv
	file copy data/sreg-annot2.sft tmp/sreg-sample2.tsv
	file copy data/sreg-annot2.sft tmp/sreg-sample3.tsv
	cg pmulticompar -split 1 tmp/temp.sft tmp/var-sample1.tsv tmp/var-sample2.tsv tmp/var-sample3.tsv
	exec diff tmp/temp.sft data/expected-multicompar-split-reannot.sft
} {} 

test pmulticompar {split reannot test diff alleles} {
	test_cleantmp
	write_tab tmp/var-sample1.tsv {
		chromosome	begin	end	type	ref	alt	sequenced	zyg alleleSeq1	alleleSeq2	name
		1	10	11	snp	G	A	r	o	G	C	1A10
		1	10	11	snp	G	C	v	t	G	C	1C10
		1	10	11	snp	G	G	r	o	G	C	1G10
		1	11	12	snp	T	A	v	m	A	A	1A11
	}
	write_tab tmp/var-sample2.tsv {
		chromosome	begin	end	type	ref	alt	sequenced	zyg alleleSeq1	alleleSeq2	name
		1	10	11	snp	G	C	v	m	C	C	2C10
	}
	write_tab tmp/sreg-sample1.tsv {
		chromosome	begin	end
		1	0	20
	}
	write_tab tmp/expected.sft {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1 alleleSeq1-sample1	alleleSeq2-sample1	name-sample1	sequenced-sample2	zyg-sample2 alleleSeq1-sample2	alleleSeq2-sample2	name-sample2
		1	10	11	snp	G	A	r	o	G	C	1A10	r	o	C	C	?
		1	10	11	snp	G	C	v	t	G	C	1C10	v	m	C	C	2C10
		1	10	11	snp	G	G	r	o	G	C	1G10	r	o	C	C	?
		1	11	12	snp	T	A	v	m	A	A	1A11	r	r	T	T	?
	}
	file copy tmp/sreg-sample1.tsv tmp/sreg-sample2.tsv
	cg pmulticompar -split 1 tmp/temp.sft tmp/var-sample1.tsv tmp/var-sample2.tsv
	exec diff tmp/temp.sft tmp/expected.sft
} {} 

test pmulticompar {basic, sequenced already present} {
	test_cleantmp
	cg pmulticompar tmp/temp.sft data/var_annot.sft data/var_annot2seq.sft
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

test pmulticompar {var and mapper naming convention} {
	test_cleantmp
	catch {file link tmp/var-varcaller1-mapper1-sample1.sft ../data/var_annot.sft}
	catch {file link tmp/var-varcaller2-mapper2-sample1.sft ../data/var_annot.sft}
	catch {file link tmp/var-varcaller1-mapper1-sample2.sft ../data/var_annot2.sft}
	catch {file link tmp/var-varcaller2-mapper2-sample2.sft ../data/var_annot2.sft}
	cg pmulticompar tmp/temp.sft tmp/var-varcaller1-mapper1-sample1.sft tmp/var-varcaller2-mapper2-sample1.sft tmp/var-varcaller1-mapper1-sample2.sft tmp/var-varcaller2-mapper2-sample2.sft
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

test pmulticompar {basic reannot varall} {
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
	catch {file delete tmp/temp.tsv}
	cg pmulticompar tmp/temp.tsv tmp/var-annot1.tsv tmp/var-annot2.tsv
	cg pmulticompar_reannot tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-multicompar_reannot_varall-var_annotvar_annot2.sft
} {} 

test pmulticompar {basic reannot varall split} {
	test_cleantmp
	cg splitalleles data/var_annot.sft tmp/var-annot1.tsv
	cg splitalleles data/var_annot2.sft tmp/var-annot2.tsv
	file copy data/sreg-annot1.sft tmp/sreg-annot1.tsv
	file copy data/sreg-annot2.sft tmp/sreg-annot2.tsv
	catch {file delete temp.tsv}
	# make tmp/varall-annot1.tsv
	cg select -f {* sequenced="v"} data/var_annot.sft tmp/temp1.tsv
	cg splitalleles tmp/temp1.tsv tmp/temp.tsv
	set f [open tmp/temp.tsv a]
	puts $f [join {chr2 4009 4010 snp C C C C teste 0.0 r} \t]
	puts $f [join {chr2 4010 4011 snp A G G C test7e 0.1 r} \t]
	close $f
	cg select -s - tmp/temp.tsv tmp/varall-annot1.tsv
	# make tmp/varall-annot2.tsv
	catch {file delete temp.tsv}
	cg select -f {* sequenced="v"} data/var_annot2.sft tmp/temp1.tsv
	cg splitalleles tmp/temp1.tsv tmp/temp.tsv
	set f [open tmp/temp.tsv a]
	puts $f [join {chr1 4050 4060 snp N T T T test3e 0.3 r} \t]
	puts $f [join {chr2 4001 4002 snp A C G C test8e 0.2 r} \t]
	close $f
	cg select -s - tmp/temp.tsv tmp/varall-annot2.tsv
	cg splitalleles data/expected-multicompar_reannot_varall-var_annotvar_annot2.sft tmp/expected.tsv
	file delete tmp/temp.tsv
	cg pmulticompar -split 1 tmp/temp.tsv tmp/var-annot1.tsv tmp/var-annot2.tsv
	cg pmulticompar_reannot tmp/temp.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

#test multicompar {merge} {
#	test_cleantmp
#	file copy data/vars1.sft tmp/mcompar.tsv
#	set c [file_read data/vars3.sft]
#	regsub -all sample1 $c sample3 c
#	regsub -all sample2 $c sample4 c
#	file_write tmp/merge.tsv $c
#	cg pmulticompar -listfields list tmp/mcompar.tsv tmp/merge.tsv
#	exec diff tmp/mcompar.tsv data/expected-multicompar-merge.tsv
#} {}

test pmulticompar {sort empty bug split} {
	test_cleantmp
	# this gave an incorrectly sorted file
	cg pmulticompar -split 1 tmp/temp.tsv data/var-compartest1.tsv data/var-compartest2.tsv
	cg checksort tmp/temp.tsv
} {}

test pmulticompar {error on split files without split option} {
	test_cleantmp
	# this gave an incorrectly sorted file
	cg pmulticompar tmp/temp.tsv data/var-compartest1.tsv data/var-compartest2.tsv
	cg checksort tmp/temp.tsv
} {*error in "*var-compartest2.tsv": file uses split alleles ("*1 207806142 207806170 sub" occurs more than once and you are not running pmulticompar with the -split option)*} match error

test pmulticompar {error on badly sorted files} {
	test_cleantmp
	# this gave an incorrectly sorted file
	cg pmulticompar tmp/temp.tsv data/vars_sorterror1.sft data/vars3.sft
} {File (*vars_sorterror1.sft) is not correctly sorted (sort correctly using "cg select -s -")
chr10:43198434-43198435:snp:G came before chr3:52847042-52847060:del:*} match error

test pmulticompar {error on badly sorted files 2} {
	test_cleantmp
	# this gave an incorrectly sorted file
	cg pmulticompar tmp/temp.tsv data/vars_sorterror1.sft data/vars_sorterror2.sft
} {File (*vars_sorterror2.sft) is not correctly sorted (sort correctly using "cg select -s -")
chr3:52847303-52847304:snp:G came before chr3:52847042-52847060:del:*} match error

test pmulticompar {split reannot split multiallelic varall,no sreg,zyg} {
	test_cleantmp
	write_tab tmp/var-sample1.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	score
		chr1	100	101	snp	T	C	t	T	C	40
	}
	write_tab tmp/var-sample2.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	score
		chr1	100	101	snp	T	G	m	G	G	30
	}
	file copy tmp/var-sample1.tsv tmp/varall-sample1.tsv
	file copy tmp/var-sample2.tsv tmp/varall-sample2.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	alleleSeq1-sample1	alleleSeq2-sample1	score-sample1	sequenced-sample2	zyg-sample2	alleleSeq1-sample2	alleleSeq2-sample2	score-sample2
		1	100	101	snp	T	C	v	t	T	C	40	r	o	G	G	30
		1	100	101	snp	T	G	r	o	T	C	40	v	m	G	G	30
	}
	catch {file delete tmp/temp.tsv}
	cg pmulticompar -split 1 tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	cg pmulticompar_reannot tmp/temp.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar {split reannot split multiallelic varall,sreg,zyg} {
	test_cleantmp
	write_tab tmp/var-sample1.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	score
		chr1	100	101	snp	T	C	t	T	C	40
	}
	write_tab tmp/var-sample2.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	score
		chr1	100	101	snp	T	G	m	G	G	30
	}
	write_tab tmp/sreg-sample1.tsv {
		chromosome	begin	end
		chr1	50	200
	}
	file copy tmp/var-sample1.tsv tmp/varall-sample1.tsv
	file copy tmp/var-sample2.tsv tmp/varall-sample2.tsv
	file copy tmp/sreg-sample1.tsv tmp/sreg-sample2.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	alleleSeq1-sample1	alleleSeq2-sample1	score-sample1	sequenced-sample2	zyg-sample2	alleleSeq1-sample2	alleleSeq2-sample2	score-sample2
		1	100	101	snp	T	C	v	t	T	C	40	r	o	G	G	30
		1	100	101	snp	T	G	r	o	T	C	40	v	m	G	G	30
	}
	catch {file delete tmp/temp.tsv}
	cg pmulticompar -split 1 tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	cg pmulticompar_reannot tmp/temp.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar {split reannot split multiallelic only sreg,zyg} {
	test_cleantmp
	write_tab tmp/var-sample1.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	score
		chr1	100	101	snp	T	C	t	T	C	40
	}
	write_tab tmp/var-sample2.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	score
		chr1	100	101	snp	T	G	m	G	G	30
	}
	write_tab tmp/sreg-sample1.tsv {
		chromosome	begin	end
		chr1	50	200
	}
	file copy tmp/sreg-sample1.tsv tmp/sreg-sample2.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	alleleSeq1-sample1	alleleSeq2-sample1	score-sample1	sequenced-sample2	zyg-sample2	alleleSeq1-sample2	alleleSeq2-sample2	score-sample2
		1	100	101	snp	T	C	v	t	T	C	40	r	o	G	G	?
		1	100	101	snp	T	G	r	o	T	C	?	v	m	G	G	30
	}
	catch {file delete tmp/temp.tsv}
	cg pmulticompar -split 1 tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	cg pmulticompar_reannot tmp/temp.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test pmulticompar {split reannot split multiallelic, only sreg} {
	test_cleantmp
	write_tab tmp/var-sample1.tsv {
		chromosome	begin	end	type	ref	alt	alleleSeq1	alleleSeq2
		chr1	100	101	snp	T	C	T	C
	}
	write_tab tmp/var-sample2.tsv {
		chromosome	begin	end	type	ref	alt	alleleSeq1	alleleSeq2
		chr1	100	101	snp	T	G	G	G
	}
	write_tab tmp/sreg-sample1.tsv {
		chromosome	begin	end
		chr1	50	200
	}
	file copy tmp/sreg-sample1.tsv tmp/sreg-sample2.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	alleleSeq1-sample1	alleleSeq2-sample1	sequenced-sample2	alleleSeq1-sample2	alleleSeq2-sample2
		1	100	101	snp	T	C	v	T	C	r	G	G
		1	100	101	snp	T	G	r	T	C	v	G	G
	}
	catch {file delete tmp/temp.tsv}
	cg pmulticompar -split 1 tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	cg pmulticompar_reannot tmp/temp.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar {split reannot split multiallelic varall,sreg,zyg, check ref indels} {
	test_cleantmp
	write_tab tmp/var-sample1.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	score
		chr1	100	100	ins	{}	C	t	{}	C	30
	}
	write_tab tmp/var-sample2.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	score
	}
	write_tab tmp/sreg-sample1.tsv {
		chromosome	begin	end
		chr1	50	200
	}
	write_tab tmp/varall-sample1.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	score
		chr1	100	100	ins	{}	C	t	{}	C	30
		chr1	100	101	snp	A	G	r	G	G	31
	}
	write_tab tmp/varall-sample2.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	score
		chr1	100	101	snp	A	G	r	G	G	40
	}
	file copy tmp/sreg-sample1.tsv tmp/sreg-sample2.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	alleleSeq1-sample1	alleleSeq2-sample1	score-sample1	sequenced-sample2	zyg-sample2	alleleSeq1-sample2	alleleSeq2-sample2	score-sample2
		1	100	100	ins	{}	C	v	t	{}	C	30	r	r	{}	{}	40
	}
	catch {file delete tmp/temp.tsv}
	cg pmulticompar -split 1 tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	cg pmulticompar_reannot tmp/temp.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar {split reannot split multiallelic ins,sreg,zyg} {
	test_cleantmp
	write_tab tmp/var-sample1.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	score
		chr1	100	100	ins	{}	C	t	{}	C	40
	}
	write_tab tmp/var-sample2.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	score
		chr1	100	100	ins	{}	CC	t	{}	CC	40
	}
	write_tab tmp/sreg-sample1.tsv {
		chromosome	begin	end
		chr1	50	200
	}
	write_tab tmp/varall-sample1.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	score
		chr1	100	100	ins	{}	C	t	{}	C	40
		chr1	100	101	snp	A	G	r	G	G	40
	}
	write_tab tmp/varall-sample2.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	score
		chr1	100	100	ins	{}	CC	t	{}	CC	40
		chr1	100	101	snp	A	G	r	G	G	40
	}
	file copy tmp/sreg-sample1.tsv tmp/sreg-sample2.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	alleleSeq1-sample1	alleleSeq2-sample1	score-sample1	sequenced-sample2	zyg-sample2	alleleSeq1-sample2	alleleSeq2-sample2	score-sample2
		1	100	100	ins	{}	C	v	t	{}	C	40	r	o	{}	CC	40
		1	100	100	ins	{}	CC	r	o	{}	C	40	v	t	{}	CC	40
	}
	catch {file delete tmp/temp.tsv}
	cg pmulticompar -split 1 tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	cg pmulticompar_reannot tmp/temp.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar {split reannot split multiallelic varall,sreg,zyg, check ref indels, overlapping snp} {
	test_cleantmp
	write_tab tmp/sreg-sample1.tsv {
		chromosome	begin	end
		chr1	50	200
	}
	file copy tmp/sreg-sample1.tsv tmp/sreg-sample2.tsv
	write_tab tmp/var-sample1.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	score
		chr1	100	100	ins	{}	C	t	{}	C	40
		chr1	151	152	snp	A	T	t	A	T	32
		chr1	160	162	del	NN	{}	m	{}	{}	33
	}
	write_tab tmp/var-sample2.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	score
		chr1	150	152	del	GT	{}	t	GT	{}	30
		chr1	160	161	snp	A	C	t	A	C	33
		chr1	160	161	snp	A	T	t	A	T	33
	}
	write_tab tmp/varall-sample1.tsv {
		chromosome	begin	end	type	ref	alt	sequenced	alleleSeq1	alleleSeq2	score
		chr1	100	100	ins	{}	C	t	{}	C	40
		chr1	100	101	snp	G	G	r	G	G	41
		chr1	150	151	snp	G	G	r	G	G	30
		chr1	151	152	snp	A	T	t	A	T	32
		chr1	152	153	snp	T	T	r	T	T	31
		chr1	160	162	del	NN	{}	m	{}	{}	33
	}
	write_tab tmp/varall-sample2.tsv {
		chromosome	begin	end	type	ref	alt	sequenced	alleleSeq1	alleleSeq2	score
		chr1	100	101	snp	G	G	r	G	G	41
		chr1	150	152	del	GT	{}	t	GT	{}	30
		chr1	150	151	snp	G	G	r	G	G	30
		chr1	151	152	snp	A	A	r	A	A	31
		chr1	152	153	snp	A	A	r	A	A	32
		chr1	160	161	snp	A	C	t	A	C	33
		chr1	160	161	snp	A	T	t	A	T	33
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	alleleSeq1-sample1	alleleSeq2-sample1	score-sample1	sequenced-sample2	zyg-sample2	alleleSeq1-sample2	alleleSeq2-sample2	score-sample2
		1	100	100	ins	{}	C	v	t	{}	C	40	r	r	{}	{}	41
		1	150	152	del	GT	{}	r	r	GT	GT	30	v	t	GT	{}	30
		1	151	152	snp	A	T	v	t	A	T	32	r	o	A	@	31
		1	160	161	snp	A	C	r	o	@	@	?	v	t	A	C	33
		1	160	161	snp	A	T	r	o	@	@	?	v	t	A	T	33
		1	160	162	del	NN	{}	v	m	{}	{}	33	r	r	NN	NN	33
	}
	catch {file delete tmp/temp.tsv}
	cg pmulticompar -split 1 tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	cg pmulticompar_reannot tmp/temp.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar {targets} {
	test_cleantmp
	write_tab tmp/sreg-sample1.tsv {
		chromosome	begin	end
		chr1	50	200
	}
	file copy tmp/sreg-sample1.tsv tmp/sreg-sample2.tsv
	write_tab tmp/var-sample1.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	score
		chr1	100	100	ins	{}	C	t	{}	C	40
		chr1	151	152	snp	A	T	t	A	T	32
		chr1	160	162	del	NN	{}	m	{}	{}	33
	}
	write_tab tmp/var-sample2.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	score
		chr1	150	152	del	GT	{}	t	GT	{}	30
		chr1	160	161	snp	A	C	t	A	C	33
		chr1	160	161	snp	A	T	t	A	T	33
	}
	write_tab tmp/varall-sample1.tsv {
		chromosome	begin	end	type	ref	alt	sequenced	alleleSeq1	alleleSeq2	score
		chr1	100	100	ins	{}	C	t	{}	C	40
		chr1	100	101	snp	G	G	r	G	G	41
		chr1	150	151	snp	G	G	r	G	G	30
		chr1	151	152	snp	A	T	t	A	T	32
		chr1	152	153	snp	T	T	r	T	T	31
		chr1	160	162	del	NN	{}	m	{}	{}	33
	}
	write_tab tmp/varall-sample2.tsv {
		chromosome	begin	end	type	ref	alt	sequenced	alleleSeq1	alleleSeq2	score
		chr1	100	101	snp	G	G	r	G	G	41
		chr1	150	152	del	GT	{}	t	GT	{}	30
		chr1	150	151	snp	G	G	r	G	G	30
		chr1	151	152	snp	A	A	r	A	A	31
		chr1	152	153	snp	A	A	r	A	A	32
		chr1	160	161	snp	A	C	t	A	C	33
		chr1	160	161	snp	A	T	t	A	T	33
	}
	write_tab tmp/targets-testsnps.tsv {
		chromosome	begin	end	type	ref	alt	name
		chr1	100	101	snp	G	C	snp1
		chr1	150	152	del	GT	{}	snp2
		chr1	160	161	snp	A	T	snp3
		chr1	165	166	snp	N	G	snp4
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	alleleSeq1-sample1	alleleSeq2-sample1	score-sample1	sequenced-sample2	zyg-sample2	alleleSeq1-sample2	alleleSeq2-sample2	score-sample2	testsnps
		1	100	100	ins	{}	C	v	t	{}	C	40	r	r	{}	{}	41	{}
		1	100	101	snp	G	C	r	r	G	G	41	r	r	G	G	41	snp1
		1	150	152	del	GT	{}	r	r	GT	GT	30	v	t	GT	{}	30	snp2
		1	151	152	snp	A	T	v	t	A	T	32	r	o	A	@	31	{}
		1	160	161	snp	A	C	r	o	@	@	?	v	t	A	C	33	{}
		1	160	161	snp	A	T	r	o	@	@	?	v	t	A	T	33	snp3
		1	160	162	del	NN	{}	v	m	{}	{}	33	r	r	NN	NN	33	{}
		1	165	166	snp	N	G	u	u	?	?	?	u	u	?	?	?	snp4
	}
	catch {file delete tmp/temp.tsv}
	cg pmulticompar -split 1 -targetsfile tmp/targets-testsnps.tsv tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	cg pmulticompar_reannot tmp/temp.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar {targets, no name} {
	test_cleantmp
	write_tab tmp/sreg-sample1.tsv {
		chromosome	begin	end
		chr1	50	200
	}
	file copy tmp/sreg-sample1.tsv tmp/sreg-sample2.tsv
	write_tab tmp/var-sample1.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	score
		chr1	100	100	ins	{}	C	t	{}	C	40
		chr1	151	152	snp	A	T	t	A	T	32
		chr1	160	162	del	NN	{}	m	{}	{}	33
	}
	write_tab tmp/var-sample2.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	score
		chr1	150	152	del	GT	{}	t	GT	{}	30
		chr1	160	161	snp	A	C	t	A	C	33
		chr1	160	161	snp	A	T	t	A	T	33
	}
	write_tab tmp/varall-sample1.tsv {
		chromosome	begin	end	type	ref	alt	sequenced	alleleSeq1	alleleSeq2	score
		chr1	100	100	ins	{}	C	t	{}	C	40
		chr1	100	101	snp	G	G	r	G	G	41
		chr1	150	151	snp	G	G	r	G	G	30
		chr1	151	152	snp	A	T	t	A	T	32
		chr1	152	153	snp	T	T	r	T	T	31
		chr1	160	162	del	NN	{}	m	{}	{}	33
	}
	write_tab tmp/varall-sample2.tsv {
		chromosome	begin	end	type	ref	alt	sequenced	alleleSeq1	alleleSeq2	score
		chr1	100	101	snp	G	G	r	G	G	41
		chr1	150	152	del	GT	{}	t	GT	{}	30
		chr1	150	151	snp	G	G	r	G	G	30
		chr1	151	152	snp	A	A	r	A	A	31
		chr1	152	153	snp	A	A	r	A	A	32
		chr1	160	161	snp	A	C	t	A	C	33
		chr1	160	161	snp	A	T	t	A	T	33
	}
	write_tab tmp/targets-testsnps.tsv {
		chromosome	begin	end	type	ref	alt
		chr1	100	101	snp	G	C
		chr1	150	152	del	GT	{}
		chr1	160	161	snp	A	T
		chr1	165	166	snp	N	G
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	alleleSeq1-sample1	alleleSeq2-sample1	score-sample1	sequenced-sample2	zyg-sample2	alleleSeq1-sample2	alleleSeq2-sample2	score-sample2	testsnps
		1	100	100	ins	{}	C	v	t	{}	C	40	r	r	{}	{}	41	{}
		1	100	101	snp	G	C	r	r	G	G	41	r	r	G	G	41	1
		1	150	152	del	GT	{}	r	r	GT	GT	30	v	t	GT	{}	30	1
		1	151	152	snp	A	T	v	t	A	T	32	r	o	A	@	31	{}
		1	160	161	snp	A	C	r	o	@	@	?	v	t	A	C	33	{}
		1	160	161	snp	A	T	r	o	@	@	?	v	t	A	T	33	1
		1	160	162	del	NN	{}	v	m	{}	{}	33	r	r	NN	NN	33	{}
		1	165	166	snp	N	G	u	u	?	?	?	u	u	?	?	?	1
	}
	catch {file delete tmp/temp.tsv}
	cg pmulticompar -split 1 -targetsfile tmp/targets-testsnps.tsv tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	cg pmulticompar_reannot tmp/temp.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar {targets only} {
	test_cleantmp
	write_tab tmp/sreg-sample1.tsv {
		chromosome	begin	end
		chr1	50	200
	}
	file copy tmp/sreg-sample1.tsv tmp/sreg-sample2.tsv
	write_tab tmp/varall-sample1.tsv {
		chromosome	begin	end	type	ref	alt	sequenced	alleleSeq1	alleleSeq2	score
		chr1	100	100	ins	{}	C	t	{}	C	40
		chr1	100	101	snp	G	G	r	G	G	41
		chr1	150	151	snp	G	G	r	G	G	30
		chr1	151	152	snp	A	T	t	A	T	32
		chr1	152	153	snp	T	T	r	T	T	31
		chr1	160	162	del	NN	{}	m	{}	{}	33
	}
	write_tab tmp/varall-sample2.tsv {
		chromosome	begin	end	type	ref	alt	sequenced	alleleSeq1	alleleSeq2	score
		chr1	100	101	snp	G	G	r	G	G	41
		chr1	150	152	del	GT	{}	t	GT	{}	30
		chr1	150	151	snp	G	G	r	G	G	30
		chr1	151	152	snp	A	A	r	A	A	31
		chr1	152	153	snp	A	A	r	A	A	32
		chr1	160	161	snp	A	C	t	A	C	33
		chr1	160	161	snp	A	T	t	A	T	33
	}
	write_tab tmp/targets-testsnps.tsv {
		chromosome	begin	end	type	ref	alt
		chr1	100	101	snp	G	C
		chr1	150	152	del	GT	{}
		chr1	160	161	snp	A	T
		chr1	165	166	snp	N	G
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	alleleSeq1-sample1	alleleSeq2-sample1	score-sample1	sequenced-sample2	zyg-sample2	alleleSeq1-sample2	alleleSeq2-sample2	score-sample2	testsnps
		1	100	100	ins	{}	C	v	t	{}	C	40	r	r	{}	{}	41	{}
		1	100	101	snp	G	C	r	r	G	G	41	r	r	G	G	41	1
		1	150	152	del	GT	{}	r	r	GT	GT	30	v	t	GT	{}	30	1
		1	151	152	snp	A	T	v	t	A	T	32	r	o	A	@	31	{}
		1	160	161	snp	A	C	r	o	@	@	?	v	t	A	C	33	{}
		1	160	161	snp	A	T	r	o	@	@	?	v	t	A	T	33	1
		1	160	162	del	NN	{}	v	m	{}	{}	33	r	r	NN	NN	33	{}
		1	165	166	snp	N	G	u	u	?	?	?	u	u	?	?	?	1
	}
	write_tab tmp/temp.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	alleleSeq1-sample1	alleleSeq2-sample1	score-sample1	sequenced-sample2	zyg-sample2	alleleSeq1-sample2	alleleSeq2-sample2	score-sample2
		1	100	100	ins	{}	C	v	t	{}	C	40	r	r	{}	{}	41
		1	150	152	del	GT	{}	r	r	GT	GT	30	v	t	GT	{}	30
		1	151	152	snp	A	T	v	t	A	T	32	r	o	A	@	31
		1	160	161	snp	A	C	r	o	@	@	?	v	t	A	C	33
		1	160	161	snp	A	T	r	o	@	@	?	v	t	A	T	33
		1	160	162	del	NN	{}	v	m	{}	{}	33	r	r	NN	NN	33
	}
	cg pmulticompar -split 1 -targetsfile tmp/targets-testsnps.tsv tmp/temp.tsv
	cg pmulticompar_reannot tmp/temp.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test_cleantmp

set ::env(PATH) $keeppath

testsummarize