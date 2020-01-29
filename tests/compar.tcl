#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl
package require genomecomb

test nat_compare {chr10 chr2} {
	nat_compare chr10 chr2
} 1

test nat_compare {chr1 1} {
	nat_compare chr1 1
} 1

test nat_compare {test-100 test-1-1} {
	nat_compare test-100 test-1-1
} 1

test nat_compare {-100 -1} {
	nat_compare -100 -1
} -1

test nat_compare {{a -100} {a -1}} {
	nat_compare {a -100} {a -1}
} -1 

test nat_compare {a-100 a-1} {
	nat_compare a-100 a-1
} 1

test nat_compare {+1 -1} {
	nat_compare +1 -1
} 1

test nat_compare {a+b aa} {
	nat_compare a+b aa
} -1

test nat_compare {a +} {
	nat_compare + a
} -1

test nat_compare {a aa} {
	nat_compare a aa
} -1

test nat_compare {multiple} {
	set error {}
	foreach {a b result} {
		+	-	1
		-	+	-1
		+1	1	1
		-1	1	-1
		-1	2	-1
		-1	+1	-1
		a-1	a+1	-1
		10	+1	1
		100	+1-1	1
		100	1+1	1
		1	+1	-1
		a+1	a1	1
		1	-1	1
		+1	-1	1
		a+1	a-1	1
		a1	a-1	1
		a+2	a1	1
	} {
		set temp [nat_compare $a $b]
		if {$temp ne $result} {lappend error "$a $b  =  $temp instead of $result"}
	}
	join $error \n
} {}

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

test loc_compare {{a-2 2} {a-1 10}} {
	loc_compare {a-1 10} {a-2 2}
} -1

test loc_compare {{{a -2} 2} {{a -1} 10}} {
	loc_compare {{a -1} 10} {{a -2} 2}
} 1

test multicompar {basic} {
	test_cleantmp
	cg multicompar -split 0 tmp/temp.sft data/var_annot.sft data/var_annot2.sft
	exec diff tmp/temp.sft data/expected-multicompar-var_annotvar_annot2.sft
} {} 

test multicompar {basic with 3} {
	test_cleantmp
	foreach sample {annot annot2 annot3} {
		file copy data/var_$sample.sft tmp/var-$sample.tsv
	}
	cg multicompar -split 0 tmp/temp.tsv tmp/var-annot.tsv tmp/var-annot2.tsv tmp/var-annot3.tsv
	exec diff tmp/temp.tsv data/expected-multicompar-var_annotvar_annot3.sft
} {} 

test multicompar {basic split} {
	test_cleantmp
	cg splitalleles data/var_annot.sft > tmp/var-sample1.tsv
	cg splitalleles data/var_annot2.sft > tmp/var-sample2.tsv
	cg splitalleles data/var_annot3.sft > tmp/var-sample3.tsv
	cg multicompar -split 1 tmp/temp.sft tmp/var-sample1.tsv tmp/var-sample2.tsv tmp/var-sample3.tsv
	exec diff tmp/temp.sft data/expected-multicompar-split.sft
} {} 

test multicompar {basic reannot} {
	test_cleantmp
	cg select -f {* zyg=zyg("")} data/var_annot.sft tmp/var-annot1.tsv
	cg select -f {* zyg=zyg("")} data/var_annot2.sft tmp/var-annot2.tsv
	file copy data/sreg-annot1.sft tmp/sreg-annot1.tsv
	file copy data/sreg-annot2.sft tmp/sreg-annot2.tsv
	cg multicompar -split 0 tmp/temp.tsv tmp/var-annot1.tsv tmp/var-annot2.tsv
	cg multicompar_reannot tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-multicompar_reannot-var_annotvar_annot2.sft
} {} 

test multicompar {basic split reannot} {
	test_cleantmp
	cg splitalleles data/var_annot.sft > tmp/var-sample1.tsv
	cg splitalleles data/var_annot2.sft > tmp/var-sample2.tsv
	cg splitalleles data/var_annot2seq.sft > tmp/prevar-sample3.tsv
	cg select -f {sequenced *} tmp/prevar-sample3.tsv tmp/var-sample3.tsv
	file copy data/sreg-annot1.sft tmp/sreg-sample1.tsv
	file copy data/sreg-annot2.sft tmp/sreg-sample2.tsv
	file copy data/sreg-annot2.sft tmp/sreg-sample3.tsv
	cg multicompar -split 1 tmp/temp.sft tmp/var-sample1.tsv tmp/var-sample2.tsv tmp/var-sample3.tsv
	cg multicompar_reannot tmp/temp.sft
	exec diff tmp/temp.sft data/expected-multicompar-split-reannot.sft
} {} 

test multicompar {basic split reannot paged} {
	test_cleantmp
	cg splitalleles data/var_annot.sft > tmp/var-sample1.tsv
	cg splitalleles data/var_annot2.sft > tmp/var-sample2.tsv
	cg splitalleles data/var_annot2seq.sft > tmp/prevar-sample3.tsv
	cg select -f {sequenced *} tmp/prevar-sample3.tsv tmp/var-sample3.tsv
	file copy data/sreg-annot1.sft tmp/sreg-sample1.tsv
	file copy data/sreg-annot2.sft tmp/sreg-sample2.tsv
	file copy data/sreg-annot2.sft tmp/sreg-sample3.tsv
	cg multicompar -split 1 tmp/temp.sft tmp/var-sample1.tsv tmp/var-sample2.tsv tmp/var-sample3.tsv
	cg multicompar_reannot -paged 2 tmp/temp.sft
	exec diff tmp/temp.sft data/expected-multicompar-split-reannot.sft
} {} 

test multicompar {basic split reannot paged with pagedstart} {
	test_cleantmp
	cg splitalleles data/var_annot.sft > tmp/var-sample1.tsv
	cg splitalleles data/var_annot2.sft > tmp/var-sample2.tsv
	cg splitalleles data/var_annot2seq.sft > tmp/prevar-sample3.tsv
	cg select -f {sequenced *} tmp/prevar-sample3.tsv tmp/var-sample3.tsv
	file copy data/sreg-annot1.sft tmp/sreg-sample1.tsv
	file copy data/sreg-annot2.sft tmp/sreg-sample2.tsv
	file copy data/sreg-annot2.sft tmp/sreg-sample3.tsv
	cg multicompar -split 1 tmp/temp.sft tmp/var-sample1.tsv tmp/var-sample2.tsv tmp/var-sample3.tsv
	cg multicompar_reannot -paged 2 -pagedstart 1 tmp/temp.sft
	exec diff tmp/temp.sft data/expected-multicompar-split-reannot-pagedstart1.sft
} {} 

test multicompar {basic, sequenced already present} {
	test_cleantmp
	cg multicompar -split 0 tmp/temp.sft data/var_annot.sft data/var_annot2seq.sft
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

#test multicompar {noalt} {
#	test_cleantmp
#	cg multicompar tmp/temp.sft data/var_annotnoalt.sft data/var_annot2noalt.sft
#	exec diff tmp/temp.sft data/expected-multicompar-var_annotvar_annot2noalt.sft
#} {} 

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
	cg multicompar -split 0 tmp/temp.sft tmp/var-varcaller1-mapper1-sample1.sft tmp/var-varcaller2-mapper2-sample1.sft tmp/var-varcaller1-mapper1-sample2.sft tmp/var-varcaller2-mapper2-sample2.sft
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

test multicompar {basic reannot varall} {
	test_cleantmp
	file copy data/var_annot.sft tmp/var-annot1.tsv
	file copy data/var_annot2.sft tmp/var-annot2.tsv
	file copy data/sreg-annot1.sft tmp/sreg-annot1.tsv
	file copy data/sreg-annot2.sft tmp/sreg-annot2.tsv
	catch {file delete tmp/temp.tsv}
	# make tmp/varall-annot1.tsv
	cg select -f {* sequenced="v"} data/var_annot.sft tmp/temp.tsv
	set f [open tmp/temp.tsv a]
	puts $f [join {chr2 4009 4010 snp C C C C teste 0.0 r} \t]
	puts $f [join {chr2 4010 4011 snp A G G C test7e 0.1 r} \t]
	close $f
	cg select -s - tmp/temp.tsv tmp/varall-annot1.tsv
	# make tmp/varall-annot2.tsv
	catch {file delete tmp/temp.tsv}
	cg select -f {* sequenced="v"} data/var_annot2.sft tmp/temp.tsv
	set f [open tmp/temp.tsv a]
	puts $f [join {chr1 4050 4060 snp G G G G test3e 0.3 r} \t]
	close $f
	cg select -s - tmp/temp.tsv tmp/varall-annot2.tsv
	catch {file delete tmp/temp.tsv}
	cg multicompar -split 0 tmp/temp.tsv tmp/var-annot1.tsv tmp/var-annot2.tsv
	cg multicompar_reannot tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-multicompar_reannot_varall-var_annotvar_annot2.sft
} {} 

test multicompar {basic reannot varall split} {
	test_cleantmp
	cg splitalleles data/var_annot.sft tmp/var-annot1.tsv
	cg splitalleles data/var_annot2.sft tmp/var-annot2.tsv
	file copy data/sreg-annot1.sft tmp/sreg-annot1.tsv
	file copy data/sreg-annot2.sft tmp/sreg-annot2.tsv
	catch {file delete tmp/temp.tsv}
	# make tmp/varall-annot1.tsv
	cg select -f {* sequenced="v"} data/var_annot.sft tmp/temp1.tsv
	cg splitalleles tmp/temp1.tsv tmp/temp.tsv
	set f [open tmp/temp.tsv a]
	puts $f [join {chr2 4009 4010 snp C C C C teste 0.0 r} \t]
	puts $f [join {chr2 4010 4011 snp A G G C test7e 0.1 r} \t]
	close $f
	cg select -s - tmp/temp.tsv tmp/varall-annot1.tsv
	# make tmp/varall-annot2.tsv
	catch {file delete tmp/temp1.tsv}
	cg select -f {* sequenced="v"} data/var_annot2.sft tmp/temp1.tsv
	cg splitalleles tmp/temp1.tsv tmp/temp.tsv
	set f [open tmp/temp.tsv a]
	puts $f [join {chr1 4050 4060 snp G G G G test3e 0.3 r} \t]
	close $f
	cg select -s - tmp/temp.tsv tmp/varall-annot2.tsv
	cg splitalleles data/expected-multicompar_reannot_varall-var_annotvar_annot2.sft tmp/expected.tsv
	file delete tmp/temp.tsv
	cg multicompar -split 1 tmp/temp.tsv tmp/var-annot1.tsv tmp/var-annot2.tsv
	cg multicompar_reannot tmp/temp.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
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
	cg multicompar -split 0 tmp/temp.tsv data/var-compartest1.tsv data/var-compartest2.tsv
	cg checksort tmp/temp.tsv
} {*error in "*var-compartest2.tsv": file uses split alleles ("*1 207806142 207806170 sub" occurs more than once and you are not running with the -split option)*} match error

test multicompar {error on badly sorted files} {
	test_cleantmp
	# this gave an incorrectly sorted file
	cg multicompar tmp/temp.tsv data/vars_sorterror1.sft data/vars3.sft
} {File (*vars_sorterror1.sft) is not correctly sorted (sort correctly using "cg select -s -")
chr10:43198434-43198435:snp:G came before chr3:52847042-52847060:del:*} match error

test multicompar {error on badly sorted files 2} {
	test_cleantmp
	# this gave an incorrectly sorted file
	cg multicompar tmp/temp.tsv data/vars_sorterror1.sft data/vars_sorterror2.sft
} {File (*vars_sorterror2.sft) is not correctly sorted (sort correctly using "cg select -s -")
chr3:52847303-52847304:snp:G came before chr3:52847042-52847060:del:*} match error

test multicompar {split reannot split multiallelic varall,no sreg,zyg} {
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
	cg multicompar -split 1 tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	cg multicompar_reannot tmp/temp.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test multicompar {split reannot split multiallelic varall,sreg,zyg} {
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
	cg multicompar -split 1 tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	cg multicompar_reannot tmp/temp.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test multicompar {split reannot split multiallelic only sreg,zyg} {
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
	cg multicompar -split 1 tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	cg multicompar_reannot tmp/temp.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test multicompar {split reannot split multiallelic, only sreg} {
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
	cg multicompar -split 1 tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	cg multicompar_reannot tmp/temp.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test multicompar {split reannot split multiallelic varall,sreg,zyg, check ref indels} {
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
	cg multicompar -split 1 tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	cg multicompar_reannot tmp/temp.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test multicompar {split reannot split multiallelic ins,sreg,zyg} {
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
	cg multicompar -split 1 tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	cg multicompar_reannot tmp/temp.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test multicompar {split reannot split multiallelic varall,sreg,zyg, check ref indels, overlapping snp} {
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
	cg multicompar -split 1 tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	cg multicompar_reannot tmp/temp.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test multicompar {targets} {
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
	cg multicompar -split 1 -targetvarsfile tmp/targets-testsnps.tsv tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	cg multicompar_reannot tmp/temp.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test multicompar {targets, no name} {
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
	cg multicompar -split 1 -targetvarsfile tmp/targets-testsnps.tsv tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	cg multicompar_reannot tmp/temp.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test multicompar {targets only} {
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
	cg multicompar -split 1 -targetvarsfile tmp/targets-testsnps.tsv tmp/temp.tsv
	cg multicompar_reannot tmp/temp.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test multicompar {different analyses} {
	test_cleantmp
	write_tab tmp/var1.tsv {
		chromosome	begin	end	type	ref	alt	zyg-gatk-sample1	alleleSeq1-gatk-sample1	alleleSeq2-gatk-sample1	score-gatk-sample1
		chr1	100	100	ins	{}	C	t	{}	C	30
	}
	write_tab tmp/var2.tsv {
		chromosome	begin	end	type	ref	alt	zyg-sam-sample1	alleleSeq1-sam-sample1	alleleSeq2-sam-sample1	score-sam-sample1
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	zyg-gatk-sample1	alleleSeq1-gatk-sample1	alleleSeq2-gatk-sample1	score-gatk-sample1	zyg-sam-sample1	alleleSeq1-sam-sample1	alleleSeq2-sam-sample1	score-sam-sample1
		1	100	100	ins	{}	C	t	{}	C	30	?	?	?	?
	}
	catch {file delete tmp/temp.tsv}
	cg multicompar -split 1 tmp/temp.tsv tmp/var1.tsv tmp/var2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test_cleantmp

set ::env(PATH) $keeppath

testsummarize
