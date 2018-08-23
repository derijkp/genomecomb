#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl
package require genomecomb

if {![info exists argv]} {set argv {}}
set pos [lsearch $argv -d]
if {$pos != -1} {
	set distribute [lindex $argv [incr pos]]
	if {$distribute eq "direct"} {
		set tests {
			_direct {}
		}
	} else {
		set tests [subst {
			"_d$distribute" "-d $distribute"
		}]
	}
} else {
	set tests {
		_direct {}
		_d4 {-d 4}
	}
}

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

proc multicompartest {num} {
	set start 10
	set files {}
	set resultheader {chromosome begin end type ref alt}
	set resultline {1 100 101 snp A C}
	for {set i 1} {$i <= $num} {incr i} {
		lappend files tmp/var-sample$i.tsv
		set f [open tmp/var-sample$i.tsv w]
		set header {chromosome begin end type ref alt alleleSeq1 alleleSeq2}
		puts $f [join $header \t]
		puts $f [join {1 100 101 snp A C A C} \t]
		close $f
		lappend resultheader sequenced-sample$i alleleSeq1-sample$i alleleSeq2-sample$i
		lappend resultline v A C
	}
	set f [open tmp/expected.tsv w]
	puts $f [join $resultheader \t]
	puts $f [join $resultline \t]
	close $f
	return $files
}

set testname _direct
set jobopts {}

foreach {testname jobopts} $tests {

test pmulticompar$testname {basic} {
	test_cleantmp
	cg pmulticompar -split 0 {*}$::jobopts tmp/temp.sft data/var_annot.sft data/var_annot2.sft
	reorder data/expected-multicompar-var_annotvar_annot2.sft tmp/expected.tsv
	exec diff tmp/temp.sft tmp/expected.tsv
} {} 

test pmulticompar$testname {no sreg error} {
	test_cleantmp
	cg pmulticompar {*}$::jobopts -split 0 -i 0 tmp/temp.sft data/var_annot.sft data/var_annot2.sft
	reorder data/expected-multicompar-var_annotvar_annot2.sft tmp/expected.tsv
	exec diff tmp/temp.sft tmp/expected.tsv
} {*no sorted region file (*/sreg-var_annot.tsv) or varallfile (*/varall-var_annot.tsv) found: not properly processed sample*} error match

test pmulticompar$testname {basic with 3} {
	test_cleantmp
	foreach sample {annot annot2 annot3} {
		file copy data/var_$sample.sft tmp/var-$sample.tsv
	}
	cg pmulticompar -split 0 {*}$::jobopts tmp/temp.tsv tmp/var-annot.tsv tmp/var-annot2.tsv tmp/var-annot3.tsv 2> tmp/warnings.log
	reorder data/expected-multicompar-var_annotvar_annot3.sft tmp/expected.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {add to existing} {
	test_cleantmp
	foreach sample {annot annot2 annot3} {
		file copy data/var_$sample.sft tmp/var-$sample.tsv
	}
	exec cg pmulticompar -split 0 {*}$::jobopts tmp/temp.tsv tmp/var-annot.tsv tmp/var-annot2.tsv
	exec cg pmulticompar -split 0 {*}$::jobopts tmp/temp.tsv tmp/var-annot3.tsv
	reorder data/expected-multicompar-var_annotvar_annot3.sft tmp/expected.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {basic split} {
	test_cleantmp
	cg splitalleles data/var_annot.sft > tmp/var-sample1.tsv
	cg splitalleles data/var_annot2.sft > tmp/var-sample2.tsv
	cg splitalleles data/var_annot3.sft > tmp/var-sample3.tsv
	cg pmulticompar {*}$::jobopts -split 1 tmp/temp.sft tmp/var-sample1.tsv tmp/var-sample2.tsv tmp/var-sample3.tsv 2> tmp/warnings.log
	exec diff tmp/temp.sft data/expected-multicompar-reannot-split.tsv
} {} 

test pmulticompar$testname {add to existing split} {
	test_cleantmp
	cg splitalleles data/var_annot.sft > tmp/var-sample1.tsv
	cg splitalleles data/var_annot2.sft > tmp/var-sample2.tsv
	cg splitalleles data/var_annot3.sft > tmp/var-sample3.tsv
	cg pmulticompar {*}$::jobopts -split 1 tmp/temp.sft tmp/var-sample1.tsv tmp/var-sample2.tsv 2> tmp/warnings.log
	cg pmulticompar {*}$::jobopts -split 1 tmp/temp.sft tmp/var-sample3.tsv 2> tmp/warnings.log
	exec diff tmp/temp.sft data/expected-multicompar-reannot-split.tsv
} {} 

test pmulticompar$testname {basic reannot} {
	test_cleantmp
	cg select -f {* zyg=zyg("")} data/var_annot.sft tmp/var-annot1.tsv
	cg select -f {* zyg=zyg("")} data/var_annot2.sft tmp/var-annot2.tsv
	file copy data/sreg-annot1.sft tmp/sreg-annot1.tsv
	file copy data/sreg-annot2.sft tmp/sreg-annot2.tsv
	cg pmulticompar -split 0 {*}$::jobopts tmp/temp.tsv tmp/var-annot1.tsv tmp/var-annot2.tsv
	reorder data/expected-multicompar_reannot-var_annotvar_annot2.sft tmp/expected.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {basic split reannot} {
	test_cleantmp
	cg splitalleles data/var_annot.sft > tmp/var-sample1.tsv
	cg splitalleles data/var_annot2.sft > tmp/var-sample2.tsv
	cg splitalleles data/var_annot2seq.sft > tmp/prevar-sample3.tsv
	cg select -f {sequenced *} tmp/prevar-sample3.tsv tmp/var-sample3.tsv
	file copy data/sreg-annot1.sft tmp/sreg-sample1.tsv
	file copy data/sreg-annot2.sft tmp/sreg-sample2.tsv
	file copy data/sreg-annot2.sft tmp/sreg-sample3.tsv
	cg pmulticompar {*}$::jobopts -split 1 tmp/temp.sft tmp/var-sample1.tsv tmp/var-sample2.tsv tmp/var-sample3.tsv
	exec diff tmp/temp.sft data/expected-multicompar-split-reannot.sft
} {} 

test pmulticompar$testname {basic split reannot lz4} {
	test_cleantmp
	cg splitalleles data/var_annot.sft > tmp/var-sample1.tsv
	cg splitalleles data/var_annot2.sft > tmp/var-sample2.tsv
	cg splitalleles data/var_annot2seq.sft > tmp/prevar-sample3.tsv
	cg select -f {sequenced *} tmp/prevar-sample3.tsv tmp/var-sample3.tsv
	file copy data/sreg-annot1.sft tmp/sreg-sample1.tsv
	file copy data/sreg-annot2.sft tmp/sreg-sample2.tsv
	file copy data/sreg-annot2.sft tmp/sreg-sample3.tsv
	cg lz4 {*}[glob tmp/*]
	cg pmulticompar {*}$::jobopts -split 1 tmp/temp.tsv tmp/var-sample1.tsv.lz4 tmp/var-sample2.tsv.lz4 tmp/var-sample3.tsv.lz4
	exec diff tmp/temp.tsv data/expected-multicompar-split-reannot.sft
} {} 

test pmulticompar$testname {basic split reannot rz} {
	test_cleantmp
	cg splitalleles data/var_annot.sft > tmp/var-sample1.tsv
	cg splitalleles data/var_annot2.sft > tmp/var-sample2.tsv
	cg splitalleles data/var_annot2seq.sft > tmp/prevar-sample3.tsv
	cg select -f {sequenced *} tmp/prevar-sample3.tsv tmp/var-sample3.tsv
	file copy data/sreg-annot1.sft tmp/sreg-sample1.tsv
	file copy data/sreg-annot2.sft tmp/sreg-sample2.tsv
	file copy data/sreg-annot2.sft tmp/sreg-sample3.tsv
	cg razip {*}[glob tmp/*]
	cg pmulticompar {*}$::jobopts -split 1 tmp/temp.tsv tmp/var-sample1.tsv.rz tmp/var-sample2.tsv.rz tmp/var-sample3.tsv.rz
	exec diff tmp/temp.tsv data/expected-multicompar-split-reannot.sft
} {} 

test pmulticompar$testname {split reannot test diff alleles} {
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
	cg pmulticompar {*}$::jobopts -split 1 tmp/temp.sft tmp/var-sample1.tsv tmp/var-sample2.tsv
	exec diff tmp/temp.sft tmp/expected.sft
} {} 

test pmulticompar$testname {basic, sequenced already present} {
	test_cleantmp
	cg pmulticompar {*}$::jobopts -split 0 tmp/temp.sft data/var_annot.sft data/var_annot2seq.sft
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

test pmulticompar$testname {var and mapper naming convention} {
	test_cleantmp
	catch {file link tmp/var-varcaller1-mapper1-sample1.sft ../data/var_annot.sft}
	catch {file link tmp/var-varcaller2-mapper2-sample1.sft ../data/var_annot.sft}
	catch {file link tmp/var-varcaller1-mapper1-sample2.sft ../data/var_annot2.sft}
	catch {file link tmp/var-varcaller2-mapper2-sample2.sft ../data/var_annot2.sft}
	cg pmulticompar {*}$::jobopts -split 0 tmp/temp.sft tmp/var-varcaller1-mapper1-sample1.sft tmp/var-varcaller2-mapper2-sample1.sft tmp/var-varcaller1-mapper1-sample2.sft tmp/var-varcaller2-mapper2-sample2.sft
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

test pmulticompar$testname {basic reannot varall} {
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
	puts $f [join {chr2 4010 4011 snp A G G C test7e 0.1 u} \t]
	close $f
	cg select -s - tmp/temp.tsv tmp/varall-annot1.tsv
	mklink data/expected-multicompar_reannot_varall-var_annotvar_annot2.sft tmp/expected.tsv
	# make tmp/varall-annot2.tsv
	cg select -overwrite 1 -f {* sequenced="v"} data/var_annot2.sft tmp/temp.tsv
	set f [open tmp/temp.tsv a]
	puts $f [join {chr1 4050 4060 snp G G G G test3e 0.3 r} \t]
	close $f
	cg select -s - tmp/temp.tsv tmp/varall-annot2.tsv
	mklink data/expected-pmulticompar_reannot_varall-var_annotvar_annot2.tsv tmp/expected.tsv
	catch {file delete tmp/temp.tsv}
	cg pmulticompar {*}$::jobopts -split 0 tmp/temp.tsv tmp/var-annot1.tsv tmp/var-annot2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {basic reannot varall lz4 compressed} {
	test_cleantmp
	file copy data/var_annot.sft tmp/var-annot1.tsv
	file copy data/var_annot2.sft tmp/var-annot2.tsv
	file copy data/sreg-annot1.sft tmp/sreg-annot1.tsv
	file copy data/sreg-annot2.sft tmp/sreg-annot2.tsv
	catch {file delete tmp/temp.tsv}
	# make tmp/varall-annot1.tsv
	cg select -overwrite 1 -f {* sequenced="v"} data/var_annot.sft tmp/temp.tsv
	set f [open tmp/temp.tsv a]
	puts $f [join {chr2 4009 4010 snp C C C C teste 0.0 r} \t]
	puts $f [join {chr2 4010 4011 snp A G G C test7e 0.1 u} \t]
	close $f
	cg select -s - tmp/temp.tsv tmp/varall-annot1.tsv
	# make tmp/varall-annot2.tsv
	cg select -overwrite 1 -f {* sequenced="v"} data/var_annot2.sft tmp/temp.tsv
	set f [open tmp/temp.tsv a]
	puts $f [join {chr1 4050 4060 snp G G G G test3e 0.3 r} \t]
	close $f
	cg select -s - tmp/temp.tsv tmp/varall-annot2.tsv
	catch {file delete tmp/temp.tsv}
	cg lz4 {*}[glob tmp/*]
	mklink data/expected-pmulticompar_reannot_varall-var_annotvar_annot2.tsv tmp/expected.tsv
	cg pmulticompar {*}$::jobopts -split 0 tmp/temp.tsv tmp/var-annot1.tsv.lz4 tmp/var-annot2.tsv.lz4
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {basic reannot varall rz compressed} {
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
	puts $f [join {chr2 4010 4011 snp A G G C test7e 0.1 u} \t]
	close $f
	cg select -s - tmp/temp.tsv tmp/varall-annot1.tsv
	# make tmp/varall-annot2.tsv
	cg select -overwrite 1 -f {* sequenced="v"} data/var_annot2.sft tmp/temp.tsv
	set f [open tmp/temp.tsv a]
	puts $f [join {chr1 4050 4060 snp G G G G test3e 0.3 r} \t]
	close $f
	cg select -s - tmp/temp.tsv tmp/varall-annot2.tsv
	catch {file delete tmp/temp.tsv}
	cg razip {*}[glob tmp/*]
	mklink data/expected-pmulticompar_reannot_varall-var_annotvar_annot2.tsv tmp/expected.tsv
	cg pmulticompar {*}$::jobopts -split 0 tmp/temp.tsv tmp/var-annot1.tsv.rz tmp/var-annot2.tsv.rz
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {basic reannot varall split} {
	test_cleantmp
	cg splitalleles data/var_annot.sft tmp/var-annot1.tsv
	cg splitalleles data/var_annot2.sft tmp/var-annot2.tsv
	file copy data/sreg-annot1.sft tmp/sreg-annot1.tsv
	file copy data/sreg-annot2.sft tmp/sreg-annot2.tsv
	# make tmp/varall-annot1.tsv
	cg select -overwrite 1 -f {* sequenced="v"} data/var_annot.sft tmp/temp1.tsv
	cg splitalleles tmp/temp1.tsv tmp/temp.tsv
	set f [open tmp/temp.tsv a]
	puts $f [join {chr2 4009 4010 snp C C C C teste 0.0 r} \t]
	puts $f [join {chr2 4010 4011 snp A G G C test7e 0.1 u} \t]
	close $f
	cg select -s - tmp/temp.tsv tmp/varall-annot1.tsv
	# make tmp/varall-annot2.tsv
	catch {file delete tmp/temp.tsv}
	cg select -overwrite 1 -f {* sequenced="v"} data/var_annot2.sft tmp/temp1.tsv
	cg splitalleles tmp/temp1.tsv tmp/temp.tsv
	set f [open tmp/temp.tsv a]
	puts $f [join {chr1 4050 4060 snp G T G G test3e 0.3 r} \t]
	# puts $f [join {chr2 4001 4002 snp A C A A test8e 0.2 r} \t]
	close $f
	cg select -s - tmp/temp.tsv tmp/varall-annot2.tsv
	cg splitalleles data/expected-pmulticompar_reannot_varall-var_annotvar_annot2.tsv tmp/expected.tsv
	file delete tmp/temp.tsv
	cg pmulticompar {*}$::jobopts -split 1 tmp/temp.tsv tmp/var-annot1.tsv tmp/var-annot2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {basic reannot varall split lz4} {
	test_cleantmp
	cg splitalleles data/var_annot.sft tmp/var-annot1.tsv
	cg splitalleles data/var_annot2.sft tmp/var-annot2.tsv
	file copy data/sreg-annot1.sft tmp/sreg-annot1.tsv
	file copy data/sreg-annot2.sft tmp/sreg-annot2.tsv
	# make tmp/varall-annot1.tsv
	cg select -overwrite 1 -f {* sequenced="v"} data/var_annot.sft tmp/temp1.tsv
	cg splitalleles tmp/temp1.tsv tmp/temp.tsv
	set f [open tmp/temp.tsv a]
	puts $f [join {chr2 4009 4010 snp C C C C teste 0.0 r} \t]
	puts $f [join {chr2 4010 4011 snp A G G C test7e 0.1 u} \t]
	close $f
	cg select -s - tmp/temp.tsv tmp/varall-annot1.tsv
	# make tmp/varall-annot2.tsv
	cg select -overwrite 1 -f {* sequenced="v"} data/var_annot2.sft tmp/temp1.tsv
	cg splitalleles tmp/temp1.tsv tmp/temp.tsv
	set f [open tmp/temp.tsv a]
	puts $f [join {chr1 4050 4060 snp G T G G test3e 0.3 r} \t]
	# puts $f [join {chr2 4001 4002 snp A C A A test8e 0.2 r} \t]
	close $f
	cg select -s - tmp/temp.tsv tmp/varall-annot2.tsv
	file delete tmp/temp.tsv
	cg lz4 {*}[glob tmp/*]
	cg splitalleles data/expected-pmulticompar_reannot_varall-var_annotvar_annot2.tsv tmp/expected.tsv
	cg pmulticompar {*}$::jobopts -split 1 tmp/temp.tsv tmp/var-annot1.tsv.lz4 tmp/var-annot2.tsv.lz4
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {basic reannot varall split rz} {
	test_cleantmp
	cg splitalleles data/var_annot.sft tmp/var-annot1.tsv
	cg splitalleles data/var_annot2.sft tmp/var-annot2.tsv
	file copy data/sreg-annot1.sft tmp/sreg-annot1.tsv
	file copy data/sreg-annot2.sft tmp/sreg-annot2.tsv
	# make tmp/varall-annot1.tsv
	cg select -overwrite 1 -f {* sequenced="v"} data/var_annot.sft tmp/temp1.tsv
	cg splitalleles tmp/temp1.tsv tmp/temp.tsv
	set f [open tmp/temp.tsv a]
	puts $f [join {chr2 4009 4010 snp C C C C teste 0.0 r} \t]
	puts $f [join {chr2 4010 4011 snp A G G C test7e 0.1 u} \t]
	close $f
	cg select -s - tmp/temp.tsv tmp/varall-annot1.tsv
	# make tmp/varall-annot2.tsv
	cg select -overwrite 1 -f {* sequenced="v"} data/var_annot2.sft tmp/temp1.tsv
	cg splitalleles tmp/temp1.tsv tmp/temp.tsv
	set f [open tmp/temp.tsv a]
	puts $f [join {chr1 4050 4060 snp G T G G test3e 0.3 r} \t]
	# puts $f [join {chr2 4001 4002 snp A C A A test8e 0.2 r} \t]
	close $f
	cg select -s - tmp/temp.tsv tmp/varall-annot2.tsv
	file delete tmp/temp.tsv
	cg razip {*}[glob tmp/*]
	cg splitalleles data/expected-pmulticompar_reannot_varall-var_annotvar_annot2.tsv tmp/expected.tsv
	cg pmulticompar {*}$::jobopts -split 1 tmp/temp.tsv tmp/var-annot1.tsv.rz tmp/var-annot2.tsv.rz
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {basic reannot varall split gz} {
	test_cleantmp
	cg splitalleles data/var_annot.sft tmp/var-annot1.tsv
	cg splitalleles data/var_annot2.sft tmp/var-annot2.tsv
	file copy data/sreg-annot1.sft tmp/sreg-annot1.tsv
	file copy data/sreg-annot2.sft tmp/sreg-annot2.tsv
	# make tmp/varall-annot1.tsv
	cg select -overwrite 1 -f {* sequenced="v"} data/var_annot.sft tmp/temp1.tsv
	cg splitalleles tmp/temp1.tsv tmp/temp.tsv
	set f [open tmp/temp.tsv a]
	puts $f [join {chr2 4009 4010 snp C C C C teste 0.0 r} \t]
	puts $f [join {chr2 4010 4011 snp A G G C test7e 0.1 u} \t]
	close $f
	cg select -s - tmp/temp.tsv tmp/varall-annot1.tsv
	# make tmp/varall-annot2.tsv
	cg select -overwrite 1 -f {* sequenced="v"} data/var_annot2.sft tmp/temp1.tsv
	cg splitalleles tmp/temp1.tsv tmp/temp.tsv
	set f [open tmp/temp.tsv a]
	puts $f [join {chr1 4050 4060 snp G T G G test3e 0.3 r} \t]
	# puts $f [join {chr2 4001 4002 snp A C A A test8e 0.2 r} \t]
	close $f
	cg select -s - tmp/temp.tsv tmp/varall-annot2.tsv
	file delete tmp/temp.tsv
	exec gzip {*}[glob tmp/*]
	cg splitalleles data/expected-pmulticompar_reannot_varall-var_annotvar_annot2.tsv tmp/expected.tsv
	cg pmulticompar {*}$::jobopts -split 1 tmp/temp.tsv tmp/var-annot1.tsv.gz tmp/var-annot2.tsv.gz
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {basic reannot varall split gvcf} {
	test_cleantmp
	file copy data/varall-gatkh-bwa-sample1.gvcf.gz tmp
	file copy data/varall-gatkh-bwa-sample2.gvcf.gz tmp
	exec cg vcf2tsv tmp/varall-gatkh-bwa-sample1.gvcf.gz | cg select -q {$genoqual >= 10} | cg regjoin > tmp/sreg-gatkh-bwa-sample1.tsv
	exec cg vcf2tsv tmp/varall-gatkh-bwa-sample2.gvcf.gz | cg select -q {$genoqual >= 10}  | cg regjoin > tmp/sreg-gatkh-bwa-sample2.tsv
	cg gatk_index tmp/varall-gatkh-bwa-sample1.gvcf.gz tmp/varall-gatkh-bwa-sample2.gvcf.gz
	cg gatk_gatk_genotypevcfs -dbdir $::refseqdir/hg19 tmp/varall-gatkh-bwa-sample1.gvcf.gz tmp/var-gatkh-bwa-sample1.vcf.gz
	cg gatk_gatk_genotypevcfs -dbdir $::refseqdir/hg19 tmp/varall-gatkh-bwa-sample2.gvcf.gz tmp/var-gatkh-bwa-sample2.vcf.gz
	cg vcf2tsv -split 1 tmp/var-gatkh-bwa-sample1.vcf.gz tmp/var-gatkh-bwa-sample1.tsv.lz4
	cg vcf2tsv -split 1 tmp/var-gatkh-bwa-sample2.vcf.gz tmp/var-gatkh-bwa-sample2.tsv.lz4
	file delete tmp/result.tsv
	cg pmulticompar {*}$::jobopts -split 1 tmp/result.tsv tmp/var-gatkh-bwa-sample1.tsv.lz4 tmp/var-gatkh-bwa-sample2.tsv.lz4
	exec diff tmp/result.tsv data/expected-gatkh-pmulticompar.tsv
} {} 

test pmulticompar$testname {sort empty bug split} {
	test_cleantmp
	# this gave an incorrectly sorted file
	cg pmulticompar {*}$::jobopts -split 1 tmp/temp.tsv data/var-compartest1.tsv data/var-compartest2.tsv
	cg checksort tmp/temp.tsv
} {}

test pmulticompar$testname {error on split files without split option} {
	test_cleantmp
	# this gave an incorrectly sorted file
	cg pmulticompar {*}$::jobopts -split 0 tmp/temp.tsv data/var-compartest1.tsv data/var-compartest2.tsv
	# no error on command for -d num, check error file
	set temp [file_read tmp/temp.tsv.index/multicompar/log_jobs/multi_merge-vars.tsv.err]
	error $temp
} {*error in "*var-compartest2.tsv": file uses split alleles ("*1 207806142 207806170 sub" occurs more than once and you are not running with the -split option)*} match error

test pmulticompar$testname {error on badly sorted files} {
	test_cleantmp
	# this gave an incorrectly sorted file
	cg pmulticompar {*}$::jobopts -split 0 tmp/temp.tsv data/vars_sorterror1.sft data/var_annot2.sft
	# no error on command for -d num, check error file
	set temp [file_read tmp/temp.tsv.index/multicompar/log_jobs/multi_merge-vars.tsv.err]
	error $temp
} {File (*vars_sorterror1.sft) is not correctly sorted (sort correctly using "cg select -s -")
chr10:43198434-43198435:snp:G came before chr3:52847042-52847060:del:*} match error

test pmulticompar$testname {error on badly sorted files 2} {
	test_cleantmp
	# this gave an incorrectly sorted file
	cg pmulticompar {*}$::jobopts -split 0 tmp/temp.tsv data/vars_sorterror1.sft data/vars_sorterror2.sft
	# no error on command for -d num, check error file
	set temp [file_read tmp/temp.tsv.index/multicompar/log_jobs/multi_merge-vars.tsv.err]
	error $temp
} {File (*vars_sorterror2.sft) is not correctly sorted (sort correctly using "cg select -s -")
chr3:52847303-52847304:snp:G came before chr3:52847042-52847060:del:*} match error

test pmulticompar$testname {error on file with missing fields} {
	test_cleantmp
	write_tab tmp/errorfile.tsv {
		name amplicon   begin   end     typr    chr     ref     alt
		{}	{}	{}	{}	{}	{}	{}	{}
		ATAD3B_1421916  1421915 1421916 snp     chr1    T       C
	}
	# exec multi_merge 1 tmp/errorfile.tsv
	cg pmulticompar -split 0 tmp/result.tsv tmp/errorfile.tsv
} {field type not found in *errorfile.tsv
child process exited abnormally} match error

test pmulticompar$testname {split reannot split multiallelic varall,no sreg,zyg} {
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
	cg pmulticompar {*}$::jobopts -split 1 tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {split reannot split multiallelic varall,sreg,zyg} {
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
	cg pmulticompar {*}$::jobopts -split 1 tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {split reannot split multiallelic only sreg,zyg} {
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
	cg pmulticompar {*}$::jobopts -split 1 tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test pmulticompar$testname {split reannot split multiallelic, only sreg} {
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
	cg pmulticompar {*}$::jobopts -split 1 tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {split reannot split multiallelic varall,sreg,zyg, check ref indels} {
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
	cg pmulticompar {*}$::jobopts -split 1 tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {split reannot split multiallelic varall ins,sreg,zyg} {
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
	cg pmulticompar {*}$::jobopts -split 1 tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {split reannot split multiallelic varall,sreg,zyg, check ref indels, overlapping snp} {
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
		chr1	150	151	snp	G	G	r	G	G	30
		chr1	150	152	del	GT	{}	t	GT	{}	30
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
	cg pmulticompar {*}$::jobopts -split 1 tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {targets} {
	test_cleantmp
	write_tab tmp/sreg-sample1.tsv {
		chromosome	begin	end
		chr1	50	165
	}
	write_tab tmp/sreg-sample2.tsv {
		chromosome	begin	end
		chr1	50	200
	}
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
		chr1	150	151	snp	G	G	r	G	G	30
		chr1	150	152	del	GT	{}	t	GT	{}	30
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
		1	165	166	snp	N	G	u	u	?	?	?	r	r	N	N	?	snp4
	}
	catch {file delete tmp/temp.tsv}
	cg pmulticompar {*}$::jobopts -split 1 -targetvarsfile tmp/targets-testsnps.tsv tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {targets, no name} {
	test_cleantmp
	write_tab tmp/sreg-sample1.tsv {
		chromosome	begin	end
		chr1	50	200
	}
	write_tab tmp/sreg-sample2.tsv {
		chromosome	begin	end
		chr1	50	165
	}
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
		chr1	150	151	snp	G	G	r	G	G	30
		chr1	150	152	del	GT	{}	t	GT	{}	30
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
		1	165	166	snp	N	G	r	r	N	N	?	u	u	?	?	?	1
	}
	catch {file delete tmp/temp.tsv}
	cg pmulticompar {*}$::jobopts -split 1 -targetvarsfile tmp/targets-testsnps.tsv tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {add targets only} {
	test_cleantmp
	write_tab tmp/sreg-sample1.tsv {
		chromosome	begin	end
		chr1	50	200
	}
	file copy tmp/sreg-sample1.tsv tmp/sreg-sample2.tsv
	write_tab tmp/var-sample1.tsv {
		chromosome	begin	end	type	ref	alt	sequenced	zyg	alleleSeq1	alleleSeq2	score
		chr1	100	100	ins	{}	C	v	t	{}	C	40
		chr1	151	152	snp	A	T	v	t	A	T	32
		chr1	160	162	del	NN	{}	v	m	{}	{}	33
	}
	write_tab tmp/var-sample2.tsv {
		chromosome	begin	end	type	ref	alt	sequenced	zyg	alleleSeq1	alleleSeq2	score
		chr1	150	152	del	GT	{}	v	t	GT	{}	30
		chr1	160	161	snp	A	C	v	t	A	C	33
		chr1	160	161	snp	A	T	v	t	A	T	33
	}
	write_tab tmp/varall-sample1.tsv {
		chromosome	begin	end	type	ref	alt	sequenced	zyg	alleleSeq1	alleleSeq2	score
		chr1	100	100	ins	{}	C	v	t	{}	C	40
		chr1	100	101	snp	G	G	r	r	G	G	41
		chr1	150	151	snp	G	G	r	r	G	G	30
		chr1	151	152	snp	A	T	v	t	A	T	32
		chr1	152	153	snp	T	T	r	r	T	T	31
		chr1	160	162	del	NN	{}	v	m	{}	{}	33
	}
	write_tab tmp/varall-sample2.tsv {
		chromosome	begin	end	type	ref	alt	sequenced	zyg	alleleSeq1	alleleSeq2	score
		chr1	100	101	snp	G	G	r	r	G	G	41
		chr1	150	151	snp	G	G	r	r	G	G	30
		chr1	150	152	del	GT	{}	v	t	GT	{}	30
		chr1	151	152	snp	A	A	r	r	A	A	31
		chr1	152	153	snp	A	A	r	r	A	A	32
		chr1	160	161	snp	A	C	v	t	A	C	33
		chr1	160	161	snp	A	T	v	t	A	T	33
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
		1	165	166	snp	N	G	r	r	N	N	?	r	r	N	N	?	1
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
	cg pmulticompar {*}$::jobopts -split 1 -targetvarsfile tmp/targets-testsnps.tsv tmp/temp.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {different header for varall} {
	test_cleantmp
	write_tab tmp/sreg-sample1.tsv {
		chromosome	begin	end
		chr1	50	165
	}
	write_tab tmp/sreg-sample2.tsv {
		chromosome	begin	end
		chr1	50	200
	}
	# added sequenced, cluster not in varall; different order
	write_tab tmp/var-sample1.tsv {
		chromosome	begin	end	type	ref	new	alt	sequenced	zyg	alleleSeq1	alleleSeq2	score	cluster
		chr1	100	100	ins	{}	0	C	v	t	{}	C	40	{}
		chr1	151	152	snp	A	0	T	v	t	A	T	32	{}
		chr1	160	162	del	NN	1	{}	v	m	{}	{}	33	1
	}
	write_tab tmp/varall-sample1.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	extra	score
		chr1	100	100	ins	{}	C	t	{}	C	1	40
		chr1	100	101	snp	G	G	r	G	G	2	41
		chr1	150	151	snp	G	G	r	G	G	3	30
		chr1	151	152	snp	A	T	t	A	T	4	32
		chr1	152	153	snp	T	T	r	T	T	5	31
		chr1	160	162	del	NN	{}	m	{}	{}	6	33
	}
	write_tab tmp/var-sample2.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	score
		chr1	150	152	del	GT	{}	t	GT	{}	30
		chr1	160	161	snp	A	C	t	A	C	33
		chr1	160	161	snp	A	T	t	A	T	33
	}
	write_tab tmp/varall-sample2.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	score
		chr1	100	101	snp	G	G	r	G	G	41
		chr1	150	151	snp	G	G	r	G	G	30
		chr1	150	152	del	GT	{}	t	GT	{}	30
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
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	alleleSeq1-sample1	alleleSeq2-sample1	new-sample1	score-sample1	cluster-sample1	sequenced-sample2	zyg-sample2	alleleSeq1-sample2	alleleSeq2-sample2	score-sample2	testsnps
		1	100	100	ins	{}	C	v	t	{}	C	0	40	{}	r	r	{}	{}	41	{}
		1	100	101	snp	G	C	r	r	G	G	?	41	?	r	r	G	G	41	snp1
		1	150	152	del	GT	{}	r	r	GT	GT	?	30	?	v	t	GT	{}	30	snp2
		1	151	152	snp	A	T	v	t	A	T	0	32	{}	r	o	A	@	31	{}
		1	160	161	snp	A	C	r	o	@	@	?	?	?	v	t	A	C	33	{}
		1	160	161	snp	A	T	r	o	@	@	?	?	?	v	t	A	T	33	snp3
		1	160	162	del	NN	{}	v	m	{}	{}	1	33	1	r	r	NN	NN	33	{}
		1	165	166	snp	N	G	u	u	?	?	?	?	?	r	r	N	N	?	snp4
	}
	catch {file delete tmp/temp.tsv}
	cg pmulticompar {*}$::jobopts -split 1 -targetvarsfile tmp/targets-testsnps.tsv tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test pmulticompar$testname {varall missing one snp} {
	test_cleantmp
	write_tab tmp/sreg-sample1.tsv {
		chromosome	begin	end
		chr1	50	200
	}
	file copy tmp/sreg-sample1.tsv tmp/sreg-sample2.tsv
	write_tab tmp/var-sample1.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	score
		chr1	151	152	snp	A	T	t	A	T	31
	}
	write_tab tmp/var-sample2.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	score
	}
	write_tab tmp/varall-sample1.tsv {
		chromosome	begin	end	type	ref	alt	sequenced	alleleSeq1	alleleSeq2	score
		chr1	150	151	snp	G	G	r	G	G	30
		chr1	151	152	snp	A	T	t	A	T	31
		chr1	152	153	snp	T	T	r	T	T	32
	}
	write_tab tmp/varall-sample2.tsv {
		chromosome	begin	end	type	ref	alt	sequenced	alleleSeq1	alleleSeq2	score
		chr1	150	151	snp	G	G	r	G	G	40
		chr1	152	153	snp	A	A	r	A	A	42
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	alleleSeq1-sample1	alleleSeq2-sample1	score-sample1	sequenced-sample2	zyg-sample2	alleleSeq1-sample2	alleleSeq2-sample2	score-sample2
		1	151	152	snp	A	T	v	t	A	T	31	r	r	A	A	?
	}
	catch {file delete tmp/temp.tsv}
	cg pmulticompar {*}$::jobopts -split 1 tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {varall missing snp between chromsomes} {
	test_cleantmp
	write_tab tmp/sreg-sample1.tsv {
		chromosome	begin	end
		chr1	50	200
	}
	write_tab tmp/sreg-sample2.tsv {
		chromosome	begin	end
		chr1	50	151
	}
	write_tab tmp/var-sample1.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	score
		chr1	151	152	snp	A	T	t	A	T	31
	}
	write_tab tmp/var-sample2.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	score
	}
	write_tab tmp/varall-sample1.tsv {
		chromosome	begin	end	type	ref	alt	sequenced	alleleSeq1	alleleSeq2	score
		chr1	150	151	snp	G	G	r	G	G	30
		chr1	151	152	snp	A	T	t	A	T	31
		chr1	152	153	snp	T	T	r	T	T	32
	}
	write_tab tmp/varall-sample2.tsv {
		chromosome	begin	end	type	ref	alt	sequenced	alleleSeq1	alleleSeq2	score
		chr1	150	151	snp	G	G	r	G	G	40
		chr2	10	20	snp	A	A	r	A	A	42
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	alleleSeq1-sample1	alleleSeq2-sample1	score-sample1	sequenced-sample2	zyg-sample2	alleleSeq1-sample2	alleleSeq2-sample2	score-sample2
		1	151	152	snp	A	T	v	t	A	T	31	u	u	?	?	?
	}
	catch {file delete tmp/temp.tsv}
	cg pmulticompar {*}$::jobopts -split 1 tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {100 samples} {
	test_cleantmp
	set files [multicompartest 100]
	cg pmulticompar {*}$::jobopts -split 1 tmp/temp.tsv {*}$files
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {20 samples, -m 5} {
	test_cleantmp
	set files [multicompartest 20]
	cg pmulticompar {*}$::jobopts -m 5 -split 1 tmp/temp.tsv {*}$files 2>@ stderr >@ stdout
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {split reannot, multiallelic, bcol, sreg,zyg, check ref indels, overlapping snp} {
	test_cleantmp
	write_tab tmp/sreg-sample1.tsv {
		chromosome	begin	end
		chr1	50	200
		chr2	50	200
	}
	file copy tmp/sreg-sample1.tsv tmp/sreg-sample2.tsv
	write_tab tmp/var-sample1.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	coverage
		chr1	100	100	ins	{}	C	t	{}	C	21
		chr1	151	152	snp	A	T	t	A	T	33
		chr1	160	162	del	NN	{}	m	{}	{}	35
		chr2	151	152	snp	G	C	m	C	C	43
	}
	write_tab tmp/coverage-sample1.tsv {
		chromosome	begin	end	coverage
		chr1	100	101	31
		chr1	150	151	32
		chr1	151	152	33
		chr1	152	153	34
		chr1	160	162	35
		chr2	100	101	41
		chr2	150	151	42
		chr2	151	152	43
	}
	cg bcol_make tmp/coverage-sample1.bcol coverage tmp/coverage-sample1.tsv
	file delete tmp/coverage-sample1.tsv
	write_tab tmp/var-sample2.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	coverage
		chr1	150	152	del	GT	{}	t	GT	{}	22
		chr1	160	161	snp	A	C	t	A	C	55
		chr1	160	161	snp	A	T	t	A	T	55
	}
	write_tab tmp/coverage-sample2.tsv {
		chromosome	begin	end	coverage
		chr1	100	101	51
		chr1	150	151	52
		chr1	151	152	53
		chr1	152	153	54
		chr1	160	161	55
		chr1	160	161	55
		chr2	100	101	61
		chr2	150	151	62
		chr2	151	152	63
	}
	cg bcol_make tmp/coverage-sample2.bcol coverage tmp/coverage-sample2.tsv
	file delete coverage tmp/coverage-sample2.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	alleleSeq1-sample1	alleleSeq2-sample1	coverage-sample1	sequenced-sample2	zyg-sample2	alleleSeq1-sample2	alleleSeq2-sample2	coverage-sample2
		1	100	100	ins	{}	C	v	t	{}	C	21	r	r	{}	{}	51
		1	150	152	del	GT	{}	r	r	GT	GT	32	v	t	GT	{}	22
		1	151	152	snp	A	T	v	t	A	T	33	r	o	A	@	53
		1	160	161	snp	A	C	r	o	@	@	35	v	t	A	C	55
		1	160	161	snp	A	T	r	o	@	@	35	v	t	A	T	55
		1	160	162	del	NN	{}	v	m	{}	{}	35	r	r	NN	NN	55
		2	151	152	snp	G	C	v	m	C	C	43	r	r	G	G	63
	}
	catch {file delete tmp/temp.tsv}
	cg pmulticompar {*}$::jobopts -split 1 tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {split, old bcol} {
	test_cleantmp
	# sample1
	file mkdir tmp/sample1
	write_tab tmp/sample1/sreg-sample1.tsv {
		chromosome	begin	end
		chr1	20	50
		chr2	20	50
	}
	write_tab tmp/sample1/var-sample1.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	coverage
		chr1	35	36	snp	A	T	m	T	T	15
		chr1	40	41	snp	A	C	m	C	C	14
		chr1	40	42	del	GT	{}	t	GT	{}	14
	}
	file mkdir tmp/sample1/coverage
	file copy data/old-chr1.bcol tmp/sample1/coverage/coverage-1-sample1.bcol
	file copy data/old-chr1.bcol.bin.rz tmp/sample1/coverage/coverage-1-sample1.bcol.bin.rz
	file copy data/old-chr2.bcol tmp/sample1/coverage/coverage-2-sample1.bcol
	file copy data/old-chr2.bcol.bin.rz tmp/sample1/coverage/coverage-2-sample1.bcol.bin.rz
	# sample2
	file mkdir tmp/sample2
	write_tab tmp/sample2/var-sample2.tsv {
		chromosome	begin	end	type	ref	alt	sequenced	zyg	alleleSeq1	alleleSeq2	coverage
		chr1	35	35	ins	{}	C	v	t	{}	C	120
		chr1	35	36	snp	A	T	v	t	A	T	120
		chr1	160	161	snp	A	T	u	t	A	T	55
		chr2	49	50	snp	A	T	v	t	A	T	10
	}
	file copy tmp/sample1/sreg-sample1.tsv tmp/sample2/sreg-sample2.tsv
	file mkdir tmp/sample2/coverage
	file copy data/old-chr2.bcol tmp/sample2/coverage/coverage-1-sample2.bcol
	file copy data/old-chr2.bcol.bin.rz tmp/sample2/coverage/coverage-1-sample2.bcol.bin.rz
	file copy data/old-chr1.bcol tmp/sample2/coverage/coverage-2-sample2.bcol
	file copy data/old-chr1.bcol.bin.rz tmp/sample2/coverage/coverage-2-sample2.bcol.bin.rz
	# old
	# cg multicompar tmp/temp.tsv tmp/sample1/var-sample1.tsv tmp/sample2/var-sample2.tsv
	# cg multicompar_reannot tmp/temp.tsv
	# test vs expected
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	alleleSeq1-sample1	alleleSeq2-sample1	coverage-sample1	sequenced-sample2	zyg-sample2	alleleSeq1-sample2	alleleSeq2-sample2	coverage-sample2
		1	35	35	ins	{}	C	r	r	{}	{}	15	v	t	{}	C	120
		1	35	36	snp	A	T	v	m	T	T	15	v	t	A	T	120
		1	40	41	snp	A	C	v	m	C	C	14	r	r	A	A	24
		1	40	42	del	GT	{}	v	t	GT	{}	14	r	r	GT	GT	24
		1	160	161	snp	A	T	u	u	?	?	0	u	t	A	T	55
		2	49	50	snp	A	T	r	r	A	A	25	v	t	A	T	10
	}
	catch {file delete tmp/temp.tsv}
	cg pmulticompar {*}$::jobopts -split 1 tmp/temp.tsv tmp/sample1/var-sample1.tsv tmp/sample2/var-sample2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {split reannot, multiallelic, regfile} {
	test_cleantmp
	write_tab tmp/sreg-sample1.tsv {
		chromosome	begin	end
		chr1	50	200
		chr2	50	200
	}
	file copy tmp/sreg-sample1.tsv tmp/sreg-sample2.tsv
	write_tab tmp/var-sample1.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	coverage
		chr1	100	100	ins	{}	C	t	{}	C	1
		chr1	151	152	snp	A	T	t	A	T	1
		chr1	160	162	del	NN	{}	m	{}	{}	{}
		chr2	150	151	snp	G	C	m	C	C	{}
		chr2	151	152	snp	G	C	m	C	C	{}
	}
	write_tab tmp/reg_coverage-sample1.tsv {
		chromosome	begin	end
		chr1	100	155
	}
	write_tab tmp/var-sample2.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	coverage
		chr1	150	152	del	GT	{}	t	GT	{}	2
		chr1	160	161	snp	A	C	t	A	C	{}
		chr1	160	161	snp	A	T	t	A	T	{}
	}
	write_tab tmp/reg_coverage-sample2.tsv {
		chromosome	begin	end	score
		chr1	101	162	2
		chr2	100	102	3
		chr2	120	151	4
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	alleleSeq1-sample1	alleleSeq2-sample1	coverage-sample1	sequenced-sample2	zyg-sample2	alleleSeq1-sample2	alleleSeq2-sample2	coverage-sample2
		1	100	100	ins	{}	C	v	t	{}	C	1	r	r	{}	{}	{}
		1	150	152	del	GT	{}	r	r	GT	GT	1	v	t	GT	{}	2
		1	151	152	snp	A	T	v	t	A	T	1	r	o	A	@	2
		1	160	161	snp	A	C	r	o	@	@	{}	v	t	A	C	{}
		1	160	161	snp	A	T	r	o	@	@	{}	v	t	A	T	{}
		1	160	162	del	NN	{}	v	m	{}	{}	{}	r	r	NN	NN	2
		2	150	151	snp	G	C	v	m	C	C	{}	r	r	G	G	4
		2	151	152	snp	G	C	v	m	C	C	{}	r	r	G	G	{}
	}
	catch {file delete tmp/temp.tsv}
	cg pmulticompar {*}$::jobopts -split 1 tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {split reannot, multiallelic, regfile} {
	test_cleantmp
	write_tab tmp/sreg-sample1.tsv {
		chromosome	begin	end
		chr1	50	200
		chr2	50	200
	}
	file copy tmp/sreg-sample1.tsv tmp/sreg-sample2.tsv
	write_tab tmp/var-sample1.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	coverage
		chr1	100	100	ins	{}	C	t	{}	C	1
		chr1	151	152	snp	A	T	t	A	T	1
		chr1	160	162	del	NN	{}	m	{}	{}	{}
		chr2	150	151	snp	G	C	m	C	C	{}
		chr2	151	152	snp	G	C	m	C	C	{}
	}
	write_tab tmp/reg_coverage-sample1.tsv {
		chromosome	begin	end
		chr1	100	155
	}
	write_tab tmp/var-sample2.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	coverage
		chr1	150	152	del	GT	{}	t	GT	{}	2
		chr1	160	161	snp	A	C	t	A	C	{}
		chr1	160	161	snp	A	T	t	A	T	{}
	}
	write_tab tmp/reg_coverage-sample2.tsv {
		chromosome	begin	end	score
		chr1	101	162	2
		chr2	100	102	3
		chr2	120	151	4
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	alleleSeq1-sample1	alleleSeq2-sample1	coverage-sample1	sequenced-sample2	zyg-sample2	alleleSeq1-sample2	alleleSeq2-sample2	coverage-sample2
		1	100	100	ins	{}	C	v	t	{}	C	1	r	r	{}	{}	{}
		1	150	152	del	GT	{}	r	r	GT	GT	1	v	t	GT	{}	2
		1	151	152	snp	A	T	v	t	A	T	1	r	o	A	@	2
		1	160	161	snp	A	C	r	o	@	@	{}	v	t	A	C	{}
		1	160	161	snp	A	T	r	o	@	@	{}	v	t	A	T	{}
		1	160	162	del	NN	{}	v	m	{}	{}	{}	r	r	NN	NN	2
		2	150	151	snp	G	C	v	m	C	C	{}	r	r	G	G	4
		2	151	152	snp	G	C	v	m	C	C	{}	r	r	G	G	{}
	}
	catch {file delete tmp/temp.tsv}
	cg pmulticompar {*}$::jobopts -split 1 tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {coverageRefScore tsv file} {
	test_cleantmp
	file mkdir tmp/sample1/coverage
	write_tab tmp/sample1/var-sample1.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	coverage	refscore
		chr21	16050007	16050008	snp	A	T	t	A	T	1	5
	}
	write_tab tmp/sample1/sreg-sample1.tsv {
		chromosome	begin	end
		chr21	16050000	16050100
		chr22	16050000	16050100
	}
	file copy data/cgtest.tsv tmp/sample1/coverage/coverageRefScore-21-sample1.tsv
	exec gzip tmp/sample1/coverage/coverageRefScore-21-sample1.tsv
	file copy data/cgtest.tsv tmp/sample1/coverage/coverageRefScore-22-sample1.tsv
	exec gzip tmp/sample1/coverage/coverageRefScore-22-sample1.tsv
	file mkdir tmp/sample2/coverage
	write_tab tmp/sample2/var-sample2.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	coverage
		chr22	16050008	16050009	snp	G	C	t	G	C	4
	}
	file copy tmp/sample1/sreg-sample1.tsv tmp/sample2/sreg-sample2.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	alleleSeq1-sample1	alleleSeq2-sample1	coverage-sample1	refscore-sample1	sequenced-sample2	zyg-sample2	alleleSeq1-sample2	alleleSeq2-sample2	coverage-sample2
		21	16050007	16050008	snp	A	T	v	t	A	T	1	5	r	r	A	A	3
		22	16050008	16050009	snp	G	C	r	r	G	G	2	10	v	t	G	C	4
	}
	file copy data/cgtest.tsv tmp/sample2/coverage/coverageRefScore-21-sample2.tsv
	exec gzip tmp/sample2/coverage/coverageRefScore-21-sample2.tsv
	file copy data/cgtest.tsv tmp/sample2/coverage/coverageRefScore-22-sample2.tsv
	exec gzip tmp/sample2/coverage/coverageRefScore-22-sample2.tsv
	catch {file delete tmp/temp.tsv}
	cg pmulticompar {*}$::jobopts -split 1 tmp/temp.tsv tmp/sample1 tmp/sample2
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {varall coverage with del} {
	test_cleantmp
	write_tab tmp/var-sample1.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	coverage	totalcoverage
		chr1	100	101	snp	T	C	t	T	C	40	41
	}
	write_tab tmp/varall-sample1.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	coverage	totalcoverage
		chr1	100	101	snp	T	C	t	T	C	40	41
	}
	write_tab tmp/var-sample2.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	coverage	totalcoverage
		chr1	100	101	del	T	{}	m	{}	{}	30	31
	}
	write_tab tmp/varall-sample2.tsv {
		chromosome	begin	end	type	ref	alt	zyg	alleleSeq1	alleleSeq2	coverage	totalcoverage
		chr1	100	101	del	T	{}	m	{}	{}	30	31
		chr1	100	101	snp	T	.	t	T	T	30	31
	}
	write_tab tmp/sreg-sample1.tsv {
		chromosome	begin	end
		chr1	50	200
	}
	file copy -force tmp/sreg-sample1.tsv tmp/sreg-sample2.tsv
	# deletions overlapping with a snp variant (sample2) are properly caught
	# and assigned zyg=o and alleles @ (meaning something different from reference)
	# for snps overlapping a deletion (sample1) this is not checked (would require lookahead)
	# so the call for the del in sample1 is r T T instead of the more correct o T C
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	alleleSeq1-sample1	alleleSeq2-sample1	coverage-sample1	totalcoverage-sample1	sequenced-sample2	zyg-sample2	alleleSeq1-sample2	alleleSeq2-sample2	coverage-sample2	totalcoverage-sample2
		1	100	101	del	T	{}	r	r	T	T	40	41	v	m	{}	{}	30	31
		1	100	101	snp	T	C	v	t	T	C	40	41	r	o	@	@	30	31
	}
	file delete tmp/temp.tsv
	cg pmulticompar {*}$::jobopts -split 1 tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test pmulticompar$testname {unsplit, alt .} {
	test_cleantmp
	write_tab tmp/var-sample1.tsv {
		chromosome	begin	end	type	ref	alt	sequenced	alleleSeq1	alleleSeq2
		chr1	100	101	snp	A	.	r	A	A
		chr1	101	102	snp	G	.	r	G	G
		chr1	104	105	snp	G	C	v	C	C
	}
	write_tab tmp/var-sample2.tsv {
		chromosome	begin	end	type	ref	alt	sequenced	alleleSeq1	alleleSeq2
		chr1	100	101	snp	A	C	v	A	C
		chr1	101	102	snp	G	.	r	G	G
	}
	write_tab tmp/sreg-sample1.tsv {
		chromosome	begin	end
		chr1	50	200
	}
	file copy tmp/sreg-sample1.tsv tmp/sreg-sample2.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	alleleSeq1-sample1	alleleSeq2-sample1	sequenced-sample2	alleleSeq1-sample2	alleleSeq2-sample2
		1	100	101	snp	A	C	r	A	A	v	A	C
		1	101	102	snp	G	.	r	G	G	r	G	G
		1	104	105	snp	G	C	v	C	C	r	G	G
	}
	catch {file delete tmp/temp.tsv}
	cg pmulticompar {*}$::jobopts -split 0 tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

test pmulticompar$testname {unsplit, target with alt .} {
	test_cleantmp
	write_tab tmp/targets.tsv {
		chromosome	begin	end	type	ref	alt
		chr1	100	101	snp	A	.
		chr1	101	102	snp	G	.
	}
	write_tab tmp/var-sample1.tsv {
		chromosome	begin	end	type	ref	alt	alleleSeq1	alleleSeq2
	}
	write_tab tmp/var-sample2.tsv {
		chromosome	begin	end	type	ref	alt	alleleSeq1	alleleSeq2
		chr1	100	101	snp	A	C	A	C
	}
	write_tab tmp/sreg-sample1.tsv {
		chromosome	begin	end
		chr1	50	200
	}
	file copy tmp/sreg-sample1.tsv tmp/sreg-sample2.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	alleleSeq1-sample1	alleleSeq2-sample1	sequenced-sample2	alleleSeq1-sample2	alleleSeq2-sample2	targets
		1	100	101	snp	A	C	r	A	A	v	A	C	{}
		1	101	102	snp	G	.	r	G	G	r	G	G	1
	}
	catch {file delete tmp/temp.tsv}
	cg pmulticompar {*}$::jobopts -split 0 -t tmp/targets.tsv tmp/temp.tsv tmp/var-sample1.tsv tmp/var-sample2.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {} 

}

test_cleantmp

set ::env(PATH) $keeppath

testsummarize
