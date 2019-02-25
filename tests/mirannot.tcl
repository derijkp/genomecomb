#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test mir_annot {basic mir annotation} {
	test_cleantmp
	file_write tmp/vars_mirna.tsv "#test\n#comments\n"
	exec cat data/vars_mirna.tsv >> tmp/vars_mirna.tsv
	exec cg annotate -dbdir $::refseqdir/hg19 tmp/vars_mirna.tsv tmp/annot_test.tsv data/mir_small.tsv
	exec diff tmp/annot_test.tsv data/expected-vars_mirna_annot.tsv
} {} 

test mir_annot {basic mir annotation without transcript col} {
	test_cleantmp
	cg select -rf transcript data/mir_small.tsv tmp/mir_small.tsv
	cg select -q {$ROW > 15} data/vars_mirna.tsv tmp/vars_mirna.tsv
	cg select -rc 1 -q {$ROW > 15} data/expected-vars_mirna_annot.tsv tmp/expected.tsv
	exec cg annotate tmp/vars_mirna.tsv tmp/annot_test.tsv tmp/mir_small.tsv
	exec diff tmp/annot_test.tsv tmp/expected.tsv
} {} 

test mir_annot {basic mir annotation same name transcipt} {
	test_cleantmp
	write_tab tmp/mir_small.tsv {
		chromosome	begin	end	strand	name	transcript	mature1start	mature1end	loopstart	loopend	mature2start	mature2end
		chr1	567704	567793	+	mir1+	mir1+	{}	{}	567749	567754	567761	567783
		chr1	567704	567793	+	mir1+	mir1+	567704	567749	567749	567754	567761	567783
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	name	small_impact	small_mir
		chr1	567603	567604	snp	G	C	pre	mir1+:upstream(a-101e)	mir1+
		chr1	567704	567705	snp	A	C	arm-5p-1	mir1+:arm5p(l-e45);mir1+:mature5p1_46(a+e1)	mir1+
		chr1	567760	567761	snp	T	C	arm-3p-end	mir1+:arm3p(m-1e)	mir1+
	}
	cg select -q {$name in {pre arm-5p-1 arm-3p-end}} data/vars_mirna.tsv tmp/vars_mirna.tsv
	exec cg annotate tmp/vars_mirna.tsv tmp/annot_test.tsv tmp/mir_small.tsv
	exec diff tmp/annot_test.tsv tmp/expected.tsv
} {}

test mir_annot {different start mixes upstream and hit} {
	test_cleantmp
	write_tab tmp/mir_small.tsv {
		chromosome	begin	end	strand	name	transcript	mature1start	mature1end	loopstart	loopend	mature2start	mature2end
		chr1	567404	567493	+	mir2+	mir2+	{}	{}	567449	567454	567461	567483
		chr1	567704	567793	+	mir1+	mir1+	{}	{}	567749	567754	567761	567783
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	name	small_impact	small_mir
		chr1	567603	567604	snp	G	C	pre	mir2+:downstream(a+111);mir1+:upstream(a-101e)	mir2+;mir1+
		chr1	567704	567705	snp	A	C	arm-5p-1	mir1+:arm5p(l-e45)	mir1+
		chr1	567760	567761	snp	T	C	arm-3p-end	mir1+:arm3p(m-1e)	mir1+
	}
	cg select -q {$name in {pre arm-5p-1 arm-3p-end}} data/vars_mirna.tsv tmp/vars_mirna.tsv
	exec cg annotate tmp/vars_mirna.tsv tmp/annot_test.tsv tmp/mir_small.tsv
	exec diff tmp/annot_test.tsv tmp/expected.tsv
} {}

test mir_annot {mir annotation with status} {
	test_cleantmp
	write_tab tmp/mir_small.tsv {
		chromosome	begin	end	strand	name	transcript	mature1start	mature1end	loopstart	loopend	mature2start	mature2end status
		chr1	567704	567793	+	mir1+	mir1+	{}	{}	567749	567754	567761	567783 v
		chr1	567704	567793	+	mir1b+	mir1b+	567704	567749	567749	567754	567761	567783 p
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	name	small_impact	small_mir small_status
		chr1	567603	567604	snp	G	C	pre	mir1+:upstream(a-101e);mir1b+:upstream(a-101e)	mir1+;mir1b+ v;p
		chr1	567704	567705	snp	A	C	arm-5p-1	mir1+:arm5p(l-e45);mir1b+:mature5p1_46(a+e1)	mir1+;mir1b+ v;p
		chr1	567760	567761	snp	T	C	arm-3p-end	mir1+:arm3p(m-1e);mir1b+:arm3p(m-1e)	mir1+;mir1b+ v;p
	}
	cg select -q {$name in {pre arm-5p-1 arm-3p-end}} data/vars_mirna.tsv tmp/vars_mirna.tsv
	exec cg annotate tmp/vars_mirna.tsv tmp/annot_test.tsv tmp/mir_small.tsv
	exec diff tmp/annot_test.tsv tmp/expected.tsv
} {}

test mir_annot {-upstream} {
	test_cleantmp
	file copy data/vars_mirna.tsv tmp/vars_mirna.tsv
	exec cg annotate --upstreamsize 100 -dbdir $::refseqdir/hg19 tmp/vars_mirna.tsv tmp/annot_test.tsv.temp data/mir_small.tsv
	cg select -rc 1 tmp/annot_test.tsv.temp tmp/annot_test.tsv
	cg select -rc 1 data/expected-vars_mirna_annot.tsv tmp/expected.tsv
	exec diff tmp/annot_test.tsv tmp/expected.tsv
} {2c2
< chr1	567603	567604	snp	G	C	pre		
---
> chr1	567603	567604	snp	G	C	pre	mir1+:upstream(a-101e);mir1-a:downstream(a+e101);mir1-b:downstream(a+e101)	mir1+;mir1-;mir1-
16c16
< chr1	567893	567894	snp	A	C	post		
---
> chr1	567893	567894	snp	A	C	post	mir1+:downstream(a+e101);mir1-a:upstream(a-101e);mir1-b:upstream(a-101e)	mir1+;mir1-;mir1-
*} error match

set ::env(PATH) $keeppath

testsummarize
