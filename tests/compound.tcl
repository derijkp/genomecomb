#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test compound {basic} {
	write_tab tmp/vars.tsv {
		chromosome	begin	end	type	ref	alt	var	sequenced-a1-sample1 quality-a1-sample1 coverage-a1-sample1 sequenced-a1-sample2 quality-a1-sample2 coverage-a1-sample2 sequenced-a2-sample2 quality-a2-sample2 coverage-a2-sample2
		chr1	11600	11601	snp	C	A	v1	v	40	20	v	15	20	v	15	20
		chr1	12000	12001	snp	T	C	v2	v	50	25	v	20	30	v	20	30
		chr1	21600	21601	snp	A	C	v3	v	60	30	r	25	40	v	25	40
		chr1	22000	22001	del	G	{}	v4	v	70	50	v	30	50	v	30	50
	}
	write_tab tmp/gene_test.tsv {
		chrom	start	end	name	strand	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	name2
		chr1	11500	12100	cds1	+	11500	12050	1	11500,11900	11800,12100	cdsgene1
		chr1	21500	21800	rna	+	{}	{}	1	21500,	21800,	rnagene
		chr1	21500	21800	cds2	+	21500	21700	1	21500	21800	cdsgene2
		chr1	21500	22100	cds3	+	21500	22050	1	21500,21900	21800,22100	cdsgene3
	}
	exec cg annotate -dbdir $::refseqdir/hg19 tmp/vars.tsv tmp/atest.tsv tmp/gene_test.tsv
	foreach {test per samples criteria impact expected} {
		test1 analysis {} {$sequenced == "v" and $quality >= 15 and $coverage >= 20} { >=RNA} {
			var	transcripts	compound-a1-sample1	compound-a1-sample2	compound-a2-sample2
			v1	cdsgene1:cds1	2	2	2
			v2	cdsgene1:cds1	2	2	2
			v3	cdsgene3:cds3	2		2
			v4	cdsgene3:cds3	2		2
		}
		test2 analysis {} {$sequenced == "v" and $quality >= 15 and $coverage >= 20} { >=CDSMIS} {
			var	transcripts	compound-a1-sample1	compound-a1-sample2	compound-a2-sample2
			v3	cdsgene3:cds3	2		2
			v4	cdsgene3:cds3	2		2
		}
		test3 analysis {} {$sequenced == "v" and $quality >= 20 and $coverage >= 20} { >=RNA} {
			var	transcripts	compound-a1-sample1	compound-a1-sample2	compound-a2-sample2
			v1	cdsgene1:cds1	2		
			v2	cdsgene1:cds1	2		
			v3	cdsgene3:cds3	2		2
			v4	cdsgene3:cds3	2		2
		}
		test3samples analysis {a1-sample2 a2-sample2} {$sequenced == "v" and $quality >= 20 and $coverage >= 20} { >=RNA} {
			var	transcripts	compound-a1-sample2	compound-a2-sample2
			v3	cdsgene3:cds3		2
			v4	cdsgene3:cds3		2
		}
		test3samplesdup analysis {a2-sample2 a2-sample2} {$sequenced == "v" and $quality >= 20 and $coverage >= 20} { >=RNA} {
			var	transcripts	compound-a2-sample2
			v3	cdsgene3:cds3	2
			v4	cdsgene3:cds3	2
		}
		test4 analysis {} {$sequenced == "v" and $quality >= 100 and $coverage >= 20} { >=CDSMIS} {
			var	transcripts	compound-a1-sample1	compound-a1-sample2	compound-a2-sample2
		}
		test1s sample {} {$sequenced-a1 == "v" and $quality-a1 >= 15 and $coverage-a1 >= 20} { >=RNA} {
			var	transcripts	compound-sample1	compound-sample2
			v1	cdsgene1:cds1	2	2
			v2	cdsgene1:cds1	2	2
			v3	cdsgene3:cds3	2	
			v4	cdsgene3:cds3	2	
		}
		test1s sample sample2 {$sequenced-a1 == "v" and $quality-a1 >= 15 and $coverage-a1 >= 20} { >=RNA} {
			var	transcripts	compound-sample2
			v1	cdsgene1:cds1	2
			v2	cdsgene1:cds1	2
		}
	} {
		puts -nonewline "$test: "
		set expected [deindent $expected]
		# putsvars expected criteria impact
		cg compound -stack 1 -v 2 -per $per -samples $samples -geneset test -criteria $criteria -impact $impact tmp/atest.tsv tmp/result.tsv >@ stdout 2>@ stderr
		set result [exec cg select -f {var transcripts compound-*} tmp/result.tsv]
		if {$result ne $expected} {
			error "different result for $test:\n$result\nshould be:\n$expected"
		} else {
			puts "ok"
		}
	}
} {}

test compound {error on wrong -samples} {
	write_tab tmp/atest.tsv {
		chromosome	begin	end	type	ref	alt	var	sequenced-a1-sample1	quality-a1-sample1	coverage-a1-sample1	sequenced-a1-sample2	quality-a1-sample2	coverage-a1-sample2	sequenced-a2-sample2	quality-a2-sample2	coverage-a2-sample2	test_impact	test_gene	test_descr
		chr1	11600	11601	snp	C	A	v1	v	40	20	v	15	20	v	15	20	CDSMIS	cdsgene1	+cds1:exon1+101:c.101C>A:p.P34H
	}
	exec cg compound -stack 1 -v 2 -samples {sample1 wrongsample} -criteria {$sequenced-a1 eq "v"} -geneset test tmp/atest.tsv tmp/result.tsv
} {some samples given in -samples are not in the file */atest.tsv: wrongsample*} error match

test compound {error on wrong -analyses} {
	write_tab tmp/atest.tsv {
		chromosome	begin	end	type	ref	alt	var	sequenced-a1-sample1	quality-a1-sample1	coverage-a1-sample1	sequenced-a1-sample2	quality-a1-sample2	coverage-a1-sample2	sequenced-a2-sample2	quality-a2-sample2	coverage-a2-sample2	test_impact	test_gene	test_descr
		chr1	11600	11601	snp	C	A	v1	v	40	20	v	15	20	v	15	20	CDSMIS	cdsgene1	+cds1:exon1+101:c.101C>A:p.P34H
	}
	exec cg compound -stack 1 -v 2 -analyses {a1-wrongsample a1-sample1} -criteria {$sequenced-a1 eq "v"} -geneset test tmp/atest.tsv tmp/result.tsv
} {some samples given in -analyses are not in the file */atest.tsv: a1-wrongsample*} error match

testsummarize

