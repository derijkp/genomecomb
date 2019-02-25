#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test reg_annot {basic} {
	file copy data/vars1.sft tmp/vars1.sft
	exec cg annotate tmp/vars1.sft tmp/temp.sft data/reg_annot.sft
	exec cg select -rf {list} tmp/temp.sft tmp/temp2.sft
	exec diff tmp/temp2.sft data/expected-vars1-reg_annot.sft
} {} 

test reg_annot {error file does not exists} {
	file copy data/vars1.sft tmp/vars1.sft
	exec cg annotate tmp/vars1.sft tmp/temp.sft data/reg_doesnotexist.tsv
} {File data/reg_doesnotexist.tsv does not exist} error

test reg_annot {compressed} {
	file copy data/vars1.sft tmp/vars1.sft
	file copy -force data/reg_annot.sft data/reg_annot.sft.opt tmp
	exec cg razip tmp/reg_annot.sft
	exec cg annotate tmp/vars1.sft tmp/temp.sft tmp/reg_annot.sft.rz
	exec cg select -rf {list} tmp/temp.sft tmp/temp2.sft
	exec diff tmp/temp2.sft data/expected-vars1-reg_annot.sft
} {} 

test reg_annot {2 compressed} {
	file copy data/vars1.sft tmp/vars1.sft
	file copy -force tmp/vars1.sft data/reg_annot.sft data/reg_annot.sft.opt tmp
	exec cg razip tmp/reg_annot.sft tmp/vars1.sft
	exec cg annotate tmp/vars1.sft.rz tmp/temp.sft tmp/reg_annot.sft.rz
	exec cg select -rf {list} tmp/temp.sft tmp/temp2.sft
	exec diff tmp/temp2.sft data/expected-vars1-reg_annot.sft
} {}

test reg_annot {basic, multiple fields} {
	file copy data/vars1.sft tmp/vars1.sft
	cg select -f {chromosome begin end type ref alt} tmp/vars1.sft tmp/vars.sft
	file copy -force data/reg_annot.sft tmp/reg_annot.sft
	file_write tmp/reg_annot.sft.opt "fields\t{type begin end}\n"
	exec cg annotate tmp/vars.sft tmp/temp.sft tmp/reg_annot.sft
	exec diff tmp/temp.sft data/expected-vars1-reg_annot-multi.sft
} {} 

test reg_annot {near} {
	file copy data/vars1.sft tmp/vars1.sft
	exec cg annotate -near 1000 tmp/vars1.sft tmp/temp.sft data/reg_annot.sft
	exec cg select -rf {list} tmp/temp.sft tmp/temp2.sft
	exec diff tmp/temp2.sft data/expected_near-vars1-reg_annot.sft
} {} 

test reg_annot {near indels} {
	file copy data/vars1.sft tmp/vars1.sft
	exec cg select -q {$type == "del" || $type == "ins"} tmp/vars1.sft tmp/indels1.sft
	exec cg annotate -near 50 tmp/vars1.sft tmp/temp.sft tmp/indels1.sft
	exec cg select -rf {list} tmp/temp.sft tmp/temp2.sft
	exec diff tmp/temp2.sft data/expected-vars1-indels.sft
} {}

test reg_annot {sort error 1 in vars} {
	file copy data/vars_sorterror1.sft tmp/vars_sorterror1.sft
	exec cg annotate tmp/vars_sorterror1.sft tmp/temp.sft data/reg_annot.sft
} {*File (*) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test reg_annot {sort error 2 in vars} {
	file copy data/vars_sorterror2.sft tmp/vars_sorterror2.sft
	exec cg annotate tmp/vars_sorterror2.sft tmp/temp.sft data/reg_annot.sft
} {*File (*) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test reg_annot {sort error 3 in vars} {
	file copy data/vars_sorterror3.sft tmp/vars_sorterror3.sft
	exec cg annotate tmp/vars_sorterror3.sft tmp/temp.sft data/reg_annot.sft
} {*File (*) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test reg_annot {no sort error 4 in vars (not relevant for reg)} {
	file copy data/vars_sorterror4.sft tmp/vars_sorterror4.sft
	catch {exec cg annotate tmp/vars_sorterror4.sft tmp/temp.sft data/reg_annot.sft}
} 0

test reg_annot {no sort error 5 in vars (not relevant for reg)} {
	file copy data/vars_sorterror5.sft tmp/vars_sorterror5.sft
	catch {exec cg annotate tmp/vars_sorterror5.sft tmp/temp.sft data/reg_annot.sft}
} 0

test reg_annot {sort error 1 in database file} {
	file copy -force data/vars_sorterror1.sft tmp/reg_annot.sft
	file copy data/vars_annottest.sft tmp/vars_annottest.sft
	exec cg annotate tmp/vars_annottest.sft tmp/temp.sft tmp/reg_annot.sft
} {*File (database file) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test reg_annot {sort error 2 in database file} {
	file copy -force data/vars_sorterror2.sft tmp/reg_annot.sft
	file copy data/vars_annottest.sft tmp/vars_annottest.sft
	exec cg annotate tmp/vars_annottest.sft tmp/temp.sft tmp/reg_annot.sft
} {*File (database file) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test reg_annot {sort error 3 in database file} {
	file copy -force data/vars_sorterror3.sft tmp/reg_annot.sft
	file copy data/vars_annottest.sft tmp/vars_annottest.sft
	exec cg annotate tmp/vars_annottest.sft tmp/temp.sft tmp/reg_annot.sft
} {*File (database file) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test reg_annot {bug fix deal with duplicate field in opt} {
	file copy data/vars1.sft tmp/vars1.sft
	cg select -f {chromosome begin end type ref alt} tmp/vars1.sft tmp/vars.sft
	file copy -force data/reg_annot.sft tmp/reg_annot.sft
	file_write tmp/reg_annot.sft.opt "fields\t{type begin end begin}\n"
	exec cg annotate tmp/vars.sft tmp/temp.sft tmp/reg_annot.sft
	exec diff tmp/temp.sft data/expected-vars1-reg_annot-multi.sft
} {} 

test var_annot {basic} {
	file copy data/vars1.sft tmp/vars1.sft
	exec cg annotate tmp/vars1.sft tmp/temp.sft data/var_annot.sft
	exec cg select -rf {list} tmp/temp.sft tmp/temp2.sft
	exec diff tmp/temp2.sft data/expected-vars1-var_annot.sft
} {} 

test var_annot {basic splitvar} {
	file copy data/vars1.sft tmp/vars1.sft
	exec cg splitalleles tmp/vars1.sft > tmp/vars1.tsv
	exec cg annotate tmp/vars1.tsv tmp/temp.tsv data/var_annot.sft
	exec cg select -rf {list} tmp/temp.tsv tmp/temp2.tsv
	exec cg splitalleles data/expected-vars1-var_annot.sft > tmp/expected.tsv
	exec diff tmp/temp2.tsv tmp/expected.tsv
} {} 

test var_annot {basic multi} {
	file mkdir tmp
	cg select -f {chromosome begin end type ref alt} data/vars1.sft tmp/vars.sft
	file copy -force data/var_annot.sft tmp/var_annot.sft
	file_write tmp/var_annot.sft.opt "fields\t{name freq alt}\n"
	exec cg annotate tmp/vars.sft tmp/temp.sft tmp/var_annot.sft
	exec diff tmp/temp.sft data/expected-vars1-var_annot-multi.sft
} {} 

test var_annot {lz4, opt, links} {
	test_cleantmp
	file mkdir tmp
	cg select -f {chromosome begin end type ref alt} data/vars1.sft tmp/vars.sft
	exec lz4c -q -c data/var_annot.sft > tmp/var_annot.sft.lz4
	file_write tmp/var_annot.sft.opt "fields\t{name freq alt}\n"
	cd tmp
	mklink var_annot.sft.lz4 var_annot.tsv.lz4
	mklink var_annot.sft.opt var_annot.tsv.opt
	cd ..
	exec cg annotate tmp/vars.sft tmp/temp.sft tmp/var_annot.tsv.lz4
	exec diff tmp/temp.sft data/expected-vars1-var_annot-multi.sft
} {} 

test var_annot {different types on same pos} {
	file copy data/vars2.tsv tmp/vars2.tsv
	exec cg annotate tmp/vars2.tsv tmp/temp.tsv data/var_annot3.tsv
	exec diff tmp/temp.tsv data/expected-vars2-var_annot3.tsv
} {}

test var_annot {multi alt, one value in vardb} {
	file mkdir tmp
	cg select -f {chromosome begin end type ref alt} data/vars1.sft tmp/vars.sft
	write_tab tmp/vars.sft {
		chromosome begin end type ref alt
		chr1 4001 4002 snp A G,C
	}
	write_tab tmp/var_annot.sft {
		chrom start end type ref alt name freq
		chr1 4001 4002 snp A G,C test2 0.2
	}
	file_write tmp/var_annot.sft.opt "fields\t{name freq alt}\n"
	exec cg annotate tmp/vars.sft tmp/temp.sft tmp/var_annot.sft
	diff_tab tmp/temp.sft {
		chromosome begin end type ref alt annot_name annot_freq annot_alt
		chr1 4001 4002 snp A G,C test2 0.2 G,C
	}
} {}

test var_annot {multi alt split, one value in vardb} {
	file mkdir tmp
	cg select -f {chromosome begin end type ref alt} data/vars1.sft tmp/vars.sft
	write_tab tmp/vars.sft {
		chromosome begin end type ref alt
		chr1 4001 4002 snp A C
		chr1 4001 4002 snp A G
	}
	write_tab tmp/var_annot.sft {
		chrom start end type ref alt name freq
		chr1 4001 4002 snp A G,C test2 0.2
	}
	file_write tmp/var_annot.sft.opt "fields\t{name freq alt}\n"
	exec cg annotate tmp/vars.sft tmp/temp.sft tmp/var_annot.sft
	diff_tab tmp/temp.sft {
		chromosome begin end type ref alt annot_name annot_freq annot_alt
		chr1 4001 4002 snp A C test2 0.2 C
		chr1 4001 4002 snp A G test2 0.2 G
	}
} {} 

test var_annot {multi alt split, one value in vardb (split)} {
	file mkdir tmp
	cg select -f {chromosome begin end type ref alt} data/vars1.sft tmp/vars.sft
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt
		chr1 4001 4002 snp A C
		chr1 4001 4002 snp A G
	}
	write_tab tmp/var_annot.tsv {
		chrom start end type ref alt name freq
		chr1 4001 4002 snp A G test2 0.2
	}
	write_tab tmp/expected.tsv {
		chromosome begin end type ref alt annot_name annot_freq annot_alt
		chr1 4001 4002 snp A C - - -
		chr1 4001 4002 snp A G test2 0.2 G
	}
	file_write tmp/var_annot.tsv.opt "fields\t{name freq alt}\n"
	exec cg annotate tmp/vars.tsv tmp/temp.tsv tmp/var_annot.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test var_annot {multi alt split, one value in vardb (split) no data} {
	test_cleantmp
	file mkdir tmp
	cg select -f {chromosome begin end type ref alt} data/vars1.sft tmp/vars.sft
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt
		chr1 4001 4002 snp A C
		chr1 4001 4002 snp A G
	}
	write_tab tmp/var_annot.tsv {
		chrom start end type ref alt
		chr1 4001 4002 snp A G
	}
	write_tab tmp/expected.tsv {
		chromosome begin end type ref alt annot
		chr1 4001 4002 snp A C -
		chr1 4001 4002 snp A G 1
	}
	exec cg annotate tmp/vars.tsv tmp/temp.tsv tmp/var_annot.tsv
	exec diff tmp/temp.tsv tmp/expected.tsv
} {}

test var_annot {sort error 1 in vars} {
	file copy data/vars_sorterror1.sft tmp/vars_sorterror1.sft
	exec cg annotate	tmp/vars_sorterror1.sft tmp/temp.sft data/var_annot.sft
} {*File (*) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test var_annot {sort error 2 in vars} {
	file copy data/vars_sorterror2.sft tmp/vars_sorterror2.sft
	exec cg annotate tmp/vars_sorterror2.sft tmp/temp.sft data/var_annot.sft
} {*File (*) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test var_annot {sort error 3 in vars} {
	file copy data/vars_sorterror3.sft tmp/vars_sorterror3.sft
	exec cg annotate tmp/vars_sorterror3.sft tmp/temp.sft data/var_annot.sft
} {*File (*) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test var_annot {sort error 4 in vars} {
	file copy data/vars_sorterror4.sft tmp/vars_sorterror4.sft
	exec cg annotate tmp/vars_sorterror4.sft tmp/temp.sft data/var_annot.sft
} {*File (*) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test var_annot {sort error 5 in vars} {
	file copy data/vars_sorterror5.sft tmp/vars_sorterror5.sft
	exec cg annotate tmp/vars_sorterror5.sft tmp/temp.sft data/var_annot.sft
} {*File (*) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test var_annot {sort error 1 in database file} {
	file copy -force data/vars_sorterror1.sft tmp/var_annot.sft
	file copy data/vars_annottest.sft tmp/vars_annottest.sft
	exec cg annotate tmp/vars_annottest.sft tmp/temp.sft tmp/var_annot.sft
} {*File (*) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test var_annot {sort error 2 in database file} {
	file copy -force data/vars_sorterror2.sft tmp/var_annot.sft
	file copy data/vars_annottest.sft tmp/vars_annottest.sft
	exec cg annotate tmp/vars_annottest.sft tmp/temp.sft tmp/var_annot.sft
} {*File (*) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test var_annot {sort error 3 in database file} {
	file copy -force data/vars_sorterror3.sft tmp/var_annot.sft
	file copy data/vars_annottest.sft tmp/vars_annottest.sft
	exec cg annotate tmp/vars_annottest.sft tmp/temp.sft tmp/var_annot.sft
} {*File (*) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test var_annot {sort error 4 in database file} {
	file copy -force data/vars_sorterror4.sft tmp/var_annot.sft
	file copy data/vars_annottest.sft tmp/vars_annottest.sft
	exec cg annotate tmp/vars_annottest.sft tmp/temp.sft tmp/var_annot.sft
} {*File (*) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test var_annot {sort error 5 in database file} {
	file copy -force data/vars_sorterror5.sft tmp/var_annot.sft
	file copy data/vars_annottest.sft tmp/vars_annottest.sft
	exec cg annotate tmp/vars_annottest.sft tmp/temp.sft tmp/var_annot.sft
} {*File (*) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test var_annot {lz4} {
	file copy data/vars1.sft tmp/vars1.sft
	exec lz4c -q -c data/var_annot.sft > tmp/var_annot.sft.lz4
	exec cg annotate tmp/vars1.sft tmp/temp.sft tmp/var_annot.sft.lz4
	exec cg select -rf {list} tmp/temp.sft tmp/temp2.sft
	exec diff tmp/temp2.sft data/expected-vars1-var_annot.sft
} {} 

test var_annot {skip var_ annots if no alt field} {
	file mkdir tmp
	write_tab tmp/vars.tsv {
		chromosome begin end type ref
		chr1 4001 4002 snp A
	}
	write_tab tmp/var_annot.tsv {
		chrom start end type ref alt name freq
		chr1 4001 4002 snp A G,C test2 0.2
	}
	write_tab tmp/reg_annot.tsv {
		chrom start end name
		chr1 4000 5000 test
	}
	exec cg annotate tmp/vars.tsv tmp/temp.tsv tmp/var_annot.tsv tmp/reg_annot.tsv
} {Skipping: */vars.tsv has no alt field} match

test var_annot {skip var_ annots if no alt field, check other annot} {
	file mkdir tmp
	write_tab tmp/vars.tsv {
		chromosome begin end type ref
		chr1 4001 4002 snp A
	}
	write_tab tmp/var_annotv.tsv {
		chrom start end type ref alt name freq
		chr1 4001 4002 snp A G,C test2 0.2
	}
	write_tab tmp/reg_annotr.tsv {
		chrom start end name
		chr1 4000 5000 test
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	annotr
		chr1	4001	4002	snp	A	test
	}
	exec cg annotate tmp/vars.tsv tmp/annot_vars.tsv tmp/var_annotv.tsv tmp/reg_annotr.tsv
	exec diff tmp/annot_vars.tsv tmp/expected.tsv
} {}

test gene_annot {variant file sort error 1} {
	file copy data/vars_sorterror1.sft tmp/vars_sorterror1.sft
	exec cg annotate -dbdir $::refseqdir/hg18 tmp/vars_sorterror1.sft tmp/temp.sft data/gene_test.tsv
} {*Cannot annotate because the variant file is not correctly sorted (sort correctly using "cg select -s -")*} error match

test gene_annot {variant file sort error 2} {
	file copy data/vars_sorterror2.sft tmp/vars_sorterror2.sft
	exec cg annotate -dbdir $::refseqdir/hg18 tmp/vars_sorterror2.sft tmp/temp.sft data/gene_test.tsv
} {*Cannot annotate because the variant file is not correctly sorted (sort correctly using "cg select -s -")*} error match

test gene_annot {variant file sort error 3} {
	file copy data/vars_sorterror3.sft tmp/vars_sorterror3.sft
	exec cg annotate -dbdir $::refseqdir/hg18 tmp/vars_sorterror3.sft tmp/temp.sft data/gene_test.tsv
} {*Cannot annotate because the variant file is not correctly sorted (sort correctly using "cg select -s -")*} error match

test gene_annot {gene wrongly sorted database file error} {
	cg select -s - data/vars_annottest.sft tmp/vars_annottest.sft
	exec cg annotate -dbdir $::refseqdir/hg18 tmp/vars_annottest.sft tmp/temp.sft data/gene_test-wrong1.tsv
} {*Cannot annotate because the database file (data/gene_test-wrong1.tsv) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test gene_annot {gene wrongly sorted database file error} {
	cg select -s - data/vars_annottest.sft tmp/vars_annottest.sft
	exec cg annotate -dbdir $::refseqdir/hg18 tmp/vars_annottest.sft tmp/temp.sft data/gene_test-wrong2.tsv
} {*Cannot annotate because the database file (data/gene_test-wrong2.tsv) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test gene_annot {gene} {
	cg select -s - data/vars_annottest.sft tmp/vars_annottest.sft
	exec cg annotate -dbdir $::refseqdir/hg18 tmp/vars_annottest.sft tmp/temp.sft data/gene_test.tsv
	catch {exec diff tmp/temp.sft data/expected-annotate-vars_annottest-gene_test.tsv} result
	set result
} {}

test gene_annot {gene -distrreg 1} {
	cg select -s - data/vars_annottest.sft tmp/vars_annottest.sft
	exec cg annotate -distrreg 1 -dbdir $::refseqdir/hg18 tmp/vars_annottest.sft tmp/temp.sft data/gene_test.tsv
	catch {exec diff tmp/temp.sft data/expected-annotate-vars_annottest-gene_test.tsv} result
	set result
} {}

test gene_annot {gene --upstreamsize option} {
	cg select -s - data/vars_annottest.sft tmp/vars_annottest.sft
	exec cg annotate --upstreamsize 1000 -dbdir $::refseqdir/hg18 tmp/vars_annottest.sft tmp/temp.sft data/gene_test.tsv
	exec diff tmp/temp.sft data/expected-annotate-vars_annottest-gene_test.tsv
} {44c44
< chr1	43198434	43198435	snp	T	G	"upstream SLC2A1"			
---
> chr1	43198434	43198435	snp	T	G	"upstream SLC2A1"	upstream	SLC2A1	-NM_006516:up-1001:c.-1526A>C
child process exited abnormally} error

test gene_annot {bug check empty _gene field with only name (used for transcript and gene)} {
	cg select -s - data/vars_annottest.sft tmp/vars_annottest.sft
	cg select -f {chromosome start end strand name cdsStart cdsEnd exonCount exonStarts exonEnds} data/gene_test.tsv tmp/gene_test.tsv
	exec cg annotate -dbdir $::refseqdir/hg18 tmp/vars_annottest.sft tmp/temp.sft tmp/gene_test.tsv
	lindex [cg select -g all -q {$test_gene ne ""} tmp/temp.sft] end
} {46} 

test gene_annot {gene exon deletion} {
	write_tab tmp/vars.tsv {
		chromosome	begin	end	type	ref	alt
		chr1	2499	2601	del	102	{}
	}
	write_tab tmp/gene_test.tsv {
		chrom	start	end	name	strand	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	name2
		chr1	2000	3000	test	+	2050	2950	3	2000,2500,2900,	2100,2600,3000,	testgene
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	test_impact	test_gene	test_descr
		chr1	2499	2601	del	102	{}	CDSSPLICE	testgene	+test:intron1+400_intron2+1:c.51-1_150+1del:p.?
	}
	exec cg annotate -dbdir $::refseqdir/hg18 tmp/vars.tsv tmp/result.tsv tmp/gene_test.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {} 

test gene_annot {gene exon deletion with no type given} {
	write_tab tmp/vars.tsv {
		chromosome	begin	end
		chr1	2499	2601
	}
	write_tab tmp/gene_test.tsv {
		chrom	start	end	name	strand	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	name2
		chr1	2000	3000	test	+	2050	2950	3	2000,2500,2900,	2100,2600,3000,	testgene
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	test_impact	test_gene	test_descr
		chr1	2499	2601	CDSSPLICE	testgene	+test:intron1+400_intron2+1:c.51-1_150+1del:p.?
	}
	exec cg annotate -dbdir $::refseqdir/hg18 tmp/vars.tsv tmp/result.tsv tmp/gene_test.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {} 

test gene_annot {gene and coding gene deletion} {
	write_tab tmp/vars.tsv {
		chromosome	begin	end	type	ref	alt
		chr1	1000	2000	del	1000	{}
	}
	write_tab tmp/gene_test.tsv {
		chrom	start	end	name	strand	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	name2
		chr1	1500	1800	test	+	{}	{}	1	1500,	1800,	testgene
		chr1	1500	1800	cdstest	+	1500	1800	1	1500,	1800,	cdstestgene
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	test_impact	test_gene	test_descr
		chr1	1000	2000	del	1000	{}	GENEDEL;GENEDEL	testgene;cdstestgene	testgene:del;cdstestgene:del
	}
	exec cg annotate -dbdir $::refseqdir/hg18 tmp/vars.tsv tmp/result.tsv tmp/gene_test.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {} 

test gene_annot {gene, extra comments} {
	file_write tmp/temp2.sft "# a comment\n# another comment\n"
	exec cg select -s - data/vars_annottest.sft >> tmp/temp2.sft
	exec cg annotate -dbdir $::refseqdir/hg18 tmp/temp2.sft tmp/temp.sft data/gene_test.tsv
	exec diff tmp/temp.sft data/expected-annotate-vars_annottest-gene_test.tsv
} {1,2d0
< # a comment
< # another comment
child process exited abnormally} error

test gene_annot {wrong nr fields} {
	# intentionally missing a column
	write_tab tmp/vars.tsv {
		chromosome	begin	end	type	ref	alt	comment
		chr1	851164	851165	snp	G	C
	}
	cg annotate -dbdir $::refseqdir/hg18 tmp/vars.tsv tmp/annot_results.tsv $::refseqdir/hg18/gene_hg18_refGene.tsv.lz4
	cg select -sh /dev/null -q {$refGene_impact eq "UTR5"} tmp/annot_results.tsv
} {chr1	851164	851165	snp	G	C		UTR5	SAMD11	+NM_152486:exon2+1:c.-20G>C}

test gene_annot {hgvs + strand gene coding} {
	# extra exon added to gene to test UT3 splice, short intron with only splice
	set dbline {chr1 850983 869932 NM_152486e + 591 851184 869396 14 850983,851164,855397,856281,861014,864282,864517,866386,867378,867652,867801,868495,868940,869150,869832, 851043,851256,855579,856332,861139,864372,864703,866549,867494,867731,868301,868620,869051,869824,869932, 0 SAMD11 cmpl cmpl -1,0,0,2,2,1,1,1,2,1,2,1,0,0,0,}
	file_write tmp/gene_part_test.tsv [join {chromosome start end name strand bin cdsStart cdsEnd exonCount exonStarts exonEnds id name2 cdsStartStat cdsEndStat exonFrames} \t]\n[join $dbline \t]\n
	if 0 {
		# just here for testing/debugging
		set genomefile $::refseqdir/hg18/genome_hg18.ifas
		catch {genome_close $genomef} ; set genomef [genome_open $genomefile]
		set dposs {0 1 2 4 6 7 8 9 10 3 12}
		set upstreamsize 2000
		set geneobj [annotatevar_gene_makegeneobj $genomef $dbline $dposs $upstreamsize]
		unset -nocomplain adata ; array set adata $geneobj
		# join $adata(ftlist) \n
	}
	cg select -s - data/annot_gene_tests_fw_coding.tsv tmp/sannot_gene_tests.tsv
	cg annotate -dbdir $::refseqdir/hg18 tmp/sannot_gene_tests.tsv tmp/annot_results.tsv tmp/gene_part_test.tsv
	set errors {}
	foreach line [split [cg select -sh /dev/null -q {$test_impact ne $expected_impact or $test_descr ne $expected_descr} tmp/annot_results.tsv] \n] {
		set line [split $line \t]
		append errors "[list set loc [lrange $line 0 5]]\n  r: [lindex $line 9] [lindex $line 11]\n  e: [lrange $line 6 7] \n"
	}	
	set errors
} {}

test gene_annot {hgvs + strand gene non-coding} {
	# extra exon added to gene to test UT3 splice, short intron with only splice
	write_tab tmp/gene_part_test.tsv {
		chromosome start end name strand bin cdsStart cdsEnd exonCount exonStarts exonEnds id name2 cdsStartStat cdsEndStat exonFrames
		chr1 850983 869932 NM_152486n + 591 851184 851184 14 850983,851164,855397,856281,861014,864282,864517,866386,867378,867652,867801,868495,868940,869150,869832, 851043,851256,855579,856332,861139,864372,864703,866549,867494,867731,868301,868620,869051,869824,869932, 0 SAMD11 cmpl cmpl -1,0,0,2,2,1,1,1,2,1,2,1,0,0,0,
	}
	cg select -s - data/annot_gene_tests_fw_noncoding.tsv tmp/sannot_gene_tests.tsv
	file delete tmp/annot_results.tsv
	cg annotate -dbdir $::refseqdir/hg18 tmp/sannot_gene_tests.tsv tmp/annot_results.tsv tmp/gene_part_test.tsv
	set errors {}
	foreach line [split [cg select -sh /dev/null -q {$test_impact ne $expected_impact or $test_descr ne $expected_descr} tmp/annot_results.tsv] \n] {
		set line [split $line \t]
		append errors "[list set loc [lrange $line 0 5]]\n  r: [lindex $line 9] [lindex $line 11]\n  e: [lrange $line 6 7] \n"
	}	
	set errors
} {}

test gene_annot {hgvs - strand gene coding} {
	# extra exon added to gene to test UT3 splice, short intron with only splice
	write_tab tmp/gene_part_test.tsv {
		chromosome start end name strand bin cdsStart cdsEnd exonCount exonStarts exonEnds id name2 cdsStartStat cdsEndStat exonFrames
		chr1	1706588	1812355	NM_002074	-	598	1708629	1746752	12	1706588,1708620,1710351,1711693,1714543,1725717,1727773,1737054,1739135,1746695,1760488,1812118,	1708352,1708736,1710568,1711895,1714610,1725880,1727837,1737161,1739174,1746798,1760537,1812355,	0	GNB1	cmpl	cmpl	-1,1,0,2,1,0,2,0,0,0,-1,-1,
	}
	cg select -s - data/annot_gene_tests_rv_coding.tsv tmp/sannot_gene_tests.tsv
	cg annotate -dbdir $::refseqdir/hg18 tmp/sannot_gene_tests.tsv tmp/annot_results.tsv tmp/gene_part_test.tsv
	set errors {}
	foreach line [split [cg select -sh /dev/null -q {$test_impact ne $expected_impact or $test_descr ne $expected_descr} tmp/annot_results.tsv] \n] {
		set line [split $line \t]
		append errors "[list set loc [lrange $line 0 5]]\n  r: [lindex $line 9] [lindex $line 11]\n  e: [lrange $line 6 7] \n"
	}	
	set errors
} {}

test gene_annot {hgvs - strand gene non-coding} {
	# extra exon added to gene to test UT3 splice, short intron with only splice
	set dbline {chr1	1706588	1812355	NM_002074	-	598	1706588	1706588	12	1706588,1708620,1710351,1711693,1714543,1725717,1727773,1737054,1739135,1746695,1760488,1812118,	1708352,1708736,1710568,1711895,1714610,1725880,1727837,1737161,1739174,1746798,1760537,1812355,	0	GNB1	cmpl	cmpl	-1,1,0,2,1,0,2,0,0,0,-1,-1,}
	file_write tmp/gene_part_test.tsv [join {chromosome start end name strand bin cdsStart cdsEnd exonCount exonStarts exonEnds id name2 cdsStartStat cdsEndStat exonFrames} \t]\n[join $dbline \t]\n
	cg select -s - data/annot_gene_tests_rv_noncoding.tsv tmp/sannot_gene_tests.tsv
	cg annotate -dbdir $::refseqdir/hg18 tmp/sannot_gene_tests.tsv tmp/annot_results.tsv tmp/gene_part_test.tsv
	set errors {}
	foreach line [split [cg select -sh /dev/null -q {$test_impact ne $expected_impact or $test_descr ne $expected_descr} tmp/annot_results.tsv] \n] {
		set line [split $line \t]
		append errors "[list set loc [lrange $line 0 5]]\n  r: [lindex $line 9] [lindex $line 11]\n  e: [lrange $line 6 7] \n"
	}	
	set errors
} {}

test gene_annot {hgvs - extra tests} {
	# extra exon added to gene to test UT3 splice, short intron with only splice
	write_tab tmp/gene_hg19s_part.tsv {
		chrom	start	end	strand	geneid	name	score	bin	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	name2	cdsStartStat	cdsEndStat	exonFrames
		chr9    35056064        35072739        -       VCP     NM_007126       0       852     35057113        35072350        17      35056064,35057372,35059060,35059489,35060309,35060797,35061011,35061573,35061999,35062213,35062974,35064150,35065247,35066671,35067887,35068247,35072333,     35057219,35057527,35059216,35059798,35060522,35060920,35061176,35061686,35062135,35062347,35063077,35064282,35065378,35066814,35068060,35068359,35072739,      VCP     cmpl    cmpl    2,0,0,0,0,0,0,1,0,1,0,0,1,2,0,2,0,
	}
	write_tab tmp/vars.tsv {
		chromosome	begin	end	type	ref	alt	comment
		9	35065348	35065349	snp	G	A	rs387906789
		9	35065359	35065360	snp	C	T	rs121909329
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	comment	part_impact	part_gene	part_descr
		9	35065348	35065349	snp	G	A	rs387906789	CDSMIS	VCP	-NM_007126:exon5+30:c.475C>T:p.R159C
		9	35065359	35065360	snp	C	T	rs121909329	CDSMIS	VCP	-NM_007126:exon5+19:c.464G>A:p.R155H
	}
	cg annotate -dbdir $::refseqdir/hg19 tmp/vars.tsv tmp/annot_results.tsv tmp/gene_hg19s_part.tsv
	exec diff tmp/annot_results.tsv tmp/expected.tsv
} {}

test gene_annot {variant error end > end chromosome} {
	write_tab tmp/vars.tsv {
		chromosome	begin	end	type	ref	alt
		Un_gl000220	117831	161803	del	43972	{}
	}
	file_write tmp/gene_test.tsv [deindent {
		chrom	start	end	strand	geneid	source	name	score	bin	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	name2	cdsStartStat	cdsEndStat	exonFrames	proteinID	alignID
		chrUn_gl000220	733	3190	-	FP236315.2	gene_hg19_gencode	ENST00000624506.1_1	0	585	733	733	9	733,781,823,985,1039,1515,1534,1833,1902,	780,822,984,1038,1514,1533,1832,1901,3190,	FP236315.2	none	none	-1,-1,-1,-1,-1,-1,-1,-1,-1,		
		chrUn_gl000220	733	3190	-	CU634019.3	gene_hg19_gencode	ENST00000623962.1_1	0	585	733	733	5	733,823,1039,1534,1902,	780,984,1514,1832,3190,	CU634019.3	none	none	-1,-1,-1,-1,-1,		
		chrUn_gl000220	2897	4713	+	FP236315.3	gene_hg19_gencode	ENST00000623587.1_1	0	585	2897	2897	4	2897,4515,4543,4703,	3006,4540,4697,4713,	FP236315.3	none	none	-1,-1,-1,-1,		
		chrUn_gl000220	4556	7010	-	FP236315.1	gene_hg19_gencode	ENST00000623723.1_1	0	585	4556	4556	3	4556,4703,6707,	4697,4710,7010,	FP236315.1	none	none	-1,-1,-1,		
		chrUn_gl000220	97128	126250	+	FP671120.1	gene_hg19_gencode	ENST00000623664.1_1	0	585	97128	97128	8	97128,100238,100409,121604,121828,122990,124777,125783,	97436,100398,100554,121701,122829,123082,124995,126250,	FP671120.1	none	none	-1,-1,-1,-1,-1,-1,-1,-1,		
		chrUn_gl000220	97128	126696	+	LOC100507412	known	uc011mfs.2			97128	97128	7	97128,100238,121604,121828,122990,124777,125783,	97436,100554,121701,122829,123082,124995,126696,						uc011mfs.2
		chrUn_gl000220	97128	126696	+	LOC100507412	ref	NR_038958	0	585	126696	126696	7	97128,100238,121604,121828,122990,124777,125783,	97436,100554,121701,122829,123082,124995,126696,	LOC100507412	unk	unk	-1,-1,-1,-1,-1,-1,-1,		
		chrUn_gl000220	100713	100791	+	AL592188.3	ens	ENST00000585135	0	585	100791	100791	1	100713,	100791,	ENSG00000263806	none	none	-1,		
		chrUn_gl000220	104731	104823	+	MIR6724-1	ref	NR_106782	0	585	104823	104823	1	104731,	104823,	MIR6724-1	unk	unk	-1,		
		chrUn_gl000220	105423	118780	+	RNA45SN5	ref	NR_046235	0	585	118780	118780	1	105423,	118780,	RNA45SN5	unk	unk	-1,		
		chrUn_gl000220	105423	156152	+	RNA5-8S5	known	uc022brd.2			105423	105423	2	105423,155996,	118780,156152,						uc022brd.2
		chrUn_gl000220	107447	107561	+	AL592188.2	ens	ENST00000488948	0	585	107561	107561	1	107447,	107561,	ENSG00000244180	none	none	-1,		
		chrUn_gl000220	107908	108088	+	AL592188.4	ens	ENST00000584459	0	585	108088	108088	1	107908,	108088,	ENSG00000264827	none	none	-1,		
		chrUn_gl000220	108279	108340	+	AL592188.8	ens	ENST00000582791	0	585	108340	108340	1	108279,	108340,	ENSG00000266219	none	none	-1,		
		chrUn_gl000220	109077	110946	+	RNA18SN5	ref	NR_003286	0	585	110946	110946	1	109077,	110946,	RNA18SN5	unk	unk	-1,		
		chrUn_gl000220	109830	110753	-	FP671120.3	gene_hg19_gencode	ENST00000631211.1_1	0	585	109830	109830	1	109830,	110753,	FP671120.3	none	none	-1,		
		chrUn_gl000220	112023	112180	+	RNA5-8SN5	ref	NR_003285	0	585	112180	112180	1	112023,	112180,	RNA5-8SN5	unk	unk	-1,		
		chrUn_gl000220	112024	112177	+	RNA5-8S5	ens	ENST00000474885	0	585	112177	112177	1	112024,	112177,	ENSG00000242716	none	none	-1,		
		chrUn_gl000220	112024	112180	+	RNA5-8SN4	ref	NR_146120	0	585	112180	112180	1	112024,	112180,	RNA5-8SN4	unk	unk	-1,		
		chrUn_gl000220	113347	118417	+	RNA28SN5	ref	NR_003287	0	585	118417	118417	1	113347,	118417,	RNA28SN5	unk	unk	-1,		
		chrUn_gl000220	114150	114242	+	RNA28S5	ens	ENST00000579027	0	585	114242	114242	1	114150,	114242,	ENSG00000266658	none	none	-1,		
		chrUn_gl000220	118196	118253	+	RNA28S5	ens	ENST00000582153	0	585	118253	118253	1	118196,	118253,	ENSG00000266658	none	none	-1,		
		chrUn_gl000220	125636	125756	+	JN872559	known	uc031tgc.1			125636	125636	1	125636,	125756,						uc031tgc.1
		chrUn_gl000220	126534	126717	+	JN872559	known	uc031tgd.1			126534	126534	1	126534,	126717,						uc031tgd.1
		chrUn_gl000220	133777	134256	+	JN872556	known	uc031tge.1			133777	133777	1	133777,	134256,						uc031tge.1
		chrUn_gl000220	148703	148795	+	MIR6724-1	ref	NR_106782	0	586	148795	148795	1	148703,	148795,	MIR6724-1	unk	unk	-1,		
		chrUn_gl000220	150505	150525	+	JA668106	known	uc022brf.1			150505	150505	1	150505,	150525,						uc022brf.1
		chrUn_gl000220	151419	151533	+	AL592188.1	ens	ENST00000459341	0	586	151533	151533	1	151419,	151533,	ENSG00000243151	none	none	-1,		
		chrUn_gl000220	151444	151466	+	JB158072	known	uc022brg.1			151444	151444	1	151444,	151466,						uc022brg.1
		chrUn_gl000220	151508	151530	+	JB074932	known	uc022brh.1			151508	151508	1	151508,	151530,						uc022brh.1
		chrUn_gl000220	151880	152060	+	AL592188.6	ens	ENST00000581046	0	586	152060	152060	1	151880,	152060,	ENSG00000265807	none	none	-1,		
		chrUn_gl000220	152251	152312	+	AL592188.5	ens	ENST00000577630	0	586	152312	152312	1	152251,	152312,	ENSG00000265525	none	none	-1,		
		chrUn_gl000220	153049	154918	+	RNA18SN5	ref	NR_003286	0	586	154918	154918	1	153049,	154918,	RNA18SN5	unk	unk	-1,		
		chrUn_gl000220	153802	154725	-	FP236383.2	gene_hg19_gencode	ENST00000625598.1_1	0	586	153802	153802	1	153802,	154725,	FP236383.2	none	none	-1,		
		chrUn_gl000220	153802	154725	-	FP671120.4	gene_hg19_gencode	ENST00000629969.1_1	0	586	153802	153802	2	153802,154623,	154617,154725,	FP671120.4	none	none	-1,-1,		
		chrUn_gl000220	155995	156152	+	RNA5-8SN5	ref	NR_003285	0	586	156152	156152	1	155995,	156152,	RNA5-8SN5	unk	unk	-1,		
		chrUn_gl000220	155996	156149	+	RNA5-8S5	ens	ENST00000459522	0	586	156149	156149	1	155996,	156149,	ENSG00000241335	none	none	-1,		
		chrUn_gl000220	155996	156152	+	RNA5-8SN4	ref	NR_146120	0	586	156152	156152	1	155996,	156152,	RNA5-8SN4	unk	unk	-1,		
		chrUn_gl000220	158122	158214	+	AL592188.7	ens	ENST00000581141	0	586	158214	158214	1	158122,	158214,	ENSG00000265830	none	none	-1,		
	}]
	file_write tmp/expected.tsv [deindent {
		chromosome	begin	end	type	ref	alt	test_impact	test_gene	test_descr
		Un_gl000220	117831	161803	del	43972		RNAEND;RNAEND;RNAEND;RNAEND;RNAEND;GENEDEL;GENEDEL;GENEDEL;GENEDEL;GENEDEL;GENEDEL;GENEDEL;GENEDEL;GENEDEL;GENEDEL;GENEDEL;GENEDEL;GENEDEL;GENEDEL;GENEDEL;GENEDEL	FP671120.1;LOC100507412;RNA45SN5;RNA5-8S5;RNA28SN5;RNA28S5;JN872559;JN872556;MIR6724-1;JA668106;AL592188.1;JB158072;JB074932;AL592188.6;AL592188.5;RNA18SN5;FP236383.2;FP671120.4;RNA5-8SN5;RNA5-8S5;RNA5-8SN4	+ENST00000623664.1_1:intron3+17278_down+35553:n.614-3773_38041del;+uc011mfs.2:intron2+17278_down+35107:n.625-3773_38052del;+NR_046235:exon1+12409_down+43023:n.12409_56380del;+uc022brd.2:exon1+12409_down+5651:n.12409_19164del;+NR_003287:exon1+4485_down+43386:n.4485_48456del;RNA28S5:del;JN872559:del;JN872556:del;MIR6724-1:del;JA668106:del;AL592188.1:del;JB158072:del;JB074932:del;AL592188.6:del;AL592188.5:del;RNA18SN5:del;FP236383.2:del;FP671120.4:del;RNA5-8SN5:del;RNA5-8S5:del;RNA5-8SN4:del
	}]\n
	exec cg annotate -dbdir $::refseqdir/hg19 tmp/vars.tsv tmp/result.tsv tmp/gene_test.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {} 

test reg_annot {basic, extra comments} {
	file_write tmp/temp2.sft "# a comment\n"
	exec cat data/vars1.sft >> tmp/temp2.sft
	exec cg annotate tmp/temp2.sft tmp/temp.sft data/reg_annot.sft
	exec cg select -rf {list} tmp/temp.sft tmp/temp3.sft
	exec diff tmp/temp3.sft data/expected-vars1-reg_annot.sft
} {1d0
< # a comment
child process exited abnormally} error

test var_annot {different types on same pos, extra comments} {
	file_write tmp/temp2.sft "# a comment\n"
	exec cat data/vars2.tsv >> tmp/temp2.sft
	exec cg annotate tmp/temp2.sft tmp/temp.tsv data/var_annot3.tsv
	exec diff tmp/temp.tsv data/expected-vars2-var_annot3.tsv
} {1d0
< # a comment
child process exited abnormally} error

test reg_annot {existing field error} {
	file copy data/vars1.sft tmp/vars1.sft
	exec cg annotate tmp/vars1.sft tmp/temp.sft data/reg_annot.sft
	exec cg annotate -near 1000 tmp/temp.sft tmp/temp2.sft data/reg_annot.sft
} {*Error: field(s) regtest already in file} match error

test reg_annot {-replace y dbfile newer} {
	file copy -force data/vars1.sft tmp/vars1.sft
	file copy -force data/reg_annot.sft tmp/reg_annot.tsv
	file copy -force data/reg_annot.sft.opt tmp/reg_annot.tsv.opt
	exec cg annotate tmp/vars1.sft tmp/temp.sft tmp/reg_annot.tsv
	after 1000
	exec touch tmp/reg_annot.tsv
	exec cg annotate -near 1000 -replace y tmp/temp.sft tmp/temp2.sft tmp/reg_annot.tsv
	exec cg select -rf {list} tmp/temp2.sft tmp/temp3.sft
	exec diff tmp/temp3.sft data/expected_near-vars1-reg_annot.sft
	split [cg select -f regtest tmp/temp2.sft] \n
} {regtest reg1 reg1 reg1 reg1 reg1 reg3 reg3 reg3 reg3 reg3 reg3 reg3 reg3 {}}

test reg_annot {-replace y dbfile older} {
	file copy -force data/vars1.sft tmp/vars1.sft
	file copy -force data/reg_annot.sft tmp/reg_annot.tsv
	file copy -force data/reg_annot.sft.opt tmp/reg_annot.tsv.opt
	exec cg annotate tmp/vars1.sft tmp/temp.sft tmp/reg_annot.tsv
	after 100
	exec touch tmp/temp.sft
	exec cg annotate -near 1000 -replace y tmp/temp.sft tmp/temp2.sft tmp/reg_annot.tsv
	exec cg select -rf {list} tmp/temp2.sft tmp/temp3.sft
	split [cg select -f regtest tmp/temp3.sft] \n
} {regtest reg1 reg1 reg1 reg1 {} {} {} {} reg3 reg3 {} {} {} {}} 

test reg_annot {-replace a} {
	file copy -force data/vars1.sft tmp/vars1.sft
	file copy -force data/reg_annot.sft tmp/reg_annot.tsv
	file copy -force data/reg_annot.sft.opt tmp/reg_annot.tsv.opt
	exec cg annotate tmp/vars1.sft tmp/temp.sft tmp/reg_annot.tsv
	after 100
	exec touch tmp/temp.sft
	exec cg annotate -near 1000 -replace a tmp/temp.sft tmp/temp2.sft tmp/reg_annot.tsv
	exec cg select -rf {list} tmp/temp2.sft tmp/temp3.sft
	exec diff tmp/temp3.sft data/expected_near-vars1-reg_annot.sft
} {} 

test reg_annot {-replace n} {
	file copy data/vars1.sft tmp/vars1.sft
	exec cg annotate tmp/vars1.sft tmp/temp.sft data/reg_annot.sft
	exec cg annotate -near 1000 -replace n tmp/temp.sft tmp/temp2.sft data/reg_annot.sft
	exec cg select -rf {list} tmp/temp2.sft tmp/temp3.sft
	split [cg select -f regtest tmp/temp3.sft] \n
} {regtest reg1 reg1 reg1 reg1 {} {} {} {} reg3 reg3 {} {} {} {}} 

test reg_annot {-replace e} {
	file copy data/vars1.sft tmp/vars1.sft
	exec cg annotate tmp/vars1.sft tmp/temp.sft data/reg_annot.sft
	exec cg annotate -near 1000 -replace e tmp/temp.sft tmp/temp2.sft data/reg_annot.sft
	exec cg select -rf {list} tmp/temp2.sft tmp/temp3.sft
	split [cg select -f regtest tmp/temp3.sft] \n
} {Error: field(s) regtest already in file} error

test bcol_annot {basic} {
	test_cleantmp
	exec cg bcol make -p pos -c chromosome tmp/bcol_coverage.tsv coverage < data/cov.tsv
	file copy data/bcol_annot-test.tsv tmp/bcol_annot-test.tsv
	exec cg annotate tmp/bcol_annot-test.tsv tmp/annot_test.tsv tmp/bcol_coverage.tsv
	exec diff tmp/annot_test.tsv data/expected-bcol_annot-test.tsv
} {} 

test bcol_annot {header error} {
	test_cleantmp
	exec cg bcol make -p pos -c chromosome tmp/bcol_coverage.tsv coverage < data/cov.tsv
	file_write tmp/bcol_coverage.tsv "chr1\ttemp-chr1.bcol\nchr2\ttemp-chr2.bcol\n"
	file copy data/bcol_annot-test.tsv tmp/bcol_annot-test.tsv
	exec cg annotate tmp/bcol_annot-test.tsv tmp/annot_test.tsv tmp/bcol_coverage.tsv
} {*bcol database (tmp/bcol_coverage.tsv) should have a header of the type: chromosome begin end, old style bcols (single chr not supported)*} error match

test bcol_annot {basic compressed bcol} {
	test_cleantmp
	exec cg bcol make -p pos -c chromosome tmp/bcol_coverage.tsv coverage < data/cov.tsv
	file copy data/bcol_annot-test.tsv tmp/bcol_annot-test.tsv
	exec cg annotate tmp/bcol_annot-test.tsv tmp/annot_test.tsv tmp/bcol_coverage.tsv
	exec diff tmp/annot_test.tsv data/expected-bcol_annot-test.tsv
} {} 

test bcol_annot {basic uncompressed bcol} {
	test_cleantmp
	exec cg bcol make -p pos -c chromosome tmp/bcol_coverage.tsv coverage < data/cov.tsv
	cg unzip {*}[glob tmp/*.rz tmp/*.lz4]
	file copy data/bcol_annot-test.tsv tmp/bcol_annot-test.tsv
	exec cg annotate tmp/bcol_annot-test.tsv tmp/annot_test.tsv tmp/bcol_coverage.tsv
	exec diff tmp/annot_test.tsv data/expected-bcol_annot-test.tsv
} {} 

test bcol_annot {only chr1} {
	test_cleantmp
	cg select -q {$chromosome eq "chr1"} data/bcol_annot-test.tsv tmp/bcol_annot-test.tsv
	cg select -q {$chromosome eq "chr1"} data/expected-bcol_annot-test.tsv tmp/expected.tsv
	exec cg bcol make -p pos -c chromosome tmp/bcol_coverage.tsv coverage < data/cov.tsv
	exec cg annotate tmp/bcol_annot-test.tsv tmp/annot_test.tsv tmp/bcol_coverage.tsv
	exec diff tmp/annot_test.tsv tmp/expected.tsv
} {} 

test bcol_annot {wrongly sorted var file} {
	file_write tmp/test.tsv [deindent {
		chromosome	begin	end
		chr1	8	9
		chr10	9	10
		chr2	19	20
	}]\n
	file_write tmp/score.tsv [deindent {
		chromosome	pos	score
		chr1	8	1
		chr2	19	2
		chr10	9	3
	}]\n
	exec cg bcol make -p pos -c chromosome tmp/bcol_score.bcol score < tmp/score.tsv
	exec cg annotate tmp/test.tsv tmp/annot_test.tsv tmp/bcol_score.bcol
} {File (*/test.tsv.index/vars.tsv) is not correctly sorted (sort correctly using "cg select -s -")
chr10:9-10:: came before chr2:19-20::
*} error match

test bcol_annot {wrongly sorted bcol annot file} {
	file_write tmp/test.tsv [deindent {
		chromosome	begin	end
		chr1	8	9
		chr2	19	20
		chr10	9	10
	}]\n
	file_write tmp/score.tsv [deindent {
		chromosome	pos	score
		chr1	8	1
		chr10	9	3
		chr2	19	2
	}]\n
	exec cg bcol make -p pos -c chromosome tmp/bcol_score.bcol score < tmp/score.tsv
	exec cg annotate tmp/test.tsv tmp/annot_test.tsv tmp/bcol_score.bcol
} {File (*/bcol_score.bcol) is not correctly sorted:
chromosome 10 came before 2
*} error match

test bcol_annot {only chr2} {
	test_cleantmp
	cg select -q {$chromosome eq "chr2"} data/bcol_annot-test.tsv tmp/bcol_annot-test.tsv
	cg select -q {$chromosome eq "chr2"} data/expected-bcol_annot-test.tsv tmp/expected.tsv
	exec cg bcol make -p pos -c chromosome tmp/bcol_coverage.tsv coverage < data/cov.tsv
	exec cg annotate tmp/bcol_annot-test.tsv tmp/annot_test.tsv tmp/bcol_coverage.tsv
	exec diff tmp/annot_test.tsv tmp/expected.tsv
} {} 

test bcol_annot {added chr3} {
	test_cleantmp
	file copy data/bcol_annot-test.tsv tmp/bcol_annot-test.tsv
	set f [open tmp/bcol_annot-test.tsv a]
	puts $f "chr3\t9\t10"
	close $f
	exec cg bcol make -p pos -c chromosome tmp/bcol_coverage.tsv coverage < data/cov.tsv
	exec cg annotate tmp/bcol_annot-test.tsv tmp/annot_test.tsv tmp/bcol_coverage.tsv
	exec diff tmp/annot_test.tsv data/expected-bcol_annot-test.tsv
} {11d10
< chr3	9	10	0
child process exited abnormally} error 

test bcol_annot {--precision} {
	test_cleantmp
	write_tab tmp/annot.tsv {
		chromosome	begin	score
		chr1	0	0.421
	}
	cg bcol make --precision 2 -t f -p begin -c chromosome tmp/bcol_annot.bcol score < tmp/annot.tsv
	write_tab tmp/vars.tsv {
		chromosome	begin	end	type
		chr1	0	1	snp
	}
	exec cg annotate tmp/vars.tsv tmp/results.tsv tmp/bcol_annot.bcol
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	annot
		chr1	0	1	snp	0.42
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test var_annot {basic from vcf} {
	file copy data/vars1.vcf tmp/vars1.vcf
	exec cg annotate tmp/vars1.vcf tmp/annot.sft data/var_annot.sft
	set fields {chromosome	begin	end	type	ref	alt	alleleSeq1-sample1	alleleSeq2-sample1	coverage-sample1	alleleSeq1-sample2	alleleSeq2-sample2	coverage-sample2	annot_name	annot_freq}
	exec cg select -rc 1 -f $fields tmp/annot.sft tmp/annot2.sft
	cg splitalleles data/expected-vars1-var_annot.sft tmp/expected.sft.temp
	exec cg select -rc 1 -f $fields tmp/expected.sft.temp tmp/expected.sft
	exec diff tmp/annot2.sft tmp/expected.sft
} {} 

test reg_annot {ins at end of reg} {
	test_cleantmp
	write_tab tmp/vars.tsv {
		chromosome begin end	type	num
		1	9	10	snp	1
		1	10	10	ins	2
		1	10	11	snp	3
		1	11	11	ins	4
		1	19	20	snp	5
		1	20	20	ins	6
		1	20	21	snp	7
		1	21	21	ins	8
	}
	write_tab tmp/reg_test.tsv {
		chromosome	begin	end	name
		1	10	20	10-20
	}
	write_tab tmp/expected.tsv {
		chromosome begin end	type	num	test
		1	9	10	snp	1	{}
		1	10	10	ins	2	10-20
		1	10	11	snp	3	10-20
		1	11	11	ins	4	10-20
		1	19	20	snp	5	10-20
		1	20	20	ins	6	10-20
		1	20	21	snp	7	{}
		1	21	21	ins	8	{}
	}
	exec cg annotate tmp/vars.tsv tmp/result.tsv tmp/reg_test.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {} 

test reg_annot {comments} {
	test_cleantmp
	write_tab tmp/vars.tsv {
		#a	b
		chromosome begin end	type	num
		1	19	20	snp	5
	}
	write_tab tmp/reg_test.tsv {
		chromosome	begin	end	name
		1	10	20	10-20
	}
	write_tab tmp/expected.tsv {
		#a	b
		chromosome begin end	type	num	test
		1	19	20	snp	5	10-20
	}
	exec cg annotate tmp/vars.tsv tmp/result.tsv tmp/reg_test.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {} 

test bcol_var_annot {basic} {
	test_cleantmp
	cg bcol make -t f --multicol alt --multilist A,C,T,G -p begin -c chromosome tmp/var_annot.bcol score < data/var-annot.tsv
	write_tab tmp/vars.tsv {
		chromosome	begin	end	type	alt
		chr1	0	1	snp	A
		chr1	0	1	snp	C
		chr1	1	2	snp	A
		chr1	1	2	snp	C
		chr1	1	2	snp	T
		chr1	4	5	snp	T
		chr1	10	20	snp	T
		chr2	22	23	snp	G
		chr2	26	27	snp	A
		chr2	29	30	snp	C
	}
	exec cg annotate tmp/vars.tsv tmp/results.tsv tmp/var_annot.bcol
	# first = 0 because type f is rounded to 6 decimal places (precision)
	# latest = 3000000 because type f cannot store 3000000.1 
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	alt	annot
		chr1	0	1	snp	A	0
		chr1	0	1	snp	C	0.1
		chr1	1	2	snp	A	2
		chr1	1	2	snp	C	0
		chr1	1	2	snp	T	2.1
		chr1	4	5	snp	T	6.1
		chr1	10	20	snp	T	0
		chr2	22	23	snp	G	0
		chr2	26	27	snp	A	0
		chr2	29	30	snp	C	3000000
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test bcol_var_annot {basic default} {
	test_cleantmp
	cg bcol make -t f --multicol alt -d - --multilist A,C,T,G -p begin -c chromosome tmp/var_annot.bcol score < data/var-annot.tsv
	write_tab tmp/vars.tsv {
		chromosome	begin	end	type	alt
		chr1	0	1	snp	A
		chr1	0	1	snp	C
		chr1	1	2	snp	A
		chr1	1	2	snp	C
		chr1	1	2	snp	T
		chr1	4	5	snp	T
		chr1	10	20	snp	T
		chr2	22	23	snp	G
		chr2	26	27	snp	A
		chr2	29	30	snp	C
	}
	exec cg annotate tmp/vars.tsv tmp/results.tsv tmp/var_annot.bcol
	# first = 0 because type f is rounded to 6 decimal places (precision)
	# latest = 3000000 because type f cannot store 3000000.1 
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	alt	annot
		chr1	0	1	snp	A	0
		chr1	0	1	snp	C	0.1
		chr1	1	2	snp	A	2
		chr1	1	2	snp	C	0
		chr1	1	2	snp	T	2.1
		chr1	4	5	snp	T	6.1
		chr1	10	20	snp	T	0
		chr2	22	23	snp	G	0
		chr2	26	27	snp	A	0
		chr2	29	30	snp	C	3000000
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test bcol_var_annot {basic type d} {
	test_cleantmp
	cg bcol make -t d --multicol alt --multilist A,C,T,G -p begin -c chromosome tmp/var_annot.bcol score < data/var-annot.tsv
	write_tab tmp/vars.tsv {
		chromosome	begin	end	type	alt
		chr1	0	1	snp	A
		chr1	0	1	snp	C
		chr1	1	2	snp	A
		chr1	1	2	snp	C
		chr1	1	2	snp	T
		chr1	4	5	snp	T
		chr1	10	20	snp	T
		chr2	22	23	snp	G
		chr2	26	27	snp	A
		chr2	29	30	snp	C
	}
	exec cg annotate tmp/vars.tsv tmp/results.tsv tmp/var_annot.bcol
	# first is ok because type d is rounded to 9 decimal places (precision)
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	alt	annot
		chr1	0	1	snp	A	0.0000001
		chr1	0	1	snp	C	0.1
		chr1	1	2	snp	A	2
		chr1	1	2	snp	C	0
		chr1	1	2	snp	T	2.1
		chr1	4	5	snp	T	6.1
		chr1	10	20	snp	T	0
		chr2	22	23	snp	G	0
		chr2	26	27	snp	A	0
		chr2	29	30	snp	C	3000000.1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test bcol_var_annot {basic uncompressed} {
	test_cleantmp
	cg bcol make -t f --multicol alt --multilist A,C,T,G --compress 0 -p begin -c chromosome tmp/var_annot.bcol score < data/var-annot.tsv
	write_tab tmp/vars.tsv {
		chromosome	begin	end	type	alt
		chr1	0	1	snp	A
		chr1	1	2	snp	A
		chr1	1	2	snp	C
		chr1	1	2	snp	T
		chr1	4	5	snp	T
		chr1	10	20	snp	T
		chr1	10	20	sub	AT
		chr2	22	23	snp	G
		chr2	26	27	snp	A
		chr2	29	30	snp	C
	}
	exec cg annotate tmp/vars.tsv tmp/results.tsv tmp/var_annot.bcol
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	alt	annot
		chr1	0	1	snp	A	0
		chr1	1	2	snp	A	2
		chr1	1	2	snp	C	0
		chr1	1	2	snp	T	2.1
		chr1	4	5	snp	T	6.1
		chr1	10	20	snp	T	0
		chr1	10	20	sub	AT	0
		chr2	22	23	snp	G	0
		chr2	26	27	snp	A	0
		chr2	29	30	snp	C	3000000
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test bcol_var_annot {split 0} {
	test_cleantmp
	cg bcol make -t d --multicol alt --multilist A,C,T,G -p begin -c chromosome tmp/var_annot.bcol score < data/var-annot.tsv
	write_tab tmp/vars.tsv {
		chromosome	begin	end	type	alt
		chr1	1	2	snp	A,C,T
		chr1	4	5	snp	T
		chr1	10	20	snp	T
		chr2	22	23	snp	C,G
		chr2	26	27	snp	A
		chr2	29	30	snp	C,T
		chr3	600000000	600000001	snp	G
	}
	exec cg annotate tmp/vars.tsv tmp/results.tsv tmp/var_annot.bcol
	# 30.200001 because of float
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	alt	annot
		chr1	1	2	snp	A,C,T	2,0,2.1
		chr1	4	5	snp	T	6.1
		chr1	10	20	snp	T	0
		chr2	22	23	snp	C,G	0,0
		chr2	26	27	snp	A	0
		chr2	29	30	snp	C,T	3000000.1,30.21
		chr3	600000000	600000001	snp	G	0
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test bcol_var_annot {split 0, precision} {
	test_cleantmp
	write_tab tmp/var_annot.tsv {
		chromosome	begin	end	type	alt	score
		chr1	1	2	snp	A,C,T	0.1,0.01,0.001
		chr1	2	3	snp	A,C,T	0.0001,0.00001,0.000001
		chr1	4	5	snp	A,C,T	0.0000001,0.00000001,0.000000001
	}
	cg bcol make -t f --multicol alt --multilist A,C,T,G -p begin -c chromosome tmp/var_score.bcol score < tmp/var_annot.tsv
	cg select -rf score tmp/var_annot.tsv tmp/vars.tsv
	exec cg annotate tmp/vars.tsv tmp/results.tsv tmp/var_score.bcol
	exec diff tmp/results.tsv tmp/var_annot.tsv
} {4c4
< chr1	4	5	snp	A,C,T	0,0,0
---
> chr1	4	5	snp	A,C,T	0.0000001,0.00000001,0.000000001
child process exited abnormally} error

test bcol_var_annot {--precision} {
	test_cleantmp
	write_tab tmp/annot.tsv {
		chromosome	begin	end	type	alt	score
		chr1	0	1	snp	A,G	0.421,0.9
	}
	cg bcol make --precision 2 -t f --multicol alt --multilist A,C,T,G -p begin -c chromosome tmp/var_annot.bcol score < tmp/annot.tsv
	write_tab tmp/vars.tsv {
		chromosome	begin	end	type	alt
		chr1	0	1	snp	A,G
	}
	exec cg annotate tmp/vars.tsv tmp/results.tsv tmp/var_annot.bcol
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	alt	annot
		chr1	0	1	snp	A,G	0.42,0.9
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test gene_annot {multiple dbs} {
	write_tab tmp/vars.tsv {
		chromosome	begin	end	type	ref	alt
		chr1	1000	2000	del	1000	{}
	}
	write_tab tmp/gene_test.tsv {
		chrom	start	end	name	strand	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	name2
		chr1	1500	1800	test	+	{}	{}	1	1500,	1800,	testgene
		chr1	1500	1800	cdstest	+	1500	1800	1	1500,	1800,	cdstestgene
	}
	write_tab tmp/reg_rtest.tsv {
		chrom	start	end	name	score
		chr1	500	1000	test1	1
		chr1	1000	2000	test2	2
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	test_impact	test_gene	test_descr	rtest_name	rtest_score
		chr1	1000	2000	del	1000	{}	GENEDEL;GENEDEL	testgene;cdstestgene	testgene:del;cdstestgene:del	test2	2
	}
	exec cg annotate -dbdir $::refseqdir/hg18 tmp/vars.tsv tmp/result.tsv tmp/gene_test.tsv tmp/reg_rtest.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {} 

test reg_annot {bugcheck overwrite of .temp src} {
	file copy data/vars1.sft tmp/vars1.tsv.temp
	exec cg annotate tmp/vars1.tsv.temp tmp/vars1.tsv data/reg_annot.sft
	exec cg select -rf {list} tmp/vars1.tsv tmp/temp2.tsv
	exec diff tmp/temp2.tsv data/expected-vars1-reg_annot.sft
} {}

test gene_annot {no refseq} {
	write_tab tmp/vars.tsv {
		chromosome	begin	end	type	ref	alt
		test	1600	1601	del	N	{}
	}
	write_tab tmp/gene_test.tsv {
		chrom	start	end	name	strand	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	name2
		test	1500	1800	test	+	{}	{}	1	1500,	1800,	testgene
		test	1500	1800	cdstest	+	1500	1800	1	1500,	1800,	cdstestgene
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	test_impact	test_gene	test_descr
		test	1600	1601	del	N	{}	RNA;CDSFRAME	testgene;cdstestgene	+test:exon1+101:n.101del;+cdstest:exon1+101:c.101del:p.X101Xfs*?
	}
	exec cg annotate -dbdir $::refseqdir/hg18 tmp/vars.tsv tmp/result.tsv tmp/gene_test.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {} 

test reg_annot {check for diff size in paste error} {
	write_tab tmp/vars1.tsv {
		chromosome	begin	end	type	ref	alt
		chr1	4000	4001	snp	G	A
		chr2	4000	4001	snp	G	A
	}
	file mkdir tmp/vars1.tsv.index
	write_tab tmp/vars1.tsv.index/vars.tsv {
		chromosome	begin	end	type	ref	alt
		chr1	4000	4001	snp	G	A
	}
	write_tab tmp/reg_annot.tsv {
		chromosome	begin	end	name
		chr1	4000	4010	A
		chr2	4000	4010	B
	}
	exec cg annotate tmp/vars1.tsv tmp/temp.tsv tmp/reg_annot.tsv
} {*file */vars.tsv.annot_annot has less lines than other files in paste*} error match

test reg_annot {bug fix, hang chromosome not in reference (do not match Ns when moving del pos)} {
	write_tab tmp/vars.tsv {
		chromosome	begin	end	type	ref	alt
		6_qbl_hap6	1191644	1191645	del	T
	}
	write_tab tmp/gene_annot.tsv {
		chrom	start	end	strand	geneid	source	name	score	bin	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	name2	cdsStartStat	cdsEndStat	exonFrames	proteinID	alignID
		chr6_qbl_hap6	1150112	1203581	+	HLA-A	known	uc031suu.1	{}	{}	1150221	1203580	2	1150112,1203311,	1150294,1203581,	{}	{}	{}	{}	O19559	uc031suu.1
		chr6_qbl_hap6	1150221	1152943	+	HLA-H	ens	ENST00000429724	0	593	1152943	1152943	7	1150221,1150418,1150929,1151790,1152167,1152724,1152899,	1150294,1150688,1151205,1152065,1152284,1152757,1152943,	ENSG00000227296	none	none	-1,-1,-1,-1,-1,-1,-1,	{}	{}
		chr6_qbl_hap6	1151792	1206454	+	HLA-A	known	uc011igl.2	{}	{}	1151807	1206025	6	1151792,1204756,1205073,1205628,1205803,1206020,	1151848,1204974,1205190,1205661,1205851,1206454,	{}	{}	{}	{}	B4DJI3	uc011igl.2
		chr6_qbl_hap6	1159017	1160149	+	HLA-T	ens	ENST00000430243	0	593	1160149	1160149	4	1159017,1159414,1159925,1160110,	1159289,1159527,1159962,1160149,	ENSG00000229552	none	none	-1,-1,-1,-1,	{}	{}
		chr6_qbl_hap6	1161752	1203991	-	AK097625	known	uc011ign.1	{}	{}	1161752	1161752	6	1161752,1187266,1189768,1190537,1190824,1203805,	1163306,1187359,1189984,1190641,1190932,1203991,	{}	{}	{}	{}	{}	uc011ign.1
		chr6_qbl_hap6	1168977	1169344	+	DDX39BP1	ens	ENST00000439001	0	593	1169344	1169344	2	1168977,1169287,	1169112,1169344,	ENSG00000236441	none	none	-1,-1,	{}	{}
		chr6_qbl_hap6	1170219	1171079	-	MCCD1P1	ens	ENST00000429760	0	593	1171079	1171079	2	1170219,1170912,	1170377,1171079,	ENSG00000229411	none	none	-1,-1,	{}	{}
		chr6_qbl_hap6	1187015	1188075	-	HCG4B	known	uc011igo.1	{}	{}	1187015	1187015	1	1187015,	1188075,	{}	{}	{}	{}	{}	uc011igo.1
		chr6_qbl_hap6	1187015	1189600	-	HCG4B	ref	NR_001317	0	594	1189600	1189600	1	1187015,	1189600,	HCG4B	unk	unk	-1,	{}	{}
		chr6_qbl_hap6	1187146	1189600	-	HCG4B	ens	ENST00000550894	0	594	1189600	1189600	1	1187146,	1189600,	ENSG00000220391	none	none	-1,	{}	{}
		chr6_qbl_hap6	1188367	1189358	-	HCG4B	ens	ENST00000406611	0	594	1189358	1189358	1	1188367,	1189358,	ENSG00000220391	none	none	-1,	{}	{}
		chr6_qbl_hap6	1188713	1192114	+	BC035647	known	uc011igp.2	{}	{}	1188713	1188713	5	1188713,1190842,1191404,1191571,1191788,	1190713,1190959,1191431,1191619,1192114,	{}	{}	{}	{}	{}	uc011igp.2
		chr6_qbl_hap6	1188843	1191618	+	HLA-K	ens	ENST00000431463	0	594	1191618	1191618	6	1188843,1189044,1189562,1190434,1190842,1191571,	1188915,1189305,1189837,1190713,1190959,1191618,	ENSG00000224526	none	none	-1,-1,-1,-1,-1,-1,	{}	{}
		chr6_qbl_hap6	1194728	1194913	+	HLA-U	ens	ENST00000442622	0	594	1194913	1194913	1	1194728,	1194913,	ENSG00000229428	none	none	-1,	{}	{}
		chr6_qbl_hap6	1201813	1206454	+	HLA-A	ens	ENST00000383605	0	594	1203108	1206025	10	1201813,1202520,1202958,1203311,1203822,1204698,1205073,1205628,1205803,1206020,	1201873,1202651,1203181,1203581,1204098,1204974,1205190,1205661,1205851,1206454,	ENSG00000206505	cmpl	cmpl	-1,-1,0,1,1,1,1,1,1,1,	{}	{}
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	annot_impact	annot_gene	annot_descr
		6_qbl_hap6	1191644	1191645	del	T	{}	intron	HLA-A;HLA-A;AK097625;BC035647	+uc031suu.1:intron1+41351:c.74-11667del;+uc011igl.2:intron1+39797:c.42-13112del;-uc011ign.1:intron1+12161:n.187-713del;+uc011igp.2:intron4+26:n.2192+26del
	}
	exec cg annotate -dbdir $::refseqdir/hg19 tmp/vars.tsv tmp/annot_vars.tsv tmp/gene_annot.tsv
	exec diff tmp/annot_vars.tsv tmp/expected.tsv
} {}

file delete -force tmp/temp.sft
file delete -force tmp/temp2.sft

set ::env(PATH) $keeppath

testsummarize
