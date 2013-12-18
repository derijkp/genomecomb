#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

file mkdir tmp

test reg_annot {basic} {
	exec cg annotate data/vars1.sft tmp/temp.sft data/reg_annot.sft 2> /dev/null
	exec cg select -rf {list} tmp/temp.sft tmp/temp2.sft
	exec diff tmp/temp2.sft data/expected-vars1-reg_annot.sft
} {} 

test reg_annot {compressed} {
	file copy -force data/reg_annot.sft data/reg_annot.sft.opt tmp
	exec cg razip tmp/reg_annot.sft 2> /dev/null
	exec cg annotate data/vars1.sft tmp/temp.sft tmp/reg_annot.sft.rz 2> /dev/null
	exec cg select -rf {list} tmp/temp.sft tmp/temp2.sft
	exec diff tmp/temp2.sft data/expected-vars1-reg_annot.sft
} {} 

test reg_annot {2 compressed} {
	file copy -force data/vars1.sft data/reg_annot.sft data/reg_annot.sft.opt tmp
	exec cg razip tmp/reg_annot.sft tmp/vars1.sft 2> /dev/null
	exec cg annotate tmp/vars1.sft.rz tmp/temp.sft tmp/reg_annot.sft.rz 2> /dev/null
	exec cg select -rf {list} tmp/temp.sft tmp/temp2.sft
	exec diff tmp/temp2.sft data/expected-vars1-reg_annot.sft
} {} 

test reg_annot {basic, multiple fields} {
	cg select -f {chromosome begin end type ref alt} data/vars1.sft tmp/vars.sft
	file copy -force data/reg_annot.sft tmp/reg_annot.sft
	file_write tmp/reg_annot.sft.opt "fields\t{type begin end}\n"
	exec cg annotate tmp/vars.sft tmp/temp.sft tmp/reg_annot.sft 2> /dev/null
	exec diff tmp/temp.sft data/expected-vars1-reg_annot-multi.sft
} {} 

test reg_annot {near} {
	exec cg annotate -near 1000 data/vars1.sft tmp/temp.sft data/reg_annot.sft 2> /dev/null
	exec cg select -rf {list} tmp/temp.sft tmp/temp2.sft
	exec diff tmp/temp2.sft data/expected_near-vars1-reg_annot.sft
} {} 

test reg_annot {near indels} {
	exec cg select -q {$type == "del" || $type == "ins"} data/vars1.sft data/indels1.sft 2> /dev/null
	exec cg annotate -near 50 data/vars1.sft tmp/temp.sft data/indels1.sft 2> /dev/null
	exec cg select -rf {list} tmp/temp.sft tmp/temp2.sft
	exec diff tmp/temp2.sft data/expected-vars1-indels.sft
} {}

test reg_annot {sort error 1 in vars} {
	exec cg annotate data/vars_sorterror1.sft tmp/temp.sft data/reg_annot.sft
} {*Cannot annotate because the variant file (*) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test reg_annot {sort error 2 in vars} {
	exec cg annotate data/vars_sorterror2.sft tmp/temp.sft data/reg_annot.sft
} {*Cannot annotate because the variant file (*) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test reg_annot {sort error 3 in vars} {
	exec cg annotate data/vars_sorterror3.sft tmp/temp.sft data/reg_annot.sft
} {*Cannot annotate because the variant file (*) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test reg_annot {no sort error 4 in vars (not relevant for reg)} {
	catch {exec cg annotate data/vars_sorterror4.sft tmp/temp.sft data/reg_annot.sft 2> /dev/null}
} 0

test reg_annot {no sort error 5 in vars (not relevant for reg)} {
	catch {exec cg annotate data/vars_sorterror5.sft tmp/temp.sft data/reg_annot.sft 2> /dev/null}
} 0

test reg_annot {sort error 1 in database file} {
	file copy -force data/vars_sorterror1.sft tmp/reg_annot.sft
	exec cg annotate data/vars_annottest.sft tmp/temp.sft tmp/reg_annot.sft
} {*Cannot annotate because the database file is not correctly sorted (sort correctly using "cg select -s -")*} error match

test reg_annot {sort error 2 in database file} {
	file copy -force data/vars_sorterror2.sft tmp/reg_annot.sft
	exec cg annotate data/vars_annottest.sft tmp/temp.sft tmp/reg_annot.sft
} {*Cannot annotate because the database file is not correctly sorted (sort correctly using "cg select -s -")*} error match

test reg_annot {sort error 3 in database file} {
	file copy -force data/vars_sorterror3.sft tmp/reg_annot.sft
	exec cg annotate data/vars_annottest.sft tmp/temp.sft tmp/reg_annot.sft
} {*Cannot annotate because the database file is not correctly sorted (sort correctly using "cg select -s -")*} error match

test var_annot {basic} {
	exec cg annotate data/vars1.sft tmp/temp.sft data/var_annot.sft 2> /dev/null
	exec cg select -rf {list} tmp/temp.sft tmp/temp2.sft
	exec diff tmp/temp2.sft data/expected-vars1-var_annot.sft
} {} 

test var_annot {basic splitvar} {
	exec cg splitalleles data/vars1.sft > tmp/vars1.tsv
	exec cg annotate tmp/vars1.tsv tmp/temp.tsv data/var_annot.sft 2> /dev/null
	exec cg select -rf {list} tmp/temp.tsv tmp/temp2.tsv
	exec cg splitalleles data/expected-vars1-var_annot.sft > tmp/expected.tsv
	exec diff tmp/temp2.tsv tmp/expected.tsv
} {} 

test var_annot {basic multi} {
	file mkdir tmp
	cg select -f {chromosome begin end type ref alt} data/vars1.sft tmp/vars.sft
	file copy -force data/var_annot.sft tmp/var_annot.sft
	file_write tmp/var_annot.sft.opt "fields\t{name freq alt}\n"
	exec cg annotate tmp/vars.sft tmp/temp.sft tmp/var_annot.sft 2> /dev/null
	exec diff tmp/temp.sft data/expected-vars1-var_annot-multi.sft
} {} 

test var_annot {different types on same pos} {
	exec cg annotate data/vars2.tsv tmp/temp.tsv data/var_annot3.tsv 2> /dev/null
	exec diff tmp/temp.tsv data/expected-vars2-var_annot3.tsv
} {} 

test var_annot {sort error 1 in vars} {
	exec cg annotate data/vars_sorterror1.sft tmp/temp.sft data/var_annot.sft
} {*Cannot annotate because the variant file (*) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test var_annot {sort error 2 in vars} {
	exec cg annotate data/vars_sorterror2.sft tmp/temp.sft data/var_annot.sft
} {*Cannot annotate because the variant file (*) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test var_annot {sort error 3 in vars} {
	exec cg annotate data/vars_sorterror3.sft tmp/temp.sft data/var_annot.sft
} {*Cannot annotate because the variant file (*) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test var_annot {sort error 4 in vars} {
	exec cg annotate data/vars_sorterror4.sft tmp/temp.sft data/var_annot.sft
} {*Cannot annotate because the variant file (*) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test var_annot {sort error 5 in vars} {
	exec cg annotate data/vars_sorterror5.sft tmp/temp.sft data/var_annot.sft
} {*Cannot annotate because the variant file (*) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test var_annot {sort error 1 in database file} {
	file copy -force data/vars_sorterror1.sft tmp/var_annot.sft
	exec cg annotate data/vars_annottest.sft tmp/temp.sft tmp/var_annot.sft
} {*Cannot annotate because the database file is not correctly sorted (sort correctly using "cg select -s -")*} error match

test var_annot {sort error 2 in database file} {
	file copy -force data/vars_sorterror2.sft tmp/var_annot.sft
	exec cg annotate data/vars_annottest.sft tmp/temp.sft tmp/var_annot.sft
} {*Cannot annotate because the database file is not correctly sorted (sort correctly using "cg select -s -")*} error match

test var_annot {sort error 3 in database file} {
	file copy -force data/vars_sorterror3.sft tmp/var_annot.sft
	exec cg annotate data/vars_annottest.sft tmp/temp.sft tmp/var_annot.sft
} {*Cannot annotate because the database file is not correctly sorted (sort correctly using "cg select -s -")*} error match

test var_annot {sort error 4 in database file} {
	file copy -force data/vars_sorterror4.sft tmp/var_annot.sft
	exec cg annotate data/vars_annottest.sft tmp/temp.sft tmp/var_annot.sft
} {*Cannot annotate because the database file is not correctly sorted (sort correctly using "cg select -s -")*} error match

test var_annot {sort error 5 in database file} {
	file copy -force data/vars_sorterror5.sft tmp/var_annot.sft
	exec cg annotate data/vars_annottest.sft tmp/temp.sft tmp/var_annot.sft
} {*Cannot annotate because the database file is not correctly sorted (sort correctly using "cg select -s -")*} error match

test gene_annot {variant file sort error 1} {
	exec cg annotate -dbdir /complgen/refseq/hg18 data/vars_sorterror1.sft tmp/temp.sft data/gene_test.tsv
} {*Cannot annotate because the variant file is not correctly sorted (sort correctly using "cg select -s -")*} error match

test gene_annot {variant file sort error 2} {
	exec cg annotate -dbdir /complgen/refseq/hg18 data/vars_sorterror2.sft tmp/temp.sft data/gene_test.tsv
} {*Cannot annotate because the variant file is not correctly sorted (sort correctly using "cg select -s -")*} error match

test gene_annot {variant file sort error 3} {
	exec cg annotate -dbdir /complgen/refseq/hg18 data/vars_sorterror3.sft tmp/temp.sft data/gene_test.tsv
} {*Cannot annotate because the variant file is not correctly sorted (sort correctly using "cg select -s -")*} error match

test gene_annot {gene wrongly sorted database file error} {
	cg select -s - data/vars_annottest.sft tmp/vars_annottest.sft
	exec cg annotate -dbdir /complgen/refseq/hg18 tmp/vars_annottest.sft tmp/temp.sft data/gene_test-wrong1.tsv
} {*Cannot annotate because the database file (data/gene_test-wrong1.tsv) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test gene_annot {gene wrongly sorted database file error} {
	cg select -s - data/vars_annottest.sft tmp/vars_annottest.sft
	exec cg annotate -dbdir /complgen/refseq/hg18 tmp/vars_annottest.sft tmp/temp.sft data/gene_test-wrong2.tsv
} {*Cannot annotate because the database file (data/gene_test-wrong2.tsv) is not correctly sorted (sort correctly using "cg select -s -")*} error match

test gene_annot {gene} {
	cg select -s - data/vars_annottest.sft tmp/vars_annottest.sft
	exec cg annotate -dbdir /complgen/refseq/hg18 tmp/vars_annottest.sft tmp/temp.sft data/gene_test.tsv 2> /dev/null
	exec diff tmp/temp.sft data/expected-annotate-vars_annottest-gene_test.tsv
} {} 

test reg_annot {basic, extra comments} {
	file_write tmp/temp2.sft "# a comment\n"
	exec cat data/vars1.sft >> tmp/temp2.sft
	exec cg annotate tmp/temp2.sft tmp/temp.sft data/reg_annot.sft 2> /dev/null
	exec cg select -rf {list} tmp/temp.sft tmp/temp2.sft
	exec diff tmp/temp2.sft data/expected-vars1-reg_annot.sft
} {1d0
< # a comment	
child process exited abnormally} error

test var_annot {different types on same pos, extra comments} {
	file_write tmp/temp2.sft "# a comment\n"
	exec cat data/vars2.tsv >> tmp/temp2.sft
	exec cg annotate tmp/temp2.sft tmp/temp.tsv data/var_annot3.tsv 2> /dev/null
	exec diff tmp/temp.tsv data/expected-vars2-var_annot3.tsv
} {1d0
< # a comment	
child process exited abnormally} error

test gene_annot {gene, extra comments} {
	file_write tmp/temp2.sft "# a comment\n# another comment\n"
	exec cg select -s - data/vars_annottest.sft >> tmp/temp2.sft
	exec cg annotate -dbdir /complgen/refseq/hg18 tmp/temp2.sft tmp/temp.sft data/gene_test.tsv 2> /dev/null
	exec diff tmp/temp.sft data/expected-annotate-vars_annottest-gene_test.tsv
} {1,2d0
< # a comment	
< # another comment	 
child process exited abnormally} error

test reg_annot {existing field error} {
	exec cg annotate data/vars1.sft tmp/temp.sft data/reg_annot.sft 2> /dev/null
	exec cg annotate -near 1000 tmp/temp.sft tmp/temp2.sft data/reg_annot.sft
} {*Error: field(s) regtest already in file} match error

test reg_annot {replace} {
	exec cg annotate data/vars1.sft tmp/temp.sft data/reg_annot.sft 2> /dev/null
	exec cg annotate -near 1000 -replace 1 tmp/temp.sft tmp/temp2.sft data/reg_annot.sft 2> /dev/null
	exec cg select -rf {list} tmp/temp2.sft tmp/temp3.sft
	exec diff tmp/temp3.sft data/expected_near-vars1-reg_annot.sft
} {} 

test bcol_annot {basic} {
	test_cleantmp
	exec cg bcol make -p pos -c chromosome tmp/temp- coverage < data/cov.tsv
	file_write tmp/bcol_coverage.tsv "chromosome\tfile\nchr1\ttemp-chr1.bcol\nchr2\ttemp-chr2.bcol\n"
	exec cg annotate data/bcol_annot-test.tsv tmp/annot_test.tsv tmp/bcol_coverage.tsv 2> /dev/null
	exec diff tmp/annot_test.tsv data/expected-bcol_annot-test.tsv
} {} 

test bcol_annot {header error} {
	test_cleantmp
	exec cg bcol make -p pos -c chromosome tmp/temp- coverage < data/cov.tsv
	file_write tmp/bcol_coverage.tsv "chr1\ttemp-chr1.bcol\nchr2\ttemp-chr2.bcol\n"
	exec cg annotate data/bcol_annot-test.tsv tmp/annot_test.tsv tmp/bcol_coverage.tsv
} {*bcol database (tmp/bcol_coverage.tsv) should have a header of the type: chromosome file*} error match

#test bcol_annot {basic compressed} {
#	test_cleantmp
#	exec cg bcol make -p pos -c chromosome tmp/temp- coverage < data/cov.tsv
#	file_write tmp/bcol_coverage.tsv "chromosome\tfile\nchr1\ttemp-chr1.bcol\nchr2\ttemp-chr2.bcol\n"
#	file copy data/bcol_annot-test.tsv tmp/bcol_annot-test.tsv
#	cg razip tmp/bcol_annot-test.tsv 2> /dev/null
#	exec cg annotate tmp/bcol_annot-test.tsv.gz tmp/annot_test.tsv tmp/bcol_coverage.tsv 2> /dev/null
#	exec diff tmp/annot_test.tsv data/expected-bcol_annot-test.tsv
#} {} 

test bcol_annot {basic uncompressed bcol} {
	test_cleantmp
	exec cg bcol make -p pos -c chromosome tmp/temp- coverage < data/cov.tsv
	cg unzip {*}[glob tmp/*.rz]
	file_write tmp/bcol_coverage.tsv "chromosome\tfile\nchr1\ttemp-chr1.bcol\nchr2\ttemp-chr2.bcol\n"
	exec cg annotate data/bcol_annot-test.tsv tmp/annot_test.tsv tmp/bcol_coverage.tsv 2> /dev/null
	exec diff tmp/annot_test.tsv data/expected-bcol_annot-test.tsv
} {} 

file delete -force tmp/temp.sft
file delete -force tmp/temp2.sft

set ::env(PATH) $keeppath

testsummarize
