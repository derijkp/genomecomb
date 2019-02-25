#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl
set keepdir [pwd]

test download {download_ucsc} {
	cg download_ucsc tmp/ucsc_hg19_tRNAs.tsv hg19 tRNAs
	list [lindex [exec md5sum tmp/ucsc_hg19_tRNAs.tsv] 0] [lindex [exec md5sum tmp/ucsc_hg19_tRNAs.tsv.info] 0]
} {bfa24aa5f51cd136bfd693ea19a679dc 249580e3df6f7cfa9cec70751aec74bb}

test download {download_genes} {
	cg download_genes tmp/gene_hg19_refGene.tsv hg19 refGene
	exec md5sum {*}[lsort -dict [glob tmp/*tsv tmp/*tsv.gz tmp/*tsv.info]]
} {2085883d03f3c2e6bb2319bd76ed89cd  tmp/gene_hg19_refGene.tsv
56b08c46b4aef22ca9f332954e9ad598  tmp/gene_hg19_refGene.tsv.gz
5c103f26a531d7e5b4ca34b51319c96d  tmp/gene_hg19_refGene.tsv.info}

test intgene {intgene} {
	exec cg intgene --stack 1 data/gene_test.tsv > tmp/gene_test.tsv
	write_tab tmp/genes_pre.tsv [deindent {
		chromosome	start	end	name	strand	bin	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	id	name2	cdsStartStat	cdsEndStat	exonFrames
		chr1	4224	19233	NR_024540a	-	585	19233	19233	11	4224,4832,5658,6469,6720,7095,7468,7777,8130,14600,19183,	4692,4901,5810,6628,6918,7231,7605,7924,8229,14754,19233,	0	WASH7P	unk	unk	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		chr1	4224	19233	NR_024540b	-	585	19233	19233	11	4225,4832,5658,6469,6720,7095,7468,7777,8130,14600,19183,	4692,4901,5810,6628,6918,7231,7605,7924,8229,14754,19233,	0	WASH7P	unk	unk	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		chr1	4224	19233	NR_024540c	-	585	19236	19233	11	4224,4832,5658,6469,6720,7095,7468,7777,8130,14600,19183,	4692,4901,5810,6628,6918,7231,7605,7924,8229,14754,19233,	0	WASH7P	unk	unk	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		chr1	4224	19233	NR_024540c2	-	585	19236	19233	11	4224,4832,5658,6469,6720,7095,7468,7777,8130,14600,19183,	4692,4901,5810,6628,6918,7231,7605,7924,8229,14754,19233,	0	WASH7P	unk	unk	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		chr1	4224	19234	NR_024540d	-	585	19233	19233	11	4224,4832,5658,6469,6720,7095,7468,7777,8130,14600,19183,	4692,4901,5810,6628,6918,7231,7605,7924,8229,14754,19233,	0	WASH7P	unk	unk	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		chr22	20229957	20235373	NM_001128635b	-	739	20230345	20235265	1	20229957,	20235373,	0	RIMBP3B	cmpl	cmpl	0,
	}]
	exec cg intgene tmp/genes_pre.tsv tmp/gene_test.tsv > tmp/gene_result.tsv
	set result [checkdiff tmp/gene_test.tsv tmp/gene_result.tsv]
	set expected [string trim [deindent {
		Files differ: tmp/gene_test.tsv tmp/gene_result.tsv
		3c3,6
		< chr1	4224	19233	NR_024540	-	gene_test	585	19233	19233	11	4224,4832,5658,6469,6720,7095,7468,7777,8130,14600,19183,	4692,4901,5810,6628,6918,7231,7605,7924,8229,14754,19233,	0	WASH7P	unk	unk	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		---
		> chr1	4224	19233	NR_024540a	-	genes_pre	585	19233	19233	11	4224,4832,5658,6469,6720,7095,7468,7777,8130,14600,19183,	4692,4901,5810,6628,6918,7231,7605,7924,8229,14754,19233,	0	WASH7P	unk	unk	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		> chr1	4224	19233	NR_024540b	-	genes_pre	585	19233	19233	11	4225,4832,5658,6469,6720,7095,7468,7777,8130,14600,19183,	4692,4901,5810,6628,6918,7231,7605,7924,8229,14754,19233,	0	WASH7P	unk	unk	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		> chr1	4224	19233	NR_024540c	-	genes_pre	585	19236	19233	11	4224,4832,5658,6469,6720,7095,7468,7777,8130,14600,19183,	4692,4901,5810,6628,6918,7231,7605,7924,8229,14754,19233,	0	WASH7P	unk	unk	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		> chr1	4224	19234	NR_024540d	-	genes_pre	585	19233	19233	11	4224,4832,5658,6469,6720,7095,7468,7777,8130,14600,19183,	4692,4901,5810,6628,6918,7231,7605,7924,8229,14754,19233,	0	WASH7P	unk	unk	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		991c994
		< chr22	20229957	20235373	NM_001128635	-	gene_test	739	20230345	20235265	1	20229957,	20235373,	0	RIMBP3B	cmpl	cmpl	0,
		---
		> chr22	20229957	20235373	NM_001128635b	-	genes_pre	739	20230345	20235265	1	20229957,	20235373,	0	RIMBP3B	cmpl	cmpl	0,
		child process exited abnormally
	}]]
	expr {$result eq $expected}
} 1

testsummarize
