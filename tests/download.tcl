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

testsummarize
