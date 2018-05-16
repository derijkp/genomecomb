#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

set keepdir [pwd]

test homwes {basic} {
	file delete -force $::bigtestdir/tmp/homwes
	file mkdir $::bigtestdir/tmp/homwes
	set varfile $::bigtestdir/ori/exomes_yri.ori/compar/annot_compar-yri_exome_test.tsv.lz4
	set dbdir $::bigtestdir/refseqtest/hg19
	cg homwes --stack 1 -dbdir $dbdir $varfile NA19240 $::bigtestdir/tmp/homwes/homwes-out.tsv
	cg tsvdiff -q 1 -x *.log $::bigtestdir/tmp/homwes $::bigtestdir/expected/yri_exome.homwes-out
} {}

test homwes {samples empty -> multiple samples} {
	file delete -force $::bigtestdir/tmp/homwes_multi
	file mkdir $::bigtestdir/tmp/homwes_multi
	set varfile $::bigtestdir/ori/exomes_yri.ori/compar/annot_compar-yri_exome_test.tsv.lz4
	set dbdir $::bigtestdir/refseqtest/hg19
	cg homwes --stack 1 -dbdir $dbdir $varfile {} $::bigtestdir/tmp/homwes_multi/homwes-out.tsv
	cg tsvdiff -q 1 -x *.log -x log_jobs $::bigtestdir/tmp/homwes_multi $::bigtestdir/expected/homwes_multi
} {}

test homwes {sample not given -> multiple samples} {
	file delete -force $::bigtestdir/tmp/homwes_multi
	file mkdir $::bigtestdir/tmp/homwes_multi
	mklink $::bigtestdir/ori/exomes_yri.ori/compar/annot_compar-yri_exome_test.tsv.lz4 $::bigtestdir/tmp/homwes_multi/annot_compar-yri_exome_test.tsv.lz4
	set varfile $::bigtestdir/tmp/homwes_multi/annot_compar-yri_exome_test.tsv.lz4
	set dbdir $::bigtestdir/refseqtest/hg19
	cg homwes --stack 1 -dbdir $dbdir $varfile {} $::bigtestdir/tmp/homwes_multi/homwes-out.tsv
	file delete $::bigtestdir/tmp/homwes_multi/annot_compar-yri_exome_test.tsv.lz4
	cg tsvdiff -q 1 -x *.log -x log_jobs -x .pversion $::bigtestdir/tmp/homwes_multi $::bigtestdir/expected/homwes_multi
} {}

test homwes {vcf source, sample not given -> multiple samples} {
	file delete -force $::bigtestdir/tmp/homwes_vcfmulti
	file mkdir $::bigtestdir/tmp/homwes_vcfmulti
	mklink $::bigtestdir/ori/exomes_yri.ori/compar/annot_compar-yri_exome_test.vcf.gz $::bigtestdir/tmp/homwes_vcfmulti/annot_compar-yri_exome_test.vcf.gz
	set varfile $::bigtestdir/tmp/homwes_vcfmulti/annot_compar-yri_exome_test.vcf.gz
	set dbdir $::bigtestdir/refseqtest/hg19
	cg homwes --stack 1 -callers {} -dbdir $dbdir $varfile {} $::bigtestdir/tmp/homwes_vcfmulti/homwes-out.tsv
	file delete $::bigtestdir/tmp/homwes_vcfmulti/annot_compar-yri_exome_test.vcf.gz
	file delete -force $::bigtestdir/tmp/homwes_vcfmulti/homwes-out.work/annot_compar-yri_exome_test.vcf.tsv.temp.index
	cg tsvdiff -q 1 -x *.log -x .pversion -x log_jobs $::bigtestdir/tmp/homwes_vcfmulti $::bigtestdir/expected/homwes_vcfmulti
} {}

test homwes {field not found error} {
	write_tab tmp/temp.tsv {
		chromosome	begin	end
		1	10	20
	}
	set dbdir $::bigtestdir/refseqtest/hg19
	cg homwes -callers {} -dbdir $dbdir tmp/temp.tsv
} {Could not find alleleSeq1 field for sample sample1 in header (checked alleleSeq1-sample1, alleleSeq1-sample1, alleleSeq1)} error

cd $keepdir

testsummarize

