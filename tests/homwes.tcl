#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

set keepdir [pwd]

#cg select -overwrite 1 \
#    -f {chromosome begin end type ref alt quality-* sequenced-* zyg-* filter-* alleleSeq* coverage-* genoqual-* snp* simpleRepeat rmsk microsat cytoBand homopolymer_* clinvar* refGene* ensGene* mir* } \
#    -q {region("chr2:100000000-180000000","chr3:30000000-50000000","chr3:164000000-170000000","chr5:100000000-150000000","chr6:40000000-44000000","chr11:80000000-90000000")} \
#     $::smalltestdir/ori/annot_compar-yri_exome_test.tsv.zst \
# data/annot_compar-yri_exome_part_test.tsv

test homwes {basic} {
	file delete -force tmp/homwes
	file mkdir tmp/homwes
	set varfile tmp/annot_compar-yri_exome_part_test.tsv
	mklink [gzfile data/annot_compar-yri_exome_part_test.tsv] $varfile
	set dbdir $::refseqdir/hg19
	cg homwes --stack 1 -dbdir $dbdir $varfile NA19240 tmp/homwes/homwes-out.tsv
	cg tsvdiff -q 1 -x *.log tmp/homwes data/homwes
} {}

test homwes {basic -variantsonly 0} {
	file delete -force tmp/homwes_variantsonly0
	file mkdir tmp/homwes_variantsonly0
	set varfile tmp/annot_compar-yri_exome_part_test.tsv
	mklink [gzfile data/annot_compar-yri_exome_part_test.tsv] $varfile
	set dbdir $::refseqdir/hg19
	cg homwes --stack 1 -variantsonly 0 -dbdir $dbdir $varfile NA19240 tmp/homwes_variantsonly0/homwes-out.tsv
	cg tsvdiff tmp/homwes_variantsonly0/homwes-out.tsv data/homwes_variantsonly0-homwes-out.tsv
} {}

test homwes {samples empty -> multiple samples} {
	file delete -force tmp/homwes_multi
	file mkdir tmp/homwes_multi
	set varfile tmp/annot_compar-yri_exome_part_test.tsv
	mklink [gzfile data/annot_compar-yri_exome_part_test.tsv] $varfile
	set dbdir $::refseqdir/hg19
	cg homwes --stack 1 -dbdir $dbdir $varfile {} tmp/homwes_multi/homwes-out.tsv
	file delete -force tmp/homwes_multi/homwes-out.work
	cg tsvdiff -q 1 -x *.log -x log_jobs tmp/homwes_multi data/homwes_multi
} {}

test homwes {sample not given -> multiple samples} {
	file delete -force tmp/homwes_multi
	file mkdir tmp/homwes_multi
	set srcfile tmp/annot_compar-yri_exome_part_test.tsv
	mklink [gzfile data/annot_compar-yri_exome_part_test.tsv] $srcfile
	set varfile tmp/homwes_multi/annot_compar-yri_exome_part_test.tsv[gzext $srcfile]
	mklink $srcfile $varfile
	set dbdir $::refseqdir/hg19
	cg homwes --stack 1 -dbdir $dbdir $varfile {} tmp/homwes_multi/homwes-out.tsv
	file delete $varfile
	file delete -force tmp/homwes_multi/homwes-out.work
	cg tsvdiff -q 1 -x *.log -x log_jobs -x .pversion tmp/homwes_multi data/homwes_multi
} {}

test homwes {vcf source, sample not given -> multiple samples} {
	file delete -force tmp/homwes_vcfmulti
	file mkdir tmp/homwes_vcfmulti
	mklink data/annot_compar-yri_exome_test.vcf tmp/homwes_vcfmulti/annot_compar-yri_exome_part_test.vcf
	set varfile tmp/homwes_vcfmulti/annot_compar-yri_exome_part_test.vcf
	set dbdir $::refseqdir/hg19
	cg homwes --stack 1 -callers {} -dbdir $dbdir $varfile {} tmp/homwes_vcfmulti/homwes-out.tsv
	file delete -force tmp/homwes_vcfmulti/annot_compar-yri_exome_part_test.vcf tmp/homwes_vcfmulti/homwes-out.work tmp/homwes_vcfmulti/log_jobs
	cg tsvdiff -q 1 \
		-x *.log -x .pversion -x log_jobs \
		-x *.analysisinfo \
		-ignorefields {annotate_cg_version} \
		tmp/homwes_vcfmulti data/homwes_vcfmulti
} {}

test homwes {field not found error} {
	write_tab tmp/temp.tsv {
		chromosome	begin	end
		1	10	20
	}
	set dbdir $::refseqdir/hg19
	cg homwes -callers {} -dbdir $dbdir tmp/temp.tsv
} {Using callers: 
Could not find alleleSeq1 field for sample sample1 in header (checked alleleSeq1-sample1, alleleSeq1-sample1, alleleSeq1)} error

cd $keepdir

testsummarize

