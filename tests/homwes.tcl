#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

set keepdir [pwd]

test homwes {basic} {
	file delete -force $::smalltestdir/tmp/homwes
	file mkdir $::smalltestdir/tmp/homwes
	set varfile [gzfile $::smalltestdir/ori/annot_compar-yri_exome_test.tsv]
	set dbdir $::smalltestdir/refseqtest/hg19
	cg homwes --stack 1 -dbdir $dbdir $varfile NA19240 $::smalltestdir/tmp/homwes/homwes-out.tsv
	cg tsvdiff -q 1 -x *.log $::smalltestdir/tmp/homwes $::smalltestdir/expected/yri_exome.homwes-out
} {}

test homwes {samples empty -> multiple samples} {
	file delete -force $::smalltestdir/tmp/homwes_multi
	file mkdir $::smalltestdir/tmp/homwes_multi
	set varfile [gzfile $::smalltestdir/ori/annot_compar-yri_exome_test.tsv]
	set dbdir $::smalltestdir/refseqtest/hg19
	cg homwes --stack 1 -dbdir $dbdir $varfile {} $::smalltestdir/tmp/homwes_multi/homwes-out.tsv
	cg tsvdiff -q 1 -x *.log -x log_jobs $::smalltestdir/tmp/homwes_multi $::smalltestdir/expected/homwes_multi
} {}

test homwes {sample not given -> multiple samples} {
	file delete -force $::smalltestdir/tmp/homwes_multi
	file mkdir $::smalltestdir/tmp/homwes_multi
	set srcfile [gzfile $::smalltestdir/ori/annot_compar-yri_exome_test.tsv]
	set varfile $::smalltestdir/tmp/homwes_multi/annot_compar-yri_exome_test.tsv[file extension $srcfile]
	mklink $srcfile $varfile
	set dbdir $::smalltestdir/refseqtest/hg19
	cg homwes --stack 1 -dbdir $dbdir $varfile {} $::smalltestdir/tmp/homwes_multi/homwes-out.tsv
	file delete $varfile
	cg tsvdiff -q 1 -x *.log -x log_jobs -x .pversion $::smalltestdir/tmp/homwes_multi $::smalltestdir/expected/homwes_multi
} {}

test homwes {vcf source, sample not given -> multiple samples} {
	file delete -force $::smalltestdir/tmp/homwes_vcfmulti
	file mkdir $::smalltestdir/tmp/homwes_vcfmulti
	mklink $::smalltestdir/ori/annot_compar-yri_exome_test.vcf.gz $::smalltestdir/tmp/homwes_vcfmulti/annot_compar-yri_exome_test.vcf.gz
	set varfile $::smalltestdir/tmp/homwes_vcfmulti/annot_compar-yri_exome_test.vcf.gz
	set dbdir $::smalltestdir/refseqtest/hg19
	cg homwes --stack 1 -callers {} -dbdir $dbdir $varfile {} $::smalltestdir/tmp/homwes_vcfmulti/homwes-out.tsv
	file delete $::smalltestdir/tmp/homwes_vcfmulti/annot_compar-yri_exome_test.vcf.gz
	file delete -force $::smalltestdir/tmp/homwes_vcfmulti/homwes-out.work/annot_compar-yri_exome_test.vcf.tsv.temp.index
	cg tsvdiff -q 1 -x *.log -x .pversion -x log_jobs $::smalltestdir/tmp/homwes_vcfmulti $::smalltestdir/expected/homwes_vcfmulti
} {}

test homwes {field not found error} {
	write_tab tmp/temp.tsv {
		chromosome	begin	end
		1	10	20
	}
	set dbdir $::smalltestdir/refseqtest/hg19
	cg homwes -callers {} -dbdir $dbdir tmp/temp.tsv
} {Using callers: 
Could not find alleleSeq1 field for sample sample1 in header (checked alleleSeq1-sample1, alleleSeq1-sample1, alleleSeq1)} error

cd $keepdir

testsummarize

