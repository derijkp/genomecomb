#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test reports {hsmetrics} {
	test_cleantmp
	cg select -q {$chromosome in "chr21 chr22"} $::refseqdir/hg19/extra/reg_hg19_exome_SeqCap_EZ_v3.tsv.lz4 tmp/regfile.tsv
	set bamfile genomecomb.testdata/ori/test-map-rdsbwa-NA19240chr2122.bam
	set regionfile tmp/regfile.tsv
	set resultfile tmp/result.hsmetrics
	cg hsmetrics $bamfile $regionfile $resultfile 2> /dev/null
	cg select -rc 1 tmp/result.hsmetrics tmp/result.hsmetrics.nocomments
	exec diff tmp/result.hsmetrics.nocomments genomecomb.testdata/expected/bam_histo-NA19240chr2122.hsmetrics
} {}

test reports {process_reports} {
	cd $::bigtestdir
	file delete -force tmp/test_reports
	file mkdir tmp/test_reports/NA19240mx2
	file copy {*}[glob expected/exomes_yri_mx2/samples/NA19240mx2/*] tmp/test_reports/NA19240mx2
	file delete -force tmp/test_reports/NA19240mx2/reports
	mklink refseqtest/hg19/extra/reg_hg19_exome_SeqCap_EZ_v3.tsv.lz4 tmp/test_reports/NA19240mx2/reg_hg19_targets.tsv.lz4
	cg process_reports -stack 1 -v 2 tmp/test_reports/NA19240mx2 refseqtest/hg19 2>@ stderr >@ stdout
	cg tsvdiff -q 1 tmp/test_reports/NA19240mx2/reports expected/exomes_yri_mx2/samples/NA19240mx2/reports
} {}

test reports {process_reports no targetfile} {
	cd $::bigtestdir
	file delete -force tmp/test_reportsnotarget
	file mkdir tmp/test_reportsnotarget/NA19240mx2
	file copy {*}[glob expected/exomes_yri_mx2/samples/NA19240mx2/*] tmp/test_reportsnotarget/NA19240mx2
	file delete -force tmp/test_reportsnotarget/NA19240mx2/reports tmp/test_reportsnotarget/NA19240mx2/reg_hg19_targets.tsv.lz4
	cg process_reports -stack 1 -v 2 tmp/test_reportsnotarget/NA19240mx2 refseqtest/hg19 2>@ stderr >@ stdout
	cg tsvdiff -q 1 tmp/test_reportsnotarget/NA19240mx2/reports expected/test_reportsnotarget/NA19240mx2/reports
} {}

testsummarize
