#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test reports {hsmetrics} {
	test_cleantmp
	cg select -q {$chromosome in "chr21 chr22"} $::refseqdir/hg19/extra/reg_hg19_exome_SeqCap_EZ_v3.tsv tmp/regfile.tsv
	set bamfile genomecomb.testdata/ori/test-map-rdsbwa-NA19240chr2122.bam
	set regionfile tmp/regfile.tsv
	set resultfile tmp/result.hsmetrics
	cg hsmetrics $bamfile $regionfile $resultfile 2> /dev/null
	cg select -rc 1 tmp/result.hsmetrics tmp/result.hsmetrics.nocomments
	exec diff tmp/result.hsmetrics.nocomments genomecomb.testdata/expected/bam_histo-NA19240chr2122.hsmetrics
} {}

testsummarize
