#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl
set keepdir [pwd]
if {[llength [get argv ""]]} {
	set dopts [get argv ""]
} else {
	set dopts {--stack 1 --verbose 2}
}
set test_cleantmp 0

test sv {manta} {
	cd $::bigtestdir
	file delete -force tmp/sv-manta
	file mkdir tmp/sv-manta
	mklink /data/genomecomb.testdata/ori/sv/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam tmp/sv-manta/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam
	mklink /data/genomecomb.testdata/ori/sv/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam.bai tmp/sv-manta/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam.bai
	cg sv_manta {*}$::dopts -refseq $::bigtestdir/refseqtest/hg19 tmp/sv-manta/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam
	cg tsvdiff -x *.xml -x svLocusGraphStats.tsv -x *.tbi tmp/sv-manta expected/sv-manta
} {}

test sv {cg sv -method manta, giving resultfile} {
	cd $::bigtestdir
	file delete -force tmp/sv-manta
	file mkdir tmp/sv-manta
	mklink /data/genomecomb.testdata/ori/sv/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam tmp/sv-manta/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam
	mklink /data/genomecomb.testdata/ori/sv/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam.bai tmp/sv-manta/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam.bai
	cg sv {*}$::dopts -method manta -refseq $::bigtestdir/refseqtest/hg19 tmp/sv-manta/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam tmp/sv-manta/resultsv.tsv
	cg zcat expected/sv-manta/sv-manta-dsbwa-ERR194147_30x_NA12878-chr21part.tsv.lz4 > tmp/sv-manta/expected.tsv
	exec diff tmp/sv-manta/resultsv.tsv tmp/sv-manta/expected.tsv | grep -v \#\#cmdline | grep -v \#\#fileDate
} {5c5
---
128c128
---
child process exited abnormally} error

testsummarize

