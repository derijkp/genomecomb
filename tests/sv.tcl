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

test sv {sniffles} {
	cd $::smalltestdir
	file delete -force tmp/sv-sniffles
	file mkdir tmp/sv-sniffles
	mklink ori/sv/ont/map-sngmlr-NA12878.bam tmp/sv-sniffles/map-sngmlr-NA12878.bam
	mklink ori/sv/ont/map-sngmlr-NA12878.bam.bai tmp/sv-sniffles/map-sngmlr-NA12878.bam.bai
	cg sv_sniffles {*}$::dopts -refseq $::smalltestdir/refseqtest/hg19 tmp/sv-sniffles/map-sngmlr-NA12878.bam
	cg tsvdiff -x *.xml -x *.vcf -x svLocusGraphStats.tsv -x *.tbi \
		tmp/sv-sniffles expected/sv-sniffles
} {}

test sv {manta} {
	cd $::smalltestdir
	file delete -force tmp/sv-manta
	file mkdir tmp/sv-manta
	mklink ori/sv/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam tmp/sv-manta/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam
	mklink ori/sv/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam.bai tmp/sv-manta/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam.bai
	cg sv_manta {*}$::dopts -refseq $::smalltestdir/refseqtest/hg19 tmp/sv-manta/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam
	cg tsvdiff -x *.xml -x svLocusGraphStats.tsv -x *.tbi -x *.py -x *.py.* -x alignmentStatsSummary.txt \
		tmp/sv-manta expected/sv-manta
} {}

test sv {cg sv -method manta, giving resultfile} {
	cd $::smalltestdir
	file delete -force tmp/sv-manta
	file mkdir tmp/sv-manta
	mklink ori/sv/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam tmp/sv-manta/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam
	mklink ori/sv/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam.bai tmp/sv-manta/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam.bai
	cg sv {*}$::dopts -method manta -refseq $::smalltestdir/refseqtest/hg19 tmp/sv-manta/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam tmp/sv-manta/resultsv.tsv
	cg zcat expected/sv-manta/sv-manta-dsbwa-ERR194147_30x_NA12878-chr21part.tsv.lz4 > tmp/sv-manta/expected.tsv
	cg tsvdiff tmp/sv-manta/resultsv.tsv tmp/sv-manta/expected.tsv
} {}

test sv {lumpy} {
	cd $::smalltestdir
	file delete -force tmp/sv-lumpy
	file mkdir tmp/sv-lumpy
	mklink ori/sv/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam tmp/sv-lumpy/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam
	mklink ori/sv/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam.bai tmp/sv-lumpy/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam.bai
	exec cg sv_lumpy {*}$::dopts -refseq $::smalltestdir/refseqtest/hg19 tmp/sv-lumpy/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam
	cg tsvdiff -x *.tbi tmp/sv-lumpy expected/sv-lumpy
} {}

test sv {gridss} {
	cd $::smalltestdir
	file delete -force tmp/sv-gridss
	file mkdir tmp/sv-gridss
	mklink ori/sv/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam tmp/sv-gridss/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam
	mklink ori/sv/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam.bai tmp/sv-gridss/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam.bai
	exec cg sv_gridss {*}$::dopts -refseq $::smalltestdir/refseqtest/hg19 tmp/sv-gridss/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam
	cg tsvdiff -x *.xml -x svLocusGraphStats.tsv -x *.tbi tmp/sv-gridss expected/sv-gridss
} {}

testsummarize

