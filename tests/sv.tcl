#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl
set keepdir [pwd]

set test_cleantmp 0


test sv {sniffles} {
	cd $::smalltestdir
	file delete -force tmp/sv-sniffles
	file mkdir tmp/sv-sniffles
	mklink ori/sv/ont/map-sngmlr-NA12878.bam tmp/sv-sniffles/map-sngmlr-NA12878.bam
	mklink ori/sv/ont/map-sngmlr-NA12878.bam.bai tmp/sv-sniffles/map-sngmlr-NA12878.bam.bai
	cg sv_sniffles {*}$::dopts \
		-refseq $::smalltestdir/refseqtest/hg19 \
		tmp/sv-sniffles/map-sngmlr-NA12878.bam \
		>& tmp/sv-sniffles.log
	cg tsvdiff -q 1 -x *.xml -x *.vcf -x svLocusGraphStats.tsv -x *.tbi \
		-x *.snf \
		-ignorefields {varcaller_cg_version} \
		tmp/sv-sniffles expected/sv-sniffles
} {}

# samtools view -b /data/projects/promethion-pub.data/pubminion/NA12878-nanopore-wgs/bwa-mem_NA12878_25FC.bam chr19:19000000-22000000 chr21:21000000-26000000 > ori/ont/bwa-mem_NA12878_25FC_part19_21.bam
test sv {sniffles2 NA12878} {
	cd $::smalltestdir
	file delete -force tmp/sv-sniffles2
	file mkdir tmp/sv-sniffles2
	mklink ori/ont/bwa-mem_NA12878_25FC_part19_21.bam tmp/sv-sniffles2/bwa-mem_NA12878_25FC_part19_21.bam
	mklink ori/ont/bwa-mem_NA12878_25FC_part19_21.bam.bai tmp/sv-sniffles2/bwa-mem_NA12878_25FC_part19_21.bam.bai
	cg sv_sniffles {*}$::dopts \
		-refseq $::smalltestdir/refseqtest/hg19 \
		tmp/sv-sniffles2/bwa-mem_NA12878_25FC_part19_21.bam \
		>& tmp/sv-sniffles2.log
	cg tsvdiff -q 1 -x *.tbi -x *.vcf \
		-x *.snf \
		-ignorefields {varcaller_cg_version} \
		tmp/sv-sniffles2 expected/sv-sniffles2
} {}

# samtools view -b /data/projects/promethion-pub.data/pubminion/NA12878-nanopore-wgs/bwa-mem_NA12878_25FC.bam chr19:19000000-22000000 chr21:21000000-26000000 > ori/ont/bwa-mem_NA12878_25FC_part19_21.bam
test sv {sniffles2 ngmlr NA12878} {
	cd $::smalltestdir
	file delete -force tmp/sv-sniffles2_ngmlr
	file mkdir tmp/sv-sniffles2_ngmlr
	mklink ori/ont/map-ngmlr-NA12878_25FC_part19_21.bam tmp/sv-sniffles2_ngmlr/map-ngmlr-NA12878_25FC_part19_21.bam
	mklink ori/ont/map-ngmlr-NA12878_25FC_part19_21.bam.bai tmp/sv-sniffles2_ngmlr/map-ngmlr-NA12878_25FC_part19_21.bam.bai
	cg sv_sniffles {*}$::dopts \
		-refseq $::smalltestdir/refseqtest/hg19 \
		tmp/sv-sniffles2_ngmlr/map-ngmlr-NA12878_25FC_part19_21.bam \
		>& tmp/sv-sniffles2_ngmlr.log
	cg tsvdiff -q 1 -x *.tbi -x *.vcf \
		-x *.snf \
		-ignorefields {varcaller_cg_version} \
		tmp/sv-sniffles2_ngmlr expected/sv-sniffles2_ngmlr
} {}

# samtools view -b /data/projects/promethion-pub.data/pubminion/NA12878-nanopore-wgs/bwa-mem_NA12878_25FC.bam chr19:19000000-22000000 chr21:21000000-26000000 > ori/ont/bwa-mem_NA12878_25FC_part19_21.bam
test sv {sniffles_minimap2 NA12878 minimap2} {
	cd $::smalltestdir
	file delete -force tmp/sv-sniffles2_minimap2
	file mkdir tmp/sv-sniffles2_minimap2
	mklink ori/ont/map-minimap2-NA12878_25FC_part19_21.bam tmp/sv-sniffles2_minimap2/map-minimap2-NA12878_25FC_part19_21.bam
	mklink ori/ont/map-minimap2-NA12878_25FC_part19_21.bam.bai tmp/sv-sniffles2_minimap2/map-minimap2-NA12878_25FC_part19_21.bam.bai
	cg sv_sniffles {*}$::dopts \
		-refseq $::smalltestdir/refseqtest/hg19 \
		tmp/sv-sniffles2_minimap2/map-minimap2-NA12878_25FC_part19_21.bam \
		>& tmp/sv-sniffles2_minimap2.log
	cg tsvdiff -q 1 -x *.tbi -x *.vcf \
		-x *.snf \
		-ignorefields {varcaller_cg_version} \
		tmp/sv-sniffles2_minimap2 expected/sv-sniffles2_minimap2
} {}

test sv {sniffles -distrreg chr} {
	cd $::smalltestdir
	file delete -force tmp/sv-sniffles2_minimap2_distrchr
	file mkdir tmp/sv-sniffles2_minimap2_distrchr
	mklink ori/ont/map-minimap2-NA12878_25FC_part19_21.bam tmp/sv-sniffles2_minimap2_distrchr/map-minimap2-NA12878_25FC_part19_21.bam
	mklink ori/ont/map-minimap2-NA12878_25FC_part19_21.bam.bai tmp/sv-sniffles2_minimap2_distrchr/map-minimap2-NA12878_25FC_part19_21.bam.bai
	cg sv -stack 1 -v 2 {*}$::dopts \
		-method sniffles -distrreg chr \
		-refseq $::smalltestdir/refseqtest/hg19 \
		tmp/sv-sniffles2_minimap2_distrchr/map-minimap2-NA12878_25FC_part19_21.bam \
		>& tmp/sv-sniffles2_minimap2_distrchr.log
	cg tsvdiff -q 1 -x *.tbi -x *.vcf -x *.zsti \
		-x *.snf \
		-ignorefields {varcaller_cg_version} \
		tmp/sv-sniffles2_minimap2_distrchr expected/sv-sniffles2_minimap2
} {}

test sv {cuteSV} {
	cd $::smalltestdir
	file delete -force tmp/sv-cuteSV
	file mkdir tmp/sv-cuteSV
	mklink ori/sv/ont/map-sngmlr-NA12878.bam tmp/sv-cuteSV/map-sngmlr-NA12878.bam
	mklink ori/sv/ont/map-sngmlr-NA12878.bam.bai tmp/sv-cuteSV/map-sngmlr-NA12878.bam.bai
	cg sv_cuteSV {*}$::dopts \
		-refseq $::smalltestdir/refseqtest/hg19 \
		tmp/sv-cuteSV/map-sngmlr-NA12878.bam \
		>& tmp/sv-cuteSV.log
	cg tsvdiff -q 1 \
		-x *.xml -x *.vcf -x svLocusGraphStats.tsv \
		-x *.zsti -x *.tbi \
		-ignorefields {varcaller_cg_version} \
		tmp/sv-cuteSV expected/sv-cuteSV
} {}

test sv {cuteSV -distrreg chr} {
	cd $::smalltestdir
	file delete -force tmp/sv-cuteSV
	file mkdir tmp/sv-cuteSV
	mklink ori/sv/ont/map-sngmlr-NA12878.bam tmp/sv-cuteSV/map-sngmlr-NA12878.bam
	mklink ori/sv/ont/map-sngmlr-NA12878.bam.bai tmp/sv-cuteSV/map-sngmlr-NA12878.bam.bai
	cg sv -method cuteSV {*}$::dopts \
		-refseq $::smalltestdir/refseqtest/hg19 \
		-distrreg chr \
		tmp/sv-cuteSV/map-sngmlr-NA12878.bam
	cg tsvdiff -q 1 -x *.xml -x *.vcf -x svLocusGraphStats.tsv -x *.tbi -x *.submitting -x *.zsti \
		-ignorefields {varcaller_cg_version} \
		tmp/sv-cuteSV expected/sv-cuteSV
} {}

test sv {cuteSV NA12878 minimap2} {
	cd $::smalltestdir
	file delete -force tmp/sv-cuteSV2_minimap2
	file mkdir tmp/sv-cuteSV2_minimap2
	mklink ori/ont/map-minimap2-NA12878_25FC_part19_21.bam tmp/sv-cuteSV2_minimap2/map-minimap2-NA12878_25FC_part19_21.bam
	mklink ori/ont/map-minimap2-NA12878_25FC_part19_21.bam.bai tmp/sv-cuteSV2_minimap2/map-minimap2-NA12878_25FC_part19_21.bam.bai
	cg sv_cuteSV {*}$::dopts \
		-refseq $::smalltestdir/refseqtest/hg19 \
		tmp/sv-cuteSV2_minimap2/map-minimap2-NA12878_25FC_part19_21.bam
	cg tsvdiff -q 1 -x *.tbi -x *.vcf \
		-ignorefields {varcaller_cg_version} \
		tmp/sv-cuteSV2_minimap2 expected/sv-cuteSV2_minimap2
} {}

test sv {npinv} {
	cd $::smalltestdir
	file delete -force tmp/sv-npinv
	file mkdir tmp/sv-npinv
	mklink ori/ont/bwa-mem_NA12878_25FC_part19_21.bam tmp/sv-npinv/bwa-mem_NA12878_25FC_part19_21.bam
	mklink ori/ont/bwa-mem_NA12878_25FC_part19_21.bam.bai tmp/sv-npinv/bwa-mem_NA12878_25FC_part19_21.bam.bai
	cg sv_npinv {*}$::dopts -refseq $::smalltestdir/refseqtest/hg19 tmp/sv-npinv/bwa-mem_NA12878_25FC_part19_21.bam
	cg tsvdiff -q 1 -x *.tbi -x *.vcf -x *.zsti \
		-ignorefields {varcaller_cg_version} \
		tmp/sv-npinv expected/sv-npinv
} {}

test sv {npinv -distrreg chr} {
	cd $::smalltestdir
	file delete -force tmp/sv-npinv
	file mkdir tmp/sv-npinv
	mklink ori/ont/bwa-mem_NA12878_25FC_part19_21.bam tmp/sv-npinv/bwa-mem_NA12878_25FC_part19_21.bam
	mklink ori/ont/bwa-mem_NA12878_25FC_part19_21.bam.bai tmp/sv-npinv/bwa-mem_NA12878_25FC_part19_21.bam.bai
	cg sv -method npinv {*}$::dopts \
		-refseq $::smalltestdir/refseqtest/hg19 \
		-distrreg chr \
		tmp/sv-npinv/bwa-mem_NA12878_25FC_part19_21.bam
	cg tsvdiff -q 1 -x *.tbi -x *.zsti \
		-ignorefields {varcaller_cg_version} \
		tmp/sv-npinv expected/sv-npinv
} {}

test sv {manta} {
	cd $::smalltestdir
	file delete -force tmp/sv-manta
	file mkdir tmp/sv-manta
	mklink ori/sv/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam tmp/sv-manta/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam
	mklink ori/sv/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam.bai tmp/sv-manta/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam.bai
	cg sv_manta {*}$::dopts \
		-refseq $::smalltestdir/refseqtest/hg19 \
		tmp/sv-manta/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam
	cg tsvdiff -q 1 -x *.xml -x svLocusGraphStats.tsv -x *.tbi -x *.py -x *.py.* -x alignmentStatsSummary.txt \
		-ignorefields {varcaller_cg_version} \
		tmp/sv-manta expected/sv-manta
} {}

#test sv {manta -regionfile} {
#	cd $::smalltestdir
#	file delete -force tmp/sv-manta-r
#	file mkdir tmp/sv-manta-r
#	mklink ori/sv/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam tmp/sv-manta-r/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam
#	mklink ori/sv/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam.bai tmp/sv-manta-r/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam.bai
#	file_write tmp/sv-manta-r/region.tsv "chromosome\tbegin\tend\nchr21\t21960600\t24444700\n"
#	exec cg sv_manta {*}$::dopts -regionfile tmp/sv-manta-r/region.tsv \
#		-refseq $::smalltestdir/refseqtest/hg19 \
#		tmp/sv-manta-r/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam >@ stdout 2>@ stderr
#	cg tsvdiff -q 1 -x *.xml -x svLocusGraphStats.tsv -x *.tbi -x *.py -x *.py.* -x alignmentStatsSummary.txt \
#		tmp/sv-manta-r expected/sv-manta-r
#} {}

test sv {cg sv -method manta, giving resultfile} {
	cd $::smalltestdir
	file delete -force tmp/sv-manta
	file mkdir tmp/sv-manta
	mklink ori/sv/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam tmp/sv-manta/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam
	mklink ori/sv/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam.bai tmp/sv-manta/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam.bai
	cg sv {*}$::dopts -method manta -refseq $::smalltestdir/refseqtest/hg19 tmp/sv-manta/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam tmp/sv-manta/resultsv.tsv
	cg zcat expected/sv-manta/sv-manta-dsbwa-ERR194147_30x_NA12878-chr21part.tsv.zst > tmp/sv-manta/expected.tsv
	cg tsvdiff -q 1 tmp/sv-manta/resultsv.tsv tmp/sv-manta/expected.tsv
} {}

# -distrreg is not actually used, because the -regionfile option in sv_manta is actually disabled (for now)
# because it causes too much problems, and is slower
test sv {cg sv -method manta -distreg 1} {
	cd $::smalltestdir
	file delete -force tmp/sv-manta
	file mkdir tmp/sv-manta
	mklink ori/sv/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam tmp/sv-manta/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam
	mklink ori/sv/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam.bai tmp/sv-manta/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam.bai
	exec cg sv -method manta -distrreg 1 {*}$::dopts -refseq $::smalltestdir/refseqtest/hg19 tmp/sv-manta/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam
	cg tsvdiff -q 1 -x *.xml -x svLocusGraphStats.tsv -x *.tbi -x *.py -x *.py.* -x alignmentStatsSummary.txt \
		-ignorefields {varcaller_cg_version} \
		tmp/sv-manta expected/sv-manta
} {}

#test sv {cg sv -method manta -distreg 1 -d 2} {
#	cd $::smalltestdir
#	file delete -force tmp/sv-manta
#	file mkdir tmp/sv-manta
#	mklink ori/sv/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam tmp/sv-manta/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam
#	mklink ori/sv/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam.bai tmp/sv-manta/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam.bai
#	exec cg sv -d 2 -method manta -distrreg 1 {*}$::dopts -refseq $::smalltestdir/refseqtest/hg19 tmp/sv-manta/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam
#	cg tsvdiff -q 1 -x *.xml -x svLocusGraphStats.tsv -x *.tbi -x *.py -x *.py.* -x alignmentStatsSummary.txt \
#		tmp/sv-manta expected/sv-manta
#} {}

test sv {lumpy} {
	cd $::smalltestdir
	file delete -force tmp/sv-lumpy
	file mkdir tmp/sv-lumpy
	mklink ori/sv/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam tmp/sv-lumpy/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam
	mklink ori/sv/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam.bai tmp/sv-lumpy/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam.bai
	cg sv_lumpy {*}$::dopts -refseq $::smalltestdir/refseqtest/hg19 tmp/sv-lumpy/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam
	cg tsvdiff -q 1 -x *.tbi \
		-ignorefields {varcaller_cg_version} \
		tmp/sv-lumpy expected/sv-lumpy
} {}

test sv {lumpy error} {
	cd $::smalltestdir
	file delete -force tmp/tmp
	file mkdir tmp/tmp
	file_write tmp/tmp/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam {}
	file_write tmp/tmp/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam.bai {}
	cg sv_lumpy {*}$::dopts -refseq $::smalltestdir/refseqtest/hg19 tmp/tmp/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam
	set result {}
} {*: ==: unary operator expected*} error match

#test sv {gridss} {
#	cd $::smalltestdir
#	file delete -force tmp/sv-gridss
#	file mkdir tmp/sv-gridss
#	mklink ori/sv/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam tmp/sv-gridss/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam
#	mklink ori/sv/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam.bai tmp/sv-gridss/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam.bai
#	cg sv_gridss {*}$::dopts -refseq $::smalltestdir/refseqtest/hg19 tmp/sv-gridss/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam
#	cg tsvdiff -q 1 -brief 1 -x *.idx -x *.tbi -x *.*_metrics -x *.pdf tmp/sv-gridss expected/sv-gridss
#} {}

testsummarize

