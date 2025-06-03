#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl
set keepdir [pwd]

set test_cleantmp 0

# gatk is not deterministic with threads, so no threads for testing and comparing to previous runs
lappend dopts -threads 1

if 0 {
	cd /data/genomes
	rs peterdr@crema:/complgen3/peterdr/projects/pubgenomes_tests.data/platinum/samples/ERR194147_30x_NA12878/map-dsbwa-ERR194147_30x_NA12878.bam .
	# create illumina test data (from platinum pubgenomes originally)
	samtools view --no-PG -hb /data/genomes/map-dsbwa-ERR194147_30x_NA12878.bam chr21 > /data/genomes/map-dsbwa-ERR194147_30x_NA12878-chr21.bam
	samtools index /data/genomes/map-dsbwa-ERR194147_30x_NA12878-chr21.bam
	samtools view --no-PG -hb /data/genomes/map-dsbwa-ERR194147_30x_NA12878-chr21.bam chr21:21800000-25100000 > /data/genomes/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam
	samtools index /data/genomes/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam
	cg bam2fastq /data/genomes/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam /data/genomes/ERR194147_30x_NA12878-chr21part_1.fq.gz /data/genomes/ERR194147_30x_NA12878-chr21part_2.fq.gz

	rs peterdr@crema:/complgen3/peterdr/projects/pubgenomes_tests.data/platinum/samples/ERR194146_30x_NA12877/map-dsbwa-ERR194146_30x_NA12877.bam* .
	# create illumina test data (from platinum pubgenomes originally)
	samtools view --no-PG -hb /data/genomes/map-dsbwa-ERR194146_30x_NA12877.bam chr21:21800000-25100000 > /data/genomes/map-dsbwa-ERR194146_30x_NA12877-chr21part.bam
	samtools index /data/genomes/map-dsbwa-ERR194146_30x_NA12877-chr21part.bam
	cg bam2fastq /data/genomes/map-dsbwa-ERR194146_30x_NA12877-chr21part.bam /data/genomes/ERR194146_30x_NA12877-chr21part_1.fq.gz /data/genomes/ERR194146_30x_NA12877-chr21part_2.fq.gz

	cp /data/genomes/*-chr21part* /data/genomecomb.testdata/ori/sv/
		
}


# tests
# =====

test process_sv {process_project ont} {
	cd $::smalltestdir
	set dest tmp/ont
	file delete -force tmp/ont
	file mkdir tmp/ont
	cg project_addsample tmp/ont NA12878 {*}[glob ori/nanopore/NA12878-nanopore-wgs/part*.fastq*]
	cg process_project {*}$::dopts -distrreg 1 -split 1 \
	  -dbdir /complgen/refseq/hg19 -reports {-fastqc predictgender fastqstats} \
	  -clip 0 -paired 0 -aligner ngmlr -removeduplicates 0 -realign 0 -svcallers {sniffles cuteSV cuteSV_pacbio} -varcallers {} \
	  tmp/ont >& tmp/ont.log
	set result {}
	lappend result [tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x *.zsti -x *_fastqc -x summary-* -x fastqc_report.html \
		-x *dupmetrics -x colinfo -x *.lz4i -x info_analysis.tsv -x *.finished -x *.index \
		-x *.analysisinfo -x *.png -x *.vcf -x *.submitting \
		-x info_analysis.tsv \
		-x *.snf \
		tmp/ont expected/ont]
	foreach file1 [glob tmp/ont/compar/info_analysis.tsv tmp/genomes_yri_mx2/samples/*/info_analysis.tsv] {
		regsub ^tmp $file1 expected file2
		lappend result [checkdiff -y --suppress-common-lines $file1 $file2 | grep -v -E {version_os|param_adapterfile|param_targetvarsfile|param_dbfiles|command|version_genomecomb}]
	}
	join [list_remove $result {}] \n
} {}

test process_sv {process_project ont_minimap2 and -extraannot AnnotSV} {
	cd $::smalltestdir
	set dest tmp/ont_minimap2
	file delete -force tmp/ont_minimap2
	file mkdir tmp/ont_minimap2
	cg project_addsample tmp/ont_minimap2 NA12878 {*}[glob ori/nanopore/NA12878-nanopore-wgs/part*.fastq*]
	cg process_project {*}$::dopts -distrreg 1 -split 1 \
	  -dbdir $::refseqdir/hg19 -reports {-fastqc predictgender fastqstats} \
	  -clip 0 -paired 0 -aligner minimap2 -removeduplicates 0 -realign 0 -svcallers {sniffles cuteSV cuteSV_pacbio} -varcallers {} \
	  -extraannot AnnotSV \
	  tmp/ont_minimap2 >& tmp/ont_minimap2.log
	set result {}
	lappend result [tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x *.zsti -x *_fastqc -x summary-* -x fastqc_report.html \
		-x *dupmetrics -x colinfo -x *.lz4i -x info_analysis.tsv -x *.finished -x *.index \
		-x *.tbi -x *.submitting -x *.xml \
		-x *.png -x *.vcf \
		-x *.snf \
		tmp/ont_minimap2 expected/ont_minimap2]
	foreach file1 [glob tmp/ont_minimap2/compar/info_analysis.tsv tmp/genomes_yri_mx2/samples/*/info_analysis.tsv] {
		regsub ^tmp $file1 expected file2
		lappend result [checkdiff -y --suppress-common-lines $file1 $file2 | grep -v -E {version_os|param_adapterfile|param_targetvarsfile|param_dbfiles|command|version_genomecomb|param_::maxopenfiles}]
	}
	join [list_remove $result {}] \n
} {}

test process_sv {manta} {
	cd $::smalltestdir
	file delete -force tmp/sv_chr21part
	file mkdir tmp/sv_chr21part
	cg project_addsample tmp/sv_chr21part ERR194146_30x_NA12877 {*}[glob ori/sv/ERR194146_30x_NA12877-chr21part_*]
	cg project_addsample tmp/sv_chr21part ERR194147_30x_NA12878 {*}[glob ori/sv/ERR194147_30x_NA12878-chr21part_*]
	mklink $::smalltestdir/ori/sv/map-rdsbwa-ERR194146_30x_NA12877-chr21part.bam tmp/sv_chr21part/samples/ERR194146_30x_NA12877/map-rdsbwa-ERR194146_30x_NA12877.bam
	mklink $::smalltestdir/ori/sv/map-rdsbwa-ERR194147_30x_NA12878-chr21part.bam tmp/sv_chr21part/samples/ERR194147_30x_NA12878/map-rdsbwa-ERR194147_30x_NA12878.bam
	cg process_project {*}$::dopts \
	  -varcallers gatk \
	  -svcallers manta -distrreg 1 -split 1 \
	  -dbdir $::refseqdir/hg19 \
	  tmp/sv_chr21part >& tmp/sv_chr21part.log
	set result {}
	lappend result [tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x *.zsti -x *_fastqc -x summary-* -x fastqc_report.html \
		-x *.tbi -x *.submitting -x *.xml \
		-x *dupmetrics -x colinfo -x *.lz4i -x *.zsti -x info_analysis.tsv -x *.finished -x *.index \
		-x *.analysisinfo -x *.png -x *.vcf \
		-x svLocusGraphStats.tsv \
		tmp/sv_chr21part expected/sv_chr21part]
	join [list_remove $result {}] \n
} {}

testsummarize


