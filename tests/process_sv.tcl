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

if 0 {
	cd /data/genomes
	rs peterdr@crema:/complgen3/peterdr/projects/pubgenomes_tests.data/platinum/samples/ERR194147_30x_NA12878/map-dsbwa-ERR194147_30x_NA12878.bam .
	# create illumina test data (from platinum pubgenomes originally)
	samtools view -hb /data/genomes/map-dsbwa-ERR194147_30x_NA12878.bam chr21 > /data/genomes/map-dsbwa-ERR194147_30x_NA12878-chr21.bam
	samtools index /data/genomes/map-dsbwa-ERR194147_30x_NA12878-chr21.bam
	samtools view -hb /data/genomes/map-dsbwa-ERR194147_30x_NA12878-chr21.bam chr21:21800000-25100000 > /data/genomes/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam
	samtools index /data/genomes/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam
	cg bam2fastq /data/genomes/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam /data/genomes/ERR194147_30x_NA12878-chr21part_1.fq.gz /data/genomes/ERR194147_30x_NA12878-chr21part_2.fq.gz

	rs peterdr@crema:/complgen3/peterdr/projects/pubgenomes_tests.data/platinum/samples/ERR194146_30x_NA12877/map-dsbwa-ERR194146_30x_NA12877.bam* .
	# create illumina test data (from platinum pubgenomes originally)
	samtools view -hb /data/genomes/map-dsbwa-ERR194146_30x_NA12877.bam chr21:21800000-25100000 > /data/genomes/map-dsbwa-ERR194146_30x_NA12877-chr21part.bam
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
	cg project_addsample tmp/ont NA12878 {*}[glob /data/nanopore/NA12878-nanopore-wgs/part*.fastq*]
	cg process_project {*}$::dopts -distrreg 1 -split 1 \
	  -dbdir /complgen/refseq/hg19 -reports {-fastqc predictgender fastqstats} \
	  -clip 0 -paired 0 -aligner ngmlr -removeduplicates 0 -realign 0 -svcallers sniffles -varcallers {} \
	  tmp/ont >& tmp/ont.log
	cg tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x *_fastqc -x summary-* -x fastqc_report.html \
		-x *dupmetrics -x colinfo -x *.lz4i -x info_analysis.tsv -x *.finished -x *.index \
		-x *.analysisinfo -x *.png -x *.vcf \
		tmp/ont expected/ont
	set errors {}
	foreach file1 [glob tmp/ont/compar/info_analysis.tsv tmp/genomes_yri_mx2/samples/*/info_analysis.tsv] {
		regsub ^tmp $file1 expected file2
		if {[catch {
			checkdiff -y --suppress-common-lines $file1 $file2 | grep -v -E {version_os|param_adapterfile|param_targetvarsfile|param_dbfiles|command|version_genomecomb}
		} msg]} {
			lappend errors $file1 $msg
		}
	}
	join $errors \n
} {}

test process_sv {manta} {
	cd $::smalltestdir
	file delete -force tmp/sv_chr21part
	file mkdir tmp/sv_chr21part
	cg project_addsample tmp/sv_chr21part ERR194146_30x_NA12877 {*}[glob ori/sv/ERR194146_30x_NA12877-chr21part_*]
	cg project_addsample tmp/sv_chr21part ERR194147_30x_NA12878 {*}[glob ori/sv/ERR194147_30x_NA12878-chr21part_*]
	mklink $::smalltestdir/ori/sv/map-dsbwa-ERR194146_30x_NA12877-chr21part.bam tmp/sv_chr21part/samples/ERR194146_30x_NA12877/map-dsbwa-ERR194146_30x_NA12877-chr21part.bam
	mklink $::smalltestdir/ori/sv/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam tmp/sv_chr21part/samples/ERR194147_30x_NA12878/map-dsbwa-ERR194147_30x_NA12878-chr21part.bam
	cg process_project {*}$::dopts \
	  -varcallers gatk \
	  -svcallers manta -distrreg 1 -split 1 \
	  -dbdir /complgen/refseq/hg19 \
	  -dbfile /complgen/refseq/hg19/extra/var_hg19_dbnsfp.tsv.lz4 \
	  -dbfile /complgen/refseq/hg19/extra/var_hg19_dbnsfp.tsv.lz4 \
	  tmp/sv_chr21part >& tmp/sv_chr21part.log
	cg tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x *_fastqc -x summary-* -x fastqc_report.html \
		-x *dupmetrics -x colinfo -x *.lz4i -x info_analysis.tsv -x *.finished -x *.index \
		-x *.analysisinfo -x *.png -x *.vcf \
		tmp/sv_chr21part expected/sv_chr21part
} {}

testsummarize

