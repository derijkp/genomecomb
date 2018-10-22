#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl
set keepdir [pwd]

set test_cleantmp 0

# tests
# =====

test process_small {process_project mastr_mx2} {
	cd $::smalltestdir
	file delete -force tmp/mastr_mx2
	file mkdir tmp/mastr_mx2/samples
	foreach sample [glob ori/mastr_mx2.start/samples/*] {
		set sdest tmp/mastr_mx2/samples/[string range [file tail $sample] 0 end-3]/fastq
		file mkdir $sdest
		foreach file [glob -nocomplain $sample/ori/*] {
			file copy $file $sdest
		}
	}
	file copy -force ori/wgs2.mastr/samplicons-wgs2.tsv tmp/mastr_mx2/samplicons-wgs2.tsv
	# file copy ori/mastr_mx2/demultiplex_stats.tsv tmp/mastr_mx2
	# if you want to see output while running
	cg process_project {*}$::dopts -split 1 -reports -predictgender \
		-minfastqreads 10 -amplicons tmp/mastr_mx2/samplicons-wgs2.tsv -extra_reports_mastr 1 \
		tmp/mastr_mx2 $::refseqdir/hg19 >& tmp/mastr_mx2.log
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 -x *log_jobs -x *hsmetrics -x *.bam -x *.bai -x *.index -x fastqc_report.html \
		-x colinfo -x mastr_mx2.html -x *.lz4i -x *.finished -x info_analysis.tsv -x *.png \
		-x *.analysisinfo -x *.png \
		tmp/mastr_mx2 expected/mastr_mx2]
	lappend result [diffanalysisinfo tmp/mastr_mx2/compar/annot_compar-mastr_mx2.tsv.analysisinfo expected/mastr_mx2/compar/annot_compar-mastr_mx2.tsv.analysisinfo]
	foreach sample {blanco2_8485 ceph1333_02_34_7220 ceph1347_02_34_7149 ceph1347_02_34_8446} {
		lappend result [checkdiff -y --suppress-common-lines tmp/mastr_mx2/samples/$sample/reports/hsmetrics-crsbwa-$sample.hsmetrics expected/mastr_mx2/samples/$sample/reports/hsmetrics-crsbwa-$sample.hsmetrics | grep -v -E "Started on|net.sf.picard.analysis.directed.CalculateHsMetrics BAIT_INT"]
	}
	lappend result [checkdiff -y --suppress-common-lines tmp/mastr_mx2/mastr_mx2.html expected/mastr_mx2/mastr_mx2.html | grep -v -E {HistogramID|htmlwidget-|^<!|^<h2>20|meta charset|script.src *=}]
	lappend result [checkdiff -y --suppress-common-lines tmp/mastr_mx2/compar/info_analysis.tsv expected/mastr_mx2/compar/info_analysis.tsv | grep -v -E {version_os|param_adapterfile|param_targetvarsfile|param_dbfiles|command|version_genomecomb}]
	join [list_remove $result {}] \n
} {}

#test process_small {process_project mastr_mx2 space in name} {
#	cd $::smalltestdir
#	set mastrdir {tmp/mastr_mx2 space}
#	set expected {expected/mastr_mx2 space}
#	file delete -force $mastrdir
#	file mkdir $mastrdir/samples
#	foreach sample [glob ori/mastr_mx2.start/samples/*] {
#		set sdest $mastrdir/samples/[string range [file tail $sample] 0 end-3]/fastq
#		file mkdir $sdest
#		foreach file [glob -nocomplain $sample/ori/*] {
#			file copy $file $sdest
#		}
#	}
#	file copy -force ori/wgs2.mastr/samplicons-wgs2.tsv $mastrdir/samplicons-wgs2.tsv
#	# file copy ori/mastr_mx2/demultiplex_stats.tsv $mastrdir
#	# if you want to see output while running
#	cg process_project {*}$::dopts -split 1 \
#		-minfastqreads 10 -amplicons $mastrdir/samplicons-wgs2.tsv -extra_reports_mastr 1 \
#		$mastrdir $::refseqdir/hg19 >& tmp/mastr_mx2_space.log
#	# check vs expected
#	set result {}
#	lappend result [tsvdiff -q 1 -x *log_jobs -x *hsmetrics -x *.bam -x *.bai -x *.index -x fastqc_report.html \
#		-x colinfo -x "mastr_mx2 space.html" -x *.lz4i -x *.finished -x info_analysis.tsv -x *.png \
#		-x *.analysisinfo \
#		$mastrdir $expected
#	lappend result [diffanalysisinfo $mastrdir/compar/annot_compar-*.tsv.analysisinfo $expected/compar/annot_compar-*.tsv.analysisinfo]
#	foreach sample {blanco2_8485 ceph1333_02_34_7220 ceph1347_02_34_7149 ceph1347_02_34_8446} {
#		lappend result [checkdiff -y --suppress-common-lines $mastrdir/samples/$sample/reports/hsmetrics-crsbwa-$sample.hsmetrics $expected/samples/$sample/reports/hsmetrics-crsbwa-$sample.hsmetrics | grep -v -E "Started on|net.sf.picard.analysis.directed.CalculateHsMetrics BAIT_INT"]
#	}
#	lappend result [checkdiff -y --suppress-common-lines "$mastrdir/mastr_mx2 space.html" "$expected/mastr_mx2 space.html" | grep -v -E {HistogramID|htmlwidget-|^<!|^<h2>20|meta charset|script.src *=}]
#	lappend result [checkdiff -y --suppress-common-lines $mastrdir/compar/info_analysis.tsv $expected/compar/info_analysis.tsv | grep -v -E {version_os|param_adapterfile|param_targetvarsfile|param_dbfiles|command|version_genomecomb}]
#	join [list_remove $result {}] \n
#} {}

test process_small {process_project -jobsample 1 mastr_mx2_js1} {
	cd $::smalltestdir
	file delete -force tmp/mastr_mx2_js1
	file mkdir tmp/mastr_mx2_js1/samples
	foreach sample [glob ori/mastr_mx2.start/samples/*] {
		set sdest tmp/mastr_mx2_js1/samples/[string range [file tail $sample] 0 end-3]/fastq
		file mkdir $sdest
		foreach file [glob -nocomplain $sample/ori/*] {
			file copy $file $sdest
		}
	}
	file copy -force ori/wgs2.mastr/samplicons-wgs2.tsv tmp/mastr_mx2_js1/samplicons-wgs2.tsv
	cg process_project {*}$::dopts -jobsample 1 -split 1 \
		-minfastqreads 10 -amplicons tmp/mastr_mx2_js1/samplicons-wgs2.tsv -extra_reports_mastr 1 \
		tmp/mastr_mx2_js1 $::refseqdir/hg19 >& tmp/mastr_mx2_js1.log
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 -x *log_jobs -x *hsmetrics -x *.bam -x *.bai -x *.index -x fastqc_report.html \
		-x colinfo -x mastr_mx2_js1.html -x *.lz4i -x *.finished -x info_analysis.tsv -x *.png \
		-x *.analysisinfo -x *.png \
		tmp/mastr_mx2_js1 expected/mastr_mx2_js1]
	lappend result [diffanalysisinfo tmp/mastr_mx2_js1/compar/annot_compar-*.tsv.analysisinfo expected/mastr_mx2_js1/compar/annot_compar-*.tsv.analysisinfo]
	foreach sample {blanco2_8485 ceph1333_02_34_7220 ceph1347_02_34_7149 ceph1347_02_34_8446} {
		lappend result [checkdiff -y --suppress-common-lines tmp/mastr_mx2_js1/samples/$sample/reports/hsmetrics-crsbwa-$sample.hsmetrics expected/mastr_mx2_js1/samples/$sample/reports/hsmetrics-crsbwa-$sample.hsmetrics | grep -v -E "Started on|net.sf.picard.analysis.directed.CalculateHsMetrics BAIT_INT"]
	}
	lappend result [checkdiff -y --suppress-common-lines tmp/mastr_mx2_js1/mastr_mx2_js1.html expected/mastr_mx2_js1/mastr_mx2_js1.html | grep -v -E {HistogramID|htmlwidget-|^<!|^<h2>20|meta charset|script.src *=}]
	lappend result [checkdiff -y --suppress-common-lines tmp/mastr_mx2_js1/compar/info_analysis.tsv expected/mastr_mx2_js1/compar/info_analysis.tsv | grep -v -E {version_os|param_adapterfile|param_targetvarsfile|param_dbfiles|command|version_genomecomb}]
	join [list_remove $result {}] \n
} {}

test process_small {process_sample one_exome_yri_mx2} {
	cd $::smalltestdir
	file delete -force tmp/one_exome_yri_mx2
	file mkdir tmp/one_exome_yri_mx2/samples
	foreach sample {
		 NA19240mx2
	} {
		file mkdir tmp/one_exome_yri_mx2/samples/$sample/fastq
		file copy {*}[glob ori/exomes_yri_mx2.start/samples/$sample/ori/*.fq.gz] tmp/one_exome_yri_mx2/samples/$sample/fastq
		mklink $::refseqdir/hg19/extra/reg_hg19_exome_SeqCap_EZ_v3.tsv tmp/one_exome_yri_mx2/samples/$sample/reg_hg19_targets.tsv
	}
	cg process_sample {*}$::dopts -split 1 \
		-dbdir $::refseqdir/hg19 \
		tmp/one_exome_yri_mx2/samples/NA19240mx2 >& tmp/one_exome_yri_mx2.log
	# cg process_sample --stack 1 --verbose 2 -d status -split 1 -dbdir $::refseqdir/hg19 tmp/one_exome_yri_mx2/samples/NA19240mx2 | less -S
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x colinfo -x fastqc_report.html \
		-x *bam.dupmetrics -x info_analysis.tsv -x *.lz4i -x *.finished -x info_analysis.tsv -x *.png \
		-x *.analysisinfo -x *.png \
		tmp/one_exome_yri_mx2/samples/NA19240mx2 expected/one_exome_yri_mx2/samples/NA19240mx2]
	lappend result [checkdiff -y --suppress-common-lines tmp/one_exome_yri_mx2/samples/NA19240mx2/map-dsbwa-NA19240mx2.bam.dupmetrics expected/one_exome_yri_mx2/samples/NA19240mx2/map-dsbwa-NA19240mx2.bam.dupmetrics | grep -v "Started on" | grep -v "net.sf.picard.sam.MarkDuplicates INPUT" | grep -v bammarkduplicates2]
	lappend result [checkdiff -y --suppress-common-lines tmp/one_exome_yri_mx2/samples/NA19240mx2/info_analysis.tsv expected/one_exome_yri_mx2/samples/NA19240mx2/info_analysis.tsv | grep -v -E {version_os|param_adapterfile|param_targetvarsfile|param_dbfiles|command|version_genomecomb}]
	join [list_remove $result {}] \n
} {}

test process_small {process_project exomes_yri_mx2} {
	cd $::smalltestdir
	file delete -force tmp/exomes_yri_mx2
	file mkdir tmp/exomes_yri_mx2/samples
	foreach sample {
		NA19238mx2  NA19239mx2  NA19240mx2
	} {
		file mkdir tmp/exomes_yri_mx2/samples/$sample/fastq
		file copy {*}[glob ori/exomes_yri_mx2.start/samples/$sample/ori/*.fq.gz] tmp/exomes_yri_mx2/samples/$sample/fastq
		mklink $::refseqdir/hg19/extra/reg_hg19_exome_SeqCap_EZ_v3.tsv tmp/exomes_yri_mx2/samples/$sample/reg_hg19_targets.tsv
	}
	cg process_project {*}$::dopts -split 1 \
		-dbdir $::refseqdir/hg19 tmp/exomes_yri_mx2 >& tmp/exomes_yri_mx2.log
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x colinfo -x fastqc_report.html \
		-x *bam.dupmetrics -x info_analysis.tsv -x *.lz4i -x *.finished -x *.index -x info_analysis.tsv \
		-x *.analysisinfo -x *.png \
		tmp/exomes_yri_mx2 expected/exomes_yri_mx2]
	lappend result [diffanalysisinfo tmp/exomes_yri_mx2/compar/annot_compar-*.tsv.analysisinfo expected/exomes_yri_mx2/compar/annot_compar-*.tsv.analysisinfo]
	lappend result [checkdiff -y --suppress-common-lines tmp/exomes_yri_mx2/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics expected/exomes_yri_mx2/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics | grep -v "Started on" | grep -v bammarkduplicates2]
	foreach file1 [glob tmp/exomes_yri_mx2/compar/info_analysis.tsv tmp/exomes_yri_mx2/samples/*/info_analysis.tsv] {
		regsub ^tmp $file1 expected file2
		lappend result [checkdiff -y --suppress-common-lines $file1 $file2 | grep -v -E {version_os|param_adapterfile|param_targetvarsfile|param_dbfiles|command|version_genomecomb}]
	}
	join [list_remove $result {}] \n
} {}

test process_small {process_project exomesfb_yri_mx2 (freebayes)} {
	cd $::smalltestdir
	file delete -force tmp/exomesfb_yri_mx2
	file mkdir tmp/exomesfb_yri_mx2/samples
	foreach sample {
		NA19238mx2  NA19239mx2  NA19240mx2
	} {
		file mkdir tmp/exomesfb_yri_mx2/samples/$sample/fastq
		file copy {*}[glob ori/exomes_yri_mx2.start/samples/$sample/ori/*.fq.gz] tmp/exomesfb_yri_mx2/samples/$sample/fastq
		mklink $::refseqdir/hg19/extra/reg_hg19_exome_SeqCap_EZ_v3.tsv tmp/exomesfb_yri_mx2/samples/$sample/reg_hg19_targets.tsv
	}
	cg process_project {*}$::dopts -split 1 -varcallers {gatk freebayes} \
		-dbdir $::refseqdir/hg19 tmp/exomesfb_yri_mx2 >& tmp/exomesfb_yri_mx2.log
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x colinfo -x fastqc_report.html \
		-x *bam.dupmetrics -x info_analysis.tsv -x *.lz4i -x *.finished -x *.index -x info_analysis.tsv \
		-x *.analysisinfo -x *.png \
		tmp/exomesfb_yri_mx2 expected/exomesfb_yri_mx2]
	lappend result [diffanalysisinfo tmp/exomesfb_yri_mx2/compar/annot_compar-*.tsv.analysisinfo expected/exomesfb_yri_mx2/compar/annot_compar-*.tsv.analysisinfo]
	lappend result [checkdiff -y --suppress-common-lines tmp/exomesfb_yri_mx2/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics expected/exomesfb_yri_mx2/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics | grep -v "Started on" | grep -v bammarkduplicates2]
	foreach file1 [glob tmp/exomesfb_yri_mx2/compar/info_analysis.tsv tmp/exomesfb_yri_mx2/samples/*/info_analysis.tsv] {
		regsub ^tmp $file1 expected file2
		lappend result [checkdiff -y --suppress-common-lines $file1 $file2 | grep -v -E {version_os|param_adapterfile|param_targetvarsfile|param_dbfiles|command|version_genomecomb}]
	}
	join [list_remove $result {}] \n
} {}

test process_small {process_project exomes_gatkh_yri_mx2 (haplotypecaller)} {
	cd $::smalltestdir
	file delete -force tmp/exomes_gatkh_yri_mx2
	file mkdir tmp/exomes_gatkh_yri_mx2/samples
	foreach sample {
		NA19238mx2  NA19239mx2  NA19240mx2
	} {
		file mkdir tmp/exomes_gatkh_yri_mx2/samples/$sample/fastq
		file copy {*}[glob ori/exomes_yri_mx2.start/samples/$sample/ori/*.fq.gz] tmp/exomes_gatkh_yri_mx2/samples/$sample/fastq
		mklink $::refseqdir/hg19/extra/reg_hg19_exome_SeqCap_EZ_v3.tsv tmp/exomes_gatkh_yri_mx2/samples/$sample/reg_hg19_targets.tsv
	}
	cg process_project {*}$::dopts -split 1 -varcallers {gatkh freebayes} \
		-dbdir $::refseqdir/hg19 tmp/exomes_gatkh_yri_mx2 >& tmp/exomes_gatkh_yri_mx2.log
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x colinfo -x fastqc_report.html \
		-x *bam.dupmetrics -x info_analysis.tsv -x *.lz4i -x *.finished -x *.index -x info_analysis.tsv \
		-x *.analysisinfo -x *.png -x *.tbi \
		tmp/exomes_gatkh_yri_mx2 expected/exomes_gatkh_yri_mx2]
	lappend result [diffanalysisinfo tmp/exomes_gatkh_yri_mx2/compar/annot_compar-*.tsv.analysisinfo expected/exomes_gatkh_yri_mx2/compar/annot_compar-*.tsv.analysisinfo]
	lappend result [checkdiff -y --suppress-common-lines tmp/exomes_gatkh_yri_mx2/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics expected/exomes_gatkh_yri_mx2/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics | grep -v "Started on" | grep -v bammarkduplicates2]
	foreach file1 [glob tmp/exomes_gatkh_yri_mx2/compar/info_analysis.tsv tmp/exomes_gatkh_yri_mx2/samples/*/info_analysis.tsv] {
		regsub ^tmp $file1 expected file2
		lappend result [checkdiff -y --suppress-common-lines $file1 $file2 | grep -v -E {version_os|param_adapterfile|param_targetvarsfile|param_dbfiles|command|version_genomecomb}]
	}
	join [list_remove $result {}] \n
} {}

#test process_small {process_project exomes yri mx2} {
#	cd $::smalltestdir
#	file delete -force tmp/exomes_yri_mx2
#	file mkdir tmp/exomes_yri_mx2/samples
#	foreach sample {
#		NA19238mx2  NA19239mx2  NA19240mx2
#	} {
#		file mkdir tmp/exomes_yri_mx2/samples/$sample/fastq
#		file copy {*}[glob ori/exomes_yri_mx2.start/samples/$sample/ori/*.fq.gz] tmp/exomes_yri_mx2/samples/$sample/fastq
#	}
#	# cg process_illumina -d 2 -split 1 -dbdir $::refseqdir/hg19 tests/yri_exome
#	cg process_project --stack 1 --verbose 2 {*}$::dopts -split 1 -dbdir $::refseqdir/hg19 tmp/exomes_yri_mx2 >& tmp/exomes_yri_mx2.log
#	# check vs expected
#	checkdiff -y --suppress-common-lines tmp/exomes_yri_mx2/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics expected/exomes_yri_mx2/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics | grep -v "Started on" | grep -v bammarkduplicates2
#	checkdiff -qr -x *log_jobs -x *.bam -x *.bai -x colinfo -x *_fastqc -x *bam.dupmetrics tmp/exomes_yri_mx2 expected/exomes_yri_mx2
#	# could have used this, but previous is faster
#	# cg tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x colinfo -x *_fastqc -x *bam.dupmetrics tmp/exomes_yri_mx2 expected/exomes_yri_mx2 > temp
#} {}

test process_small {process_sample one_genome_yri_mx2} {
	cd $::smalltestdir
	set ref $::refseqdir/hg19
	file delete -force tmp/one_genome_yri_mx2
	file mkdir tmp/one_genome_yri_mx2/samples/NA19240cgmx2
	mklink ori/genomes_yri_mx2.start/samples/NA19240cgmx2/ori tmp/one_genome_yri_mx2/samples/NA19240cgmx2/ori
	cg process_sample {*}$::dopts -split 1 \
		-dbdir $ref tmp/one_genome_yri_mx2/samples/NA19240cgmx2 >& tmp/one_genome_yri_mx2.log
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 -x info_analysis.tsv -x log_jobs -x summary-NA19240cgmx2.txt -x *.finished -x info_analysis.tsv \
		-x *.analysisinfo -x *.png \
		tmp/one_genome_yri_mx2/samples/NA19240cgmx2 expected/genomes_yri_mx2/samples/NA19240cgmx2]
	join [list_remove $result {}] \n
} {}

test process_small {process_project genomes_yri_mx2} {
	cd $::smalltestdir
	set dest tmp/genomes_yri_mx2
	file delete -force tmp/genomes_yri_mx2
	file mkdir tmp/genomes_yri_mx2
	foreach sample {
		NA19238cgmx2 NA19239cgmx2 NA19240cgmx2
		NA19240ilmx2
	} {
		file mkdir tmp/genomes_yri_mx2/samples/$sample
		mklink ori/genomes_yri_mx2.start/samples/$sample/ori tmp/genomes_yri_mx2/samples/$sample/ori
	}
	# cg process_project --stack 1 --verbose 2 -d 2 -split 1 -dbdir /complgen/refseq/testdb2/hg19 tmp/genomes_yri_mx2
	cg process_project {*}$::dopts -split 1 \
		-dbdir $::refseqdir/hg19 tmp/genomes_yri_mx2 >& tmp/genomes_yri_mx2.log
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x *_fastqc -x summary-* -x fastqc_report.html \
		-x *dupmetrics -x colinfo -x *.lz4i -x info_analysis.tsv -x *.finished -x *.index -x info_analysis.tsv \
		-x *.analysisinfo -x *.png \
		tmp/genomes_yri_mx2 expected/genomes_yri_mx2]
	lappend result [diffanalysisinfo tmp/genomes_yri_mx2/compar/annot_compar-*.tsv.analysisinfo expected/genomes_yri_mx2/compar/annot_compar-*.tsv.analysisinfo]
	foreach cgsample {NA19238cgmx2 NA19239cgmx2 NA19240cgmx2} {
		lappend result [checkdiff -y --suppress-common-lines tmp/genomes_yri_mx2/samples/$cgsample/summary-$cgsample.txt expected/genomes_yri_mx2/samples/$cgsample/summary-$cgsample.txt | grep -v "finished.*finished"]
	}
	lappend result [checkdiff -y --suppress-common-lines tmp/genomes_yri_mx2/samples/NA19240ilmx2/map-dsbwa-NA19240ilmx2.bam.dupmetrics expected/genomes_yri_mx2/samples/NA19240ilmx2/map-dsbwa-NA19240ilmx2.bam.dupmetrics | grep -v "Started on" | grep -v bammarkduplicates2]
	foreach file1 [glob tmp/genomes_yri_mx2/compar/info_analysis.tsv tmp/genomes_yri_mx2/samples/*/info_analysis.tsv] {
		regsub ^tmp $file1 expected file2
		lappend result [checkdiff -y --suppress-common-lines $file1 $file2 | grep -v -E {version_os|param_adapterfile|param_targetvarsfile|param_dbfiles|command|version_genomecomb|Linux}]
	}
	join [list_remove $result {}] \n
} {}

test process_small {process_project mixed_yri_mx2} {
	cd $::smalltestdir
	set dest tmp/mixed_yri_mx2
	file delete -force tmp/mixed_yri_mx2
	file mkdir tmp/mixed_yri_mx2
	cg project_addsample tmp/mixed_yri_mx2 cgNA19240mx2 ori/mixed_yri_mx2/cgNA19240mx2
	cg project_addsample tmp/mixed_yri_mx2 gilNA19240mx2 {*}[glob ori/mixed_yri_mx2/gilNA19240mx2/*.fq.gz]
	cg project_addsample -targetfile ori/mixed_yri_mx2/reg_hg19_exome_SeqCap_EZ_v3.tsv.lz4 tmp/mixed_yri_mx2 exNA19239mx2 {*}[glob ori/mixed_yri_mx2/exNA19239mx2/*.fq.gz]
	cg project_addsample -targetfile ori/mixed_yri_mx2/reg_hg19_exome_SeqCap_EZ_v3.tsv.lz4  tmp/mixed_yri_mx2 exNA19240mx2 ori/mixed_yri_mx2/exNA19240mx2
	cg process_project {*}$::dopts -split 1 \
	  -dbdir /complgen/refseq/hg19 \
	  -dbfile /complgen/refseq/hg19/extra/var_hg19_dbnsfp.tsv.lz4 \
	  -dbfile /complgen/refseq/hg19/extra/var_hg19_dbnsfp.tsv.lz4 \
	  tmp/mixed_yri_mx2 >& tmp/mixed_yri_mx2.log
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x *_fastqc -x summary-* -x fastqc_report.html \
		-x *dupmetrics -x colinfo -x *.lz4i -x info_analysis.tsv -x *.finished -x *.index \
		-x *.analysisinfo -x *.png \
		tmp/mixed_yri_mx2 expected/mixed_yri_mx2]
	foreach cgsample {NA19238cgmx2 NA19239cgmx2 NA19240cgmx2} {
		lappend result [checkdiff -y --suppress-common-lines tmp/genomes_yri_mx2/samples/$cgsample/summary-$cgsample.txt expected/genomes_yri_mx2/samples/$cgsample/summary-$cgsample.txt | grep -v "finished.*finished"]
	}
	lappend result [diffanalysisinfo tmp/mixed_yri_mx2/compar/annot_compar-*.tsv.analysisinfo expected/mixed_yri_mx2/compar/annot_compar-*.tsv.analysisinfo]
	lappend result [checkdiff -y --suppress-common-lines tmp/mixed_yri_mx2/samples/gilNA19240mx2/map-dsbwa-gilNA19240mx2.bam.dupmetrics expected/mixed_yri_mx2/samples/gilNA19240mx2/map-dsbwa-gilNA19240mx2.bam.dupmetrics | grep -v "Started on" | grep -v bammarkduplicates2]
	foreach file1 [glob tmp/genomes_yri_mx2/compar/info_analysis.tsv tmp/genomes_yri_mx2/samples/*/info_analysis.tsv] {
		regsub ^tmp $file1 expected file2
		lappend result [checkdiff -y --suppress-common-lines $file1 $file2 | grep -v -E {version_os|param_adapterfile|param_targetvarsfile|param_dbfiles|command|version_genomecomb}]
	}
	join [list_remove $result {}] \n
} {}

test process_small {process_project -distrreg 1 mixed_yri_mx2_distrreg} {
	cd $::smalltestdir
	set dest tmp/mixed_yri_mx2_distrreg
	file delete -force tmp/mixed_yri_mx2_distrreg
	file mkdir tmp/mixed_yri_mx2_distrreg
	cg project_addsample tmp/mixed_yri_mx2_distrreg cgNA19240mx2 ori/mixed_yri_mx2/cgNA19240mx2
	cg project_addsample tmp/mixed_yri_mx2_distrreg gilNA19240mx2 {*}[glob ori/mixed_yri_mx2/gilNA19240mx2/*.fq.gz]
	cg project_addsample -targetfile ori/mixed_yri_mx2/reg_hg19_exome_SeqCap_EZ_v3.tsv.lz4 tmp/mixed_yri_mx2_distrreg exNA19239mx2 {*}[glob ori/mixed_yri_mx2/exNA19239mx2/*.fq.gz]
	cg project_addsample -targetfile ori/mixed_yri_mx2/reg_hg19_exome_SeqCap_EZ_v3.tsv.lz4 tmp/mixed_yri_mx2_distrreg exNA19240mx2 ori/mixed_yri_mx2/exNA19240mx2
	cg process_project {*}$::dopts -distrreg 1 -split 1 \
	  -dbdir /complgen/refseq/hg19 \
	  -dbfile /complgen/refseq/hg19/extra/var_hg19_dbnsfp.tsv.lz4 \
	  -dbfile /complgen/refseq/hg19/extra/var_hg19_dbnsfp.tsv.lz4 \
	  tmp/mixed_yri_mx2_distrreg >& tmp/mixed_yri_mx2_distrreg.log
	# check vs expected
	foreach cgsample {NA19238cgmx2 NA19239cgmx2 NA19240cgmx2} {
		checkdiff -y --suppress-common-lines tmp/genomes_yri_mx2/samples/$cgsample/summary-$cgsample.txt expected/genomes_yri_mx2/samples/$cgsample/summary-$cgsample.txt | grep -v "finished.*finished"
	}
	set result {}
	lappend result [tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x *_fastqc -x summary-* -x fastqc_report.html \
		-x *dupmetrics -x colinfo -x *.lz4i -x info_analysis.tsv -x *.finished -x *.index \
		-x *.analysisinfo -x *.png \
		tmp/mixed_yri_mx2_distrreg expected/mixed_yri_mx2_distrreg]
	lappend result [diffanalysisinfo tmp/mixed_yri_mx2_distrreg/compar/annot_compar-*.tsv.analysisinfo expected/mixed_yri_mx2_distrreg/compar/annot_compar-*.tsv.analysisinfo]
	lappend result [checkdiff -y --suppress-common-lines tmp/mixed_yri_mx2_distrreg/samples/gilNA19240mx2/map-dsbwa-gilNA19240mx2.bam.dupmetrics expected/mixed_yri_mx2_distrreg/samples/gilNA19240mx2/map-dsbwa-gilNA19240mx2.bam.dupmetrics | grep -v "Started on" | grep -v bammarkduplicates2]
	foreach file1 [glob tmp/genomes_yri_mx2/compar/info_analysis.tsv tmp/genomes_yri_mx2/samples/*/info_analysis.tsv] {
		regsub ^tmp $file1 expected file2
		lappend result [checkdiff -y --suppress-common-lines $file1 $file2 | grep -v -E {version_os|param_adapterfile|param_targetvarsfile|param_dbfiles|command|version_genomecomb}]
	}
	join [list_remove $result {}] \n
} {}

if 0 {

test process_small {annotate refseqbuild/hg19} {
	cd $::smalltestdir
	file delete tmp/annot_exomes_yri_mx2.tsv
	cg annotate --stack 1 --verbose 2 expected/exomes_yri_mx2/compar/compar-exomes_yri_mx2.tsv tmp/annot_exomes_yri_mx2.tsv /data/genomecomb.testdata/refseqbuild/hg19 /data/genomecomb.testdata/refseqbuild/hg19/extra >& tmp/annot_exomes_yri_mx2.log
	checkdiff tmp/annot_exomes_yri_mx2.tsv expected/annot_exomes_yri_mx2.tsv
} {}

test process_small {annotate refseqbuild/hg38} {
	cd $::smalltestdir
	# compar-exomes_yri_mx2.tsv is actually hg19, but we use it here anyway just to test if the hg38 works
	file delete tmp/annothg38_exomes_yri_mx2.tsv
	cg annotate --stack 1 --verbose 2 expected/exomes_yri_mx2/compar/compar-exomes_yri_mx2.tsv tmp/annothg38_exomes_yri_mx2.tsv /data/genomecomb.testdata/refseqbuild/hg38 /data/genomecomb.testdata/refseqbuild/hg38/extra >& tmp/annothg38_exomes_yri_mx2.log
	chekdiff tmp/annothg38_exomes_yri_mx2.tsv expected/annothg38_exomes_yri_mx2.tsv
}

}

testsummarize

