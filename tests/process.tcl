#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl
set keepdir [pwd]

set test_cleantmp 0

if 0 {
	# make exome testdata (done using seqcap v2?)
	file mkdir /data/testdata/ori/exomes_yri.start/samples/NA19238/fastq
	cd /data/testdata/ori/exomes_yri.start/samples/NA19238/fastq
	exec wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19238/sequence_read/SRR071173_2.filt.fastq.gz 2>@ stderr
	exec wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19238/sequence_read/SRR071173_1.filt.fastq.gz 2>@ stderr
	cd .. ; 	exec echo seqcapv3 > info_capture.txt
	file mkdir /data/testdata/ori/exomes_yri.start/samples/NA19239/fastq
	cd /data/testdata/ori/exomes_yri.start/samples/NA19239/fastq
	exec wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19239/sequence_read/SRR792097_1.filt.fastq.gz 2>@ stderr
	exec wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19239/sequence_read/SRR792097_2.filt.fastq.gz 2>@ stderr
	cd .. ; exec ln -s /complgen/refseq/hg19/extra/reg_hg19_exome_SeqCap_EZ_v3.tsv reg_hg19_targets.tsv
	file mkdir /data/testdata/ori/exomes_yri.start/samples/NA19240/fastq
	cd /data/testdata/ori/exomes_yri.start/samples/NA19240/fastq
	exec wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19240/sequence_read/SRR792091_1.filt.fastq.gz 2>@ stderr
	exec wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19240/sequence_read/SRR792091_2.filt.fastq.gz 2>@ stderr
	cd .. ; exec ln -s /complgen/refseq/hg19/extra/reg_hg19_exome_SeqCap_EZ_v3.tsv reg_hg19_targets.tsv

}

test process {mastr mastr_116068_116083} {
	cd $::bigtestdir	
	file delete -force tmp/mastr_116068_116083
#	file copy ori/mastr_116068_116083 tmp/
	file mkdir tmp/mastr_116068_116083
	file copy {*}[glob ori/mastr_116068_116083/samples/*] tmp/mastr_116068_116083
	foreach dir [glob tmp/mastr_116068_116083/*] {
		file mkdir $dir/fastq
		foreach file [glob -nocomplain $dir/ori/*] {
			file copy $file $dir/fastq
		}
	}
	file delete -force tmp/wgs2.mastr
	file mkdir tmp/wgs2.mastr
	file copy -force ori/wgs2.mastr/amplicons-wgs2.tsv tmp/wgs2.mastr
	# file copy ori/mastr_116068_116083/demultiplex_stats.tsv tmp/mastr_116068_116083
	# if you want to see output while running
	 cg process_mastr {*}$::dopts -split 1 tmp/wgs2.mastr tmp/mastr_116068_116083 refseqtest/hg19 2>@ stderr >@ stdout
	# no output while running
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x *hsmetrics -x colinfo -x *.index -x *.zsti -x *.lz4i  \
		-x mastr_116068_116083.html -x *.finished \
		tmp/mastr_116068_116083 expected/mastr_116068_116083]
	lappend result [checkdiff -y --suppress-common-lines tmp/mastr_116068_116083/mastr_116068_116083.html expected/mastr_116068_116083/mastr_116068_116083.html | grep -v -E {HistogramID|htmlwidget-|^<!|^<h2>20}]
	foreach sample [dirglob tmp/mastr_116068_116083 ceph*] {
		lappend result [checkdiff -y --suppress-common-lines tmp/mastr_116068_116083/$sample/crsbwa-$sample.hsmetrics expected/mastr_116068_116083/$sample/crsbwa-$sample.hsmetrics | grep -v -E "Started on|net.sf.picard.analysis.directed.CalculateHsMetrics BAIT_INT"]
	}
	join [list_remove $result {}] \n
} {}

test process {process_project illumina exomes yri chr2122} {
	cd $::bigtestdir
	file delete -force tmp/exomes_yri_chr2122
	file mkdir tmp/exomes_yri_chr2122/samples
	foreach sample {
		NA19238chr2122  NA19239chr2122  NA19240chr2122
	} {
		cg project_addsample tmp/exomes_yri_chr2122 $sample {*}[glob ori/exomes_yri_chr2122.start/samples/$sample/fastq/*]
	}
	if 0 {
		foreach sample {NA19238chr2122 NA19239chr2122 NA19240chr2122} {
			exec cp -a expected/exomes_yri_chr2122/samples/${sample}/map-rdsbwa-${sample}.bam tmp/exomes_yri_chr2122/samples/${sample}/
			exec cp -a expected/exomes_yri_chr2122/samples/${sample}/map-rdsbwa-${sample}.bam.bai tmp/exomes_yri_chr2122/samples/${sample}/
		}
		exec touch {*}[glob tmp/exomes_yri_chr2122/samples/*/map-*.bam*]
	}
	cg process_project {*}$::dopts \
	  -targetfile ori/mixed_yri_mx2/reg_hg19_exome_SeqCap_EZ_v3.tsv.zst \
	  -split 1 -dbdir refseqtest/hg19 tmp/exomes_yri_chr2122 2>@ stderr >@ stdout
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x colinfo -x fastqc_report.html \
		-x *bam.dupmetrics -x info_analysis.tsv -x *.zsti -x *.lz4i -x *.finished -x *.index \
		tmp/exomes_yri_chr2122 expected/exomes_yri_chr2122]
	lappend result [checkdiff -y --suppress-common-lines tmp/exomes_yri_chr2122/samples/NA19238chr2122/map-dsbwa-NA19238chr2122.bam.dupmetrics expected/exomes_yri_chr2122/samples/NA19238chr2122/map-dsbwa-NA19238chr2122.bam.dupmetrics | grep -v "Started on"]
	foreach file1 [glob tmp/exomes_yri_chr2122/compar/info_analysis.tsv tmp/exomes_yri_chr2122/samples/*/info_analysis.tsv] {
		regsub ^tmp $file1 expected file2
		lappend result [checkdiff -y --suppress-common-lines $file1 $file2 | grep -v -E {version_os}]
	}
	lappend result [checkdiff -y --suppress-common-lines tmp/exomes_yri_chr2122/samples/NA19238chr2122/map-dsbwa-NA19238chr2122.bam.dupmetrics expected/exomes_yri_chr2122/samples/NA19238chr2122/map-dsbwa-NA19238chr2122.bam.dupmetrics | grep -v "Started on"]
	join [list_remove $result {}] \n
} {}

test process {process_project exomes_gatkh_strelka_yri_chr2122 (haplotypecaller + strelka)} {
	cd $::bigtestdir
	set testdir tmp/exomes_gatkh_strelka_yri_chr2122
	set src ori/exomes_yri_chr2122.start/
	file delete -force $testdir
	file mkdir $testdir/samples
	foreach sample {
		NA19238chr2122  NA19239chr2122  NA19240chr2122
	} {
		cg project_addsample $testdir $sample {*}[glob $src/samples/$sample/fastq/*]
	}
	cg process_project {*}$::dopts -split 1 -varcallers {gatkh strelka} \
		-dbdir refseqtest/hg19 $testdir >& $testdir.log
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x colinfo -x fastqc_report.html \
		-x *bam.dupmetrics -x info_analysis.tsv -x *.zsti -x *.lz4i -x *.finished -x *.index -x info_analysis.tsv \
		-x *.analysisinfo -x *.png \
		tmp/exomes_gatkh_strelka_yri_chr2122 expected/exomes_gatkh_strelka_yri_chr2122]
	lappend result [diffanalysisinfo tmp/exomes_gatkh_strelka_yri_chr2122/compar/annot_compar-*.tsv.analysisinfo expected/exomes_gatkh_strelka_yri_chr2122/compar/annot_compar-*.tsv.analysisinfo]
	checkdiff -y --suppress-common-lines tmp/exomes_gatkh_strelka_yri_chr2122/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics expected/exomes_gatkh_strelka_yri_chr2122/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics | grep -v "Started on" | grep -v bammarkduplicates2
	foreach file1 [glob tmp/exomes_gatkh_strelka_yri_chr2122/compar/info_analysis.tsv tmp/exomes_gatkh_strelka_yri_chr2122/samples/*/info_analysis.tsv] {
		regsub ^tmp $file1 expected file2
		lappend result [checkdiff -y --suppress-common-lines $file1 $file2 | grep -v -E {version_os|param_adapterfile|param_targetvarsfile|param_dbfiles|command|version_genomecomb}]
	}
	join [list_remove $result {}] \n
} {}

test process {process_project exomes_gatkh_yri_chr2122 (haplotypecaller)} {
	cd $::bigtestdir
	set testdir tmp/exomes_gatkh_yri_chr2122
	set src ori/exomes_yri_chr2122.start/
	file delete -force $testdir
	file mkdir $testdir/samples
	foreach sample {
		NA19238chr2122  NA19239chr2122  NA19240chr2122
	} {
		cg project_addsample $testdir $sample {*}[glob $src/samples/$sample/fastq/*]
	}
	cg process_project {*}$::dopts -split 1 -varcallers {gatkh freebayes} \
		-dbdir refseqtest/hg19 $testdir >& $testdir.log
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x colinfo -x fastqc_report.html \
		-x *bam.dupmetrics -x info_analysis.tsv -x *.zsti -x *.finished -x *.index -x info_analysis.tsv \
		-x *.analysisinfo -x *.png \
		tmp/exomes_gatkh_yri_chr2122 expected/exomes_gatkh_yri_chr2122]
	lappend result [diffanalysisinfo tmp/exomes_gatkh_yri_chr2122/compar/annot_compar-*.tsv.analysisinfo expected/exomes_gatkh_yri_chr2122/compar/annot_compar-*.tsv.analysisinfo]
	checkdiff -y --suppress-common-lines tmp/exomes_gatkh_yri_chr2122/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics expected/exomes_gatkh_yri_chr2122/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics | grep -v "Started on" | grep -v bammarkduplicates2
	foreach file1 [glob tmp/exomes_gatkh_yri_chr2122/compar/info_analysis.tsv tmp/exomes_gatkh_yri_chr2122/samples/*/info_analysis.tsv] {
		regsub ^tmp $file1 expected file2
		lappend result [checkdiff -y --suppress-common-lines $file1 $file2 | grep -v -E {version_os|param_adapterfile|param_targetvarsfile|param_dbfiles|command|version_genomecomb}]
	}
	join [list_remove $result {}] \n
} {}

test process {genomes yri chr2122} {
	cd $::bigtestdir	
	file delete -force tmp/genomes_yri_chr2122
	file mkdir tmp/genomes_yri_chr2122
	cg project_addsample tmp/genomes_yri_chr2122 testNA19238chr2122cg ori/genomes_yri_chr2122.start/samples/testNA19238chr2122cg.ori
	cg project_addsample tmp/genomes_yri_chr2122 testNA19239chr2122cg ori/genomes_yri_chr2122.start/samples/testNA19239chr2122cg.ori
	cg project_addsample tmp/genomes_yri_chr2122 testNA19240chr2122cg ori/genomes_yri_chr2122.start/samples/testNA19240chr2122cg.ori
	cg project_addsample tmp/genomes_yri_chr2122 testNA19240chr21il {*}[glob ori/genomes_yri_chr2122.start/samples/testNA19240chr21il.ori/*.fq.gz]
	cg process_project {*}$::dopts -split 1 -dbdir refseqtest/hg19 tmp/genomes_yri_chr2122 2>@ stderr >@ stdout
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x *_fastqc -x summary-* -x fastqc_report.html \
		-x *dupmetrics -x colinfo -x *.zsti -x info_analysis.tsv -x *.finished -x *.index \
		tmp/genomes_yri_chr2122 expected/genomes_yri_chr2122]
	lappend result [checkdiff -y --suppress-common-lines tmp/genomes_yri_chr2122/samples/testNA19240chr21il/map-dsbwa-testNA19240chr21il.bam.dupmetrics expected/genomes_yri_chr2122/samples/testNA19240chr21il/map-dsbwa-testNA19240chr21il.bam.dupmetrics | grep -v "Started on"]
	foreach file1 [glob tmp/genomes_yri_chr2122/compar/info_analysis.tsv tmp/genomes_yri_chr2122/samples/*/info_analysis.tsv] {
		regsub ^tmp $file1 expected file2
		lappend result [checkdiff -y --suppress-common-lines $file1 $file2 | grep -v -E {version_os}]
	}
	# file_write temp $e
	join [list_remove $result {}] \n
} {}

cd $keepdir

testsummarize
