#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl
set keepdir [pwd]

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

test process {process_sample exome yri chr2122} {
	cd $::bigtestdir
	file delete -force tmp/one_exome_yri_chr2122
	file mkdir tmp/one_exome_yri_chr2122/samples
	foreach sample {
		NA19240chr2122
	} {
		file mkdir tmp/one_exome_yri_chr2122/samples/$sample/fastq
		foreach file [glob ori/exomes_yri_chr2122.start/samples/$sample/fastq/*] {
			file copy $file tmp/one_exome_yri_chr2122/samples/$sample/fastq/[file tail $file]
		}
	}
	cg process_sample --stack 1 --verbose 2 -split 1 -dbdir refseqtest/hg19 tmp/one_exome_yri_chr2122/samples/NA19240chr2122 2>@ stderr >@ stdout
	# check vs expected
	checkdiff -y --suppress-common-lines tmp/one_exome_yri_chr2122/samples/NA19240chr2122/map-dsbwa-NA19240chr2122.bam.dupmetrics expected/one_exome_yri_chr2122/samples/NA19240chr2122/map-dsbwa-NA19240chr2122.bam.dupmetrics | grep -v "Started on" | grep -v "net.sf.picard.sam.MarkDuplicates INPUT"
	checkdiff -qr -x *log_jobs -x colinfo -x *_fastqc -x *bam.dupmetrics tmp/one_exome_yri_chr2122/samples/NA19240chr2122 expected/one_exome_yri_chr2122/samples/NA19240chr2122
} {}

test process {process_illumina exomes yri chr2122} {
	cd $::bigtestdir
	file delete -force tmp/exomes_yri_chr2122
	file mkdir tmp/exomes_yri_chr2122/samples
	foreach sample {
		NA19238chr2122  NA19239chr2122  NA19240chr2122
	} {
		file mkdir tmp/exomes_yri_chr2122/samples/$sample/fastq
		foreach file [glob ori/exomes_yri_chr2122.start/samples/$sample/fastq/*] {
			file copy $file tmp/exomes_yri_chr2122/samples/$sample/fastq/[file tail $file]
		}
	}
	if 0 {
		foreach sample {NA19238chr2122 NA19239chr2122 NA19240chr2122} {
			exec cp -a expected/exomes_yri_chr2122/samples/${sample}/map-rdsbwa-${sample}.bam tmp/exomes_yri_chr2122/samples/${sample}/
			exec cp -a expected/exomes_yri_chr2122/samples/${sample}/map-rdsbwa-${sample}.bam.bai tmp/exomes_yri_chr2122/samples/${sample}/
		}
		exec touch {*}[glob tmp/exomes_yri_chr2122/samples/*/map-*.bam*]
	}
	# cg process_illumina --stack 1 --verbose 2 -d 2 -split 1 -dbdir refseqtest/hg19 tests/yri_exome
	cg process_illumina --stack 1 --verbose 2 -split 1 -dbdir refseqtest/hg19 tmp/exomes_yri_chr2122 2>@ stderr >@ stdout
	# check vs expected
	checkdiff -y --suppress-common-lines tmp/exomes_yri_chr2122/samples/NA19238chr2122/map-dsbwa-NA19238chr2122.bam.dupmetrics expected/exomes_yri_chr2122/samples/NA19238chr2122/map-dsbwa-NA19238chr2122.bam.dupmetrics | grep -v "Started on"
	checkdiff -qr -x *log_jobs -x colinfo -x *_fastqc -x *bam.dupmetrics tmp/exomes_yri_chr2122 expected/exomes_yri_chr2122
	# could have used this, but previous is faster
	# cg tsvdiff -q 1 -x log_jobs -x colinfo -x _fastqc -x bam.dupmetrics tmp/exomes_yri_chr2122 expected/exomes_yri_chr2122
} {}

if 0 {
cg select -f 'chromosome begin end type ref alt zyg-gatk-* zyg-sam-*' tmp/exomes_yri_chr2122/compar/compar-exomes_yri_chr2122.tsv tmp/temp1
cg select -f 'chromosome begin end type ref alt zyg-gatk-* zyg-sam-*' expected/exomes_yri_chr2122/compar/compar-exomes_yri_chr2122.tsv tmp/temp2
kdiff3 tmp/temp1 tmp/temp2
tdiff tmp/temp1 tmp/temp2 | less
}

test process {process_sample genome yri chr2122} {
	cd $::bigtestdir
	set ref $::bigtestdir/refseqtest/hg19
	file delete -force tmp/one_genome_yri_chr2122
	file mkdir tmp/one_genome_yri_chr2122/samples/testNA19240chr2122cg
	mklink ori/genomes_yri_chr2122.start/samples/testNA19240chr2122cg.ori tmp/one_genome_yri_chr2122/samples/testNA19240chr2122cg/ori
	cg process_sample --stack 1 --verbose 2 -split 1 -dbdir $ref tmp/one_genome_yri_chr2122/samples/testNA19240chr2122cg 2>@ stderr >@ stdout
	# check vs expected
	checkdiff -y --suppress-common-lines tmp/one_genome_yri_chr2122/samples/testNA19240chr2122cg/summary-testNA19240chr2122cg.txt expected/genomes_yri_chr2122/samples/testNA19240chr2122cg/summary-testNA19240chr2122cg.txt | grep -v "finished.*finished"
	checkdiff -y --suppress-common-lines tmp/one_genome_yri_chr2122/samples/testNA19240chr2122cg/info_analysis.tsv expected/genomes_yri_chr2122/samples/testNA19240chr2122cg/info_analysis.tsv | grep -v "command"
	checkdiff -qr -x info_analysis.tsv -x log_jobs -x summary-testNA19240chr2122cg.txt tmp/one_genome_yri_chr2122/samples/testNA19240chr2122cg expected/genomes_yri_chr2122/samples/testNA19240chr2122cg
	# file_write tmp/temp $e
} {}

test process {genomes yri chr2122} {
	cd $::bigtestdir	
	file delete -force tmp/genomes_yri_chr2122
	file mkdir tmp/genomes_yri_chr2122
	foreach sample {
		testNA19238chr2122cg testNA19239chr2122cg testNA19240chr2122cg
		testNA19240chr21il
	} {
		file mkdir tmp/genomes_yri_chr2122/samples/$sample
		mklink ori/genomes_yri_chr2122.start/samples/$sample.ori tmp/genomes_yri_chr2122/samples/$sample/ori
	}
	# cg process_project --stack 1 --verbose 2 -d 2 -split 1 -dbdir /complgen/refseq/testdb2/hg19 tmp/genomes_yri_chr2122
	cg process_project --stack 1 --verbose 2 -split 1 -dbdir refseqtest/hg19 tmp/genomes_yri_chr2122 2>@ stderr >@ stdout
	# check vs expected
	foreach cgsample {testNA19238chr2122cg testNA19239chr2122cg testNA19240chr2122cg} {
		checkdiff -y --suppress-common-lines tmp/genomes_yri_chr2122/samples/$cgsample/summary-$cgsample.txt expected/genomes_yri_chr2122/samples/$cgsample/summary-$cgsample.txt | grep -v "finished.*finished"
	}
	checkdiff -y --suppress-common-lines tmp/genomes_yri_chr2122/samples/testNA19240chr21il/map-dsbwa-testNA19240chr21il.bam.dupmetrics expected/genomes_yri_chr2122/samples/testNA19240chr21il/map-dsbwa-testNA19240chr21il.bam.dupmetrics | grep -v "Started on"
	checkdiff -qr -x *log_jobs -x *_fastqc -x summary-* -x *dupmetrics -x colinfo tmp/genomes_yri_chr2122 expected/genomes_yri_chr2122
	# file_write temp $e
} {}

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
	file copy ori/wgs2.mastr/amplicons-wgs2.tsv tmp/wgs2.mastr
	# file copy ori/mastr_116068_116083/demultiplex_stats.tsv tmp/mastr_116068_116083
	# if you want to see output while running
	 cg process_mastr --stack 1 --verbose 2 -split 1 tmp/wgs2.mastr tmp/mastr_116068_116083 refseqtest/hg19 2>@ stderr >@ stdout
	# no output while running
	# cg process_mastr --stack 1 --verbose 2 -split 1 tmp/wgs2.mastr tmp/mastr_116068_116083 refseqtest/hg19
	# check vs expected
	checkdiff -qr -x *log_jobs -x *hsmetrics -x colinfo -x mastr_116068_116083.html tmp/mastr_116068_116083 expected/mastr_116068_116083
	checkdiff -y --suppress-common-lines tmp/mastr_116068_116083/mastr_116068_116083.html expected/mastr_116068_116083/mastr_116068_116083.html | grep -v -E {HistogramID|htmlwidget-|^<!|^<h2>20}
	foreach sample [dirglob tmp/mastr_116068_116083 ceph*] {
		checkdiff -y --suppress-common-lines tmp/mastr_116068_116083/$sample/crsbwa-$sample.hsmetrics expected/mastr_116068_116083/$sample/crsbwa-$sample.hsmetrics | grep -v -E "Started on|net.sf.picard.analysis.directed.CalculateHsMetrics BAIT_INT"
	}
	# could have used this, but previous is faster
	# cg tsvdiff -q 1 -x log_jobs -x mastr_116068_116083.html tmp/mastr_116068_116083 expected/mastr_116068_116083
} {}

cd $keepdir

testsummarize
