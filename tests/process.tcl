#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl
set keepdir [pwd]
cd /complgen/testdata

if 0 {
	# make testdata
	file mkdir /complgen/testdata/yri_exome/samples/NA19238/fastq
	cd /complgen/testdata/yri_exome/samples/NA19238/fastq
	exec wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19238/sequence_read/SRR071173_2.filt.fastq.gz 2>@ stderr
	exec wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19238/sequence_read/SRR071173_1.filt.fastq.gz 2>@ stderr
	file mkdir /complgen/testdata/yri_exome/samples/NA19239/fastq
	cd /complgen/testdata/yri_exome/samples/NA19239/fastq
	exec wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19239/sequence_read/SRR792097_1.filt.fastq.gz 2>@ stderr
	exec wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19239/sequence_read/SRR792097_2.filt.fastq.gz 2>@ stderr
	file mkdir /complgen/testdata/yri_exome/samples/NA19240/fastq
	cd /complgen/testdata/yri_exome/samples/NA19240/fastq
	exec wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19240/sequence_read/SRR792091_1.filt.fastq.gz 2>@ stderr
	exec wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19240/sequence_read/SRR792091_2.filt.fastq.gz 2>@ stderr
}

test process {process_illumina} {
	file delete -force yri_exome_test
	exec cp -a yri_exome yri_exome_test
	cg process_illumina -split 1 -dbdir /complgen/refseq/hg19 yri_exome_test	
} {}

cd $keepdir

testsummarize
