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

test process {process_illumina exomes yri} {
	cd /data/testdata/	
	file delete -force tmp/exomes_yri
	exec cp -al ori/exomes_yri.start tmp/exomes_yri
	# cg process_illumina --stack 1 --verbose 2 -d 2 -split 1 -dbdir refseqtest/hg19 tests/yri_exome
	cg process_illumina --stack 1 --verbose 2 -split 1 -dbdir refseqtest/hg19 tmp/exomes_yri 2>@ stderr >@ stdout
	# check vs expected
	checkdiff -y --suppress-common-lines tmp/exomes_yri/samples/NA19238/map-dsbwa-NA19238.bam.dupmetrics expected/exomes_yri/samples/NA19238/map-dsbwa-NA19238.bam.dupmetrics | grep -v "Started on"
	checkdiff -qr -x *log_jobs -x *_fastqc -x *bam.dupmetrics tmp/exomes_yri expected/exomes_yri
} {}

test process {process_illumina exomes yri} {
	cd /data/testdata/	
	file delete -force tmp/exomes_yrit
	exec cp -al ori/exomes_yri.start tmp/exomes_yrit
	# cg process_illumina --stack 1 --verbose 2 -d 2 -split 1 -dbdir refseqtest/hg19 tests/yri_exome
	cg process_illumina --stack 1 --verbose 2 -split 1 -dbdir /data/genomecomb.testdata/refseqtest/hg19 tmp/exomes_yrit 2>@ stderr >@ stdout
	# check vs expected
	checkdiff -y --suppress-common-lines tmp/exomes_yrit/samples/NA19238/map-dsbwa-NA19238.bam.dupmetrics expected/exomes_yri/samples/NA19238/map-dsbwa-NA19238.bam.dupmetrics | grep -v "Started on"
	checkdiff -qr -x *log_jobs -x *_fastqc -x *bam.dupmetrics tmp/exomes_yrit expected/exomes_yri
} {}

if 0 {
cg select -f 'chromosome begin end type ref alt zyg-gatk-* zyg-sam-*' tmp/exomes_yri/compar/compar-exomes_yri.tsv tmp/temp1
cg select -f 'chromosome begin end type ref alt zyg-gatk-* zyg-sam-*' expected/exomes_yri/compar/compar-exomes_yri.tsv tmp/temp2
kdiff3 tmp/temp1 tmp/temp2
tdiff tmp/temp1 tmp/temp2 | less
}

test process {process_sample genome yri} {
	cd $::bigtestdir
	set ref $::bigtestdir/refseqtest/hg19
	set dest tmp/genomes_yri_one
	file delete -force $dest
	file mkdir $dest/samples/testNA19240cg
	mklink ori/genomes_yritrio.start/samples/testNA19240cg.ori $dest/samples/testNA19240cg/ori
	cg process_sample --stack 1 --verbose 2 -split 1 -dbdir $ref $dest/samples/testNA19240cg 2>@ stderr >@ stdout
	# check vs expected
	checkdiff -y --suppress-common-lines tmp/genomes_yri_one/samples/testNA19240cg/summary-testNA19240cg.txt expected/genomes_yri/samples/testNA19240cg/summary-testNA19240cg.txt | grep -v "finished.*finished"
	checkdiff -qr -x log_jobs -x summary-testNA19240cg.txt tmp/genomes_yri_one/samples/testNA19240cg expected/genomes_yri/samples/testNA19240cg
	# file_write tmp/temp $e
} {}

test process {genomes yri } {
	cd $::bigtestdir	
	set dest tmp/genomes_yri
	file delete -force $dest
	file mkdir $dest
	foreach sample {
		testNA19238cg testNA19239cg testNA19240cg
		testNA19240chr21il
	} {
		file mkdir $dest/samples/$sample
		mklink ori/genomes_yritrio.start/samples/$sample.ori $dest/samples/$sample/ori
	}
#	file mkdir $dest/samples/testNA19240chr21il/fastq
#	foreach file [glob ori/genomes_yritrio.start/samples/testNA19240chr21il.ori/*.fq*] {
#		mklink $file $dest/samples/testNA19240chr21il/fastq/[file tail $file]
#	}
	# mklink ori/genomes_yritrio.start/samples/testNA19240chr21il.ori/NA19240_GAIIx_100_chr21.bam $dest/samples/testNA19240chr21il/map-rdsbwa-testNA19240chr21il.bam
	# cg process_project --stack 1 --verbose 2 -d 2 -split 1 -dbdir /complgen/refseq/testdb2/hg19 tmp/genomes_yri
	cg process_project --stack 1 --verbose 2 -split 1 -dbdir refseqtest/hg19 $dest 2>@ stderr >@ stdout
	# check vs expected
	foreach cgsample {testNA19238cg testNA19239cg testNA19240cg} {
		checkdiff -y --suppress-common-lines tmp/genomes_yri/samples/$cgsample/summary-$cgsample.txt expected/genomes_yri/samples/$cgsample/summary-$cgsample.txt | grep -v "finished.*finished"
	}
	checkdiff -y --suppress-common-lines tmp/genomes_yri/samples/testNA19240chr21il/map-dsbwa-testNA19240chr21il.bam.dupmetrics expected/genomes_yri/samples/testNA19240chr21il/map-dsbwa-testNA19240chr21il.bam.dupmetrics | grep -v "Started on"
	checkdiff -qr -x *log_jobs -x *_fastqc -x summary-* -x *dupmetrics -x colinfo tmp/genomes_yri expected/genomes_yri
	# file_write tmp/temp $e
} {}

cd $keepdir

testsummarize
