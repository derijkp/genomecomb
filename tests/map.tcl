#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test map_bwa {map_bwa basic} {
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	cg map_bwa -stack 1 -v 2 -paired 1 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[lsort -dict [glob tmp/*.fq.gz]] >@ stdout 2>@ stderr
	# chr21:42730799-42762826
	exec samtools view tmp/ali.bam > tmp/ali.sam
	exec samtools view data/bwa.bam > tmp/expected.sam
	exec diff tmp/ali.sam tmp/expected.sam
} {}

test map_bowtie2 {map_bowtie2 basic} {
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	cg map_bowtie2 -stack 1 -v 2 -paired 1 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[lsort -dict [glob tmp/*.fq.gz]] >@ stdout 2>@ stderr
	exec samtools view tmp/ali.bam > tmp/ali.sam
	exec samtools view data/bowtie2.bam > tmp/expected.sam
} {}

testsummarize
