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

#test map_minimap2 {map_minimap2 basic} {
#	test_cleantmp
#	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
#	cg map_minimap2 -stack 1 -v 2 -preset sr tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[lsort -dict [glob tmp/*.fq.gz]] >@ stdout 2>@ stderr
#	# chr21:42730799-42762826
#	exec samtools view tmp/ali.bam > tmp/ali.sam
#	exec samtools view data/minimap2.bam > tmp/expected.sam
#	exec diff tmp/ali.sam tmp/expected.sam
#} {}

test map_bwa {map_minimap2 paired} {
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	cg map_minimap2 -stack 1 -v 2 -paired 1 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[lsort -dict [glob tmp/*.fq.gz]] >@ stdout 2>@ stderr
	# chr21:42730799-42762826
	exec samtools view tmp/ali.bam > tmp/ali.sam
	exec samtools view data/bwa.bam > tmp/expected.sam
	exec diff tmp/ali.sam tmp/expected.sam
} {}

test realign {realign_gatk basic} {
	file copy data/bwa.bam tmp/bwa.bam
	cg realign_gatk -stack 1 tmp/bwa.bam tmp/ratest.bam $::refseqdir/hg19
	exec samtools view tmp/ratest.bam > tmp/ratest.sam
	exec samtools view data/ratest-gatk.bam > tmp/ratest-gatk.sam
	exec diff tmp/ratest.sam tmp/ratest-gatk.sam
} {}

test realign {realign_abra basic} {
	cg realign_abra -stack 1 data/bwa.bam tmp/ratest.bam $::refseqdir/hg19
	cg sam2tsv tmp/ratest.bam tmp/ratest.tsv
	cg sam2tsv data/ratest-abra.bam tmp/expected.tsv
	exec diff tmp/ratest.tsv tmp/expected.tsv
} {}

testsummarize
