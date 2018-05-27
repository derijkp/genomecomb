#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test map_bwa {map_bwa basic} {
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	cg map_bwa -stack 1 -paired 1 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[lsort -dict [glob tmp/*.fq.gz]]
	# chr21:42730799-42762826
	exec samtools view tmp/ali.bam > tmp/ali.sam
	exec samtools view data/bwa.bam > tmp/expected.sam
	exec diff tmp/ali.sam tmp/expected.sam
} {}

test map_bowtie2 {map_bowtie2 basic} {
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	cg map_bowtie2 -stack 1 -paired 1 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[lsort -dict [glob tmp/*.fq.gz]]
	exec samtools view tmp/ali.bam > tmp/ali.sam
	exec samtools view data/bowtie2.bam > tmp/expected.sam
} {}

#test map_minimap2 {map_minimap2 basic} {
#	test_cleantmp
#	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
#	cg map_minimap2 -stack 1 -preset sr tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[lsort -dict [glob tmp/*.fq.gz]]
#	# chr21:42730799-42762826
#	exec samtools view tmp/ali.bam > tmp/ali.sam
#	exec samtools view data/minimap2.bam > tmp/expected.sam
#	exec diff tmp/ali.sam tmp/expected.sam
#} {}

test map_minimap2 {map_minimap2 paired} {
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	cg map_minimap2 -stack 1 -paired 1 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[lsort -dict [glob tmp/*.fq.gz]]
	# chr21:42730799-42762826
	exec samtools view tmp/ali.bam > tmp/ali.sam
	exec samtools view data/minimap2-p.bam > tmp/expected.sam
	exec diff tmp/ali.sam tmp/expected.sam
} {}

test map_ngmlr {map_ngmlr basic} {
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	cg map_ngmlr -stack 1 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[lsort -dict [glob tmp/*.fq.gz]]
	# using samtools merge may result in differently (although still correctly) ordered bam each run
	# so first check if sorted correctly (only chromosome and start)
	cg sam2tsv tmp/ali.bam | cg select -f {qname chromosome=$refname begin e=$end strand mapquality ref2 begin2 strand2 tlen unmapped mateunmapped read secondary qcfail duplicate supplementary cigar seqlen seq quality} > tmp/ali.tsv
	cg checktsv tmp/ali.tsv
	# now check contents (compare sorted)
	cg sam2tsv tmp/ali.bam | cg select -s {refname begin end qname} -f {qname refname begin end strand mapquality ref2 begin2 strand2 tlen unmapped mateunmapped read secondary qcfail duplicate supplementary cigar seqlen seq quality} > tmp/ali.tsv
	cg sam2tsv data/ngmlr.bam | cg select -s {refname begin end qname} -f {qname refname begin end strand mapquality ref2 begin2 strand2 tlen unmapped mateunmapped read secondary qcfail duplicate supplementary cigar seqlen seq quality} > tmp/expected.tsv
	cg tsvdiff tmp/ali.tsv tmp/expected.tsv
} {}

test map_ngmlr {map_ngmlr 7 files -m 2} {
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	for {set i 3} {$i < 8} {incr i} {
		file copy data/seq_R1.fq.gz tmp/seq_R$i.fq.gz
	}
	cg map_ngmlr -stack 1 -maxopenfiles 2 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[lsort -dict [glob tmp/*.fq.gz]]
	# using samtools merge may result in differently (although still correctly) ordered bam each run
	# so first check if sorted correctly (only chromosome and start)
	cg sam2tsv tmp/ali.bam | cg select -f {qname chromosome=$refname begin e=$end strand mapquality ref2 begin2 strand2 tlen unmapped mateunmapped read secondary qcfail duplicate supplementary cigar seqlen seq quality} > tmp/ali.tsv
	cg checktsv tmp/ali.tsv
	# now check contents (compare sorted)
	cg sam2tsv tmp/ali.bam | cg select -s {refname begin end qname} -f {qname refname begin end strand mapquality ref2 begin2 strand2 tlen unmapped mateunmapped read secondary qcfail duplicate supplementary cigar seqlen seq quality} | uniq > tmp/ali.tsv
	cg sam2tsv data/ngmlr.bam | cg select -s {refname begin end qname} -f {qname refname begin end strand mapquality ref2 begin2 strand2 tlen unmapped mateunmapped read secondary qcfail duplicate supplementary cigar seqlen seq quality} > tmp/expected.tsv
	cg tsvdiff tmp/ali.tsv tmp/expected.tsv
	lindex [cg sam2tsv tmp/ali.bam | cg select -g all] end
} {700}

test map_ngmlr {map_ngmlr 4 files -m 2} {
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	for {set i 3} {$i < 5} {incr i} {
		file copy data/seq_R1.fq.gz tmp/seq_R$i.fq.gz
	}
	cg map_ngmlr -stack 1 -maxopenfiles 2 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[lsort -dict [glob tmp/*.fq.gz]]
	# using samtools merge may result in differently (although still correctly) ordered bam each run
	# so first check if sorted correctly (only chromosome and start)
	cg sam2tsv tmp/ali.bam | cg select -f {qname chromosome=$refname begin e=$end strand mapquality ref2 begin2 strand2 tlen unmapped mateunmapped read secondary qcfail duplicate supplementary cigar seqlen seq quality} > tmp/ali.tsv
	cg checktsv tmp/ali.tsv
	# now check contents (compare sorted)
	cg sam2tsv tmp/ali.bam | cg select -s {refname begin end qname} -f {qname refname begin end strand mapquality ref2 begin2 strand2 tlen unmapped mateunmapped read secondary qcfail duplicate supplementary cigar seqlen seq quality} | uniq > tmp/ali.tsv
	cg sam2tsv data/ngmlr.bam | cg select -s {refname begin end qname} -f {qname refname begin end strand mapquality ref2 begin2 strand2 tlen unmapped mateunmapped read secondary qcfail duplicate supplementary cigar seqlen seq quality} > tmp/expected.tsv
	cg tsvdiff tmp/ali.tsv tmp/expected.tsv
	lindex [cg sam2tsv tmp/ali.bam | cg select -g all] end
} {400}

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
