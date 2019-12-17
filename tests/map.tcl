#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test map_bwa {map_bwa basic} {
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	cg map_bwa -stack 1 -paired 1 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[lsort -dict [glob tmp/*.fq.gz]]
	# chr21:42730799-42762826
	exec samtools view -h tmp/ali.bam > tmp/ali.sam
	catch {exec diff -I {@PG	ID:bwa	PN:bwa} tmp/ali.sam data/bwa.sam}
} 0

test map_bwa {map_bwa multiple} {
	test_cleantmp
	set temp [split [string trim [exec zcat data/seq_R1.fq.gz]] \n]
	file_write tmp/seq1_R1.fq [join [lrange $temp 0 199] \n]\n
	file_write tmp/seq2_R1.fq [join [lrange $temp 200 end] \n]\n
	set temp [split [string trim [exec zcat data/seq_R2.fq.gz]] \n]
	file_write tmp/seq1_R2.fq [join [lrange $temp 0 199] \n]\n
	file_write tmp/seq2_R2.fq [join [lrange $temp 200 end] \n]\n
	cg map_bwa -stack 1 -paired 1 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[lsort -dict [glob tmp/*.fq]]
	# chr21:42730799-42762826
	exec samtools view -h tmp/ali.bam > tmp/ali.sam
	catch {exec diff -I {@PG	ID:bwa	PN:bwa} tmp/ali.sam data/bwa.sam}
} 0

test map_bwa {map_bwa cram} {
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	cg map_bwa -stack 1 -paired 1 tmp/ali.cram $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[lsort -dict [glob tmp/*.fq.gz]]
	# chr21:42730799-42762826
	dbdir $::refseqdir/hg19
	exec samtools view -h tmp/ali.cram > tmp/ali.sam
	cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA} tmp/ali.sam tmp/ali.sam.tsv
	cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA} data/bwa.sam tmp/bwa.sam.tsv
	catch {cg tsvdiff tmp/ali.sam.tsv tmp/bwa.sam.tsv}
} 0

test map_bowtie2 {map_bowtie2 basic} {
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	cg map_bowtie2 -stack 1 -paired 1 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[lsort -dict [glob tmp/*.fq.gz]]
	exec samtools view tmp/ali.bam | cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA} | cg select -rf {MC MQ} > tmp/ali.tsv
	exec samtools view data/bowtie2.bam | cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA} | cg select -rf {MC MQ} > tmp/expected.tsv
	catch {cg tsvdiff tmp/ali.tsv tmp/expected.tsv}
} 0

#test map_minimap2 {map_minimap2 basic} {
#	test_cleantmp
#	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
#	cg map_minimap2 -stack 1 -paired 0 -preset sr tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[lsort -dict [glob tmp/*.fq.gz]]
#	# chr21:42730799-42762826
#	exec samtools view tmp/ali.bam > tmp/ali.sam
#	exec samtools view data/minimap2.bam > tmp/expected.sam
#	exec diff tmp/ali.sam tmp/expected.sam
#} {}

test map_minimap2 {map_minimap2 paired} {
	if {![file exists $::refseqdir/hg19/genome_hg19.ifas.minimap2.sr]} {
		error "minimap2 sr index not made"
	}
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	cg map_minimap2 -stack 1 -paired 1 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[lsort -dict [glob tmp/*.fq.gz]]
	# chr21:42730799-42762826
	exec samtools view -h tmp/ali.bam > tmp/ali.sam
	cg sam2tsv -fields {RG NM AS nn tp cm s1 s2 MD MQ MC ms de rl} tmp/ali.sam tmp/ali.tsv
	cg select -rf {de rl} tmp/ali.tsv tmp/alis.tsv
	cg sam2tsv -fields {RG NM AS nn tp cm s1 s2 MD MQ MC ms} data/minimap2-p.sam tmp/expected.tsv
	catch {cg tsvdiff tmp/alis.tsv tmp/expected.tsv}
} 0

test map_minimap2 {error dir as refseq} {
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	file mkdir tmp/ref
	cg map_minimap2 -paired 1 tmp/ali.bam tmp/ref NA19240m {*}[lsort -dict [glob tmp/*.fq.gz]]
} {could not properly index */tmp/ref: contains no sequences} error match

test map_minimap2 {error missing fastq} {
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	mklink tmp/doesnotexists.fq tmp/bla.fq
	cg map_minimap2 -stack 1 -paired 0 -preset sr tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m tmp/bla.fq
} {*ERROR: failed to open file *tmp/bla.fq*} match error

test map_ngmlr {map_ngmlr basic} {
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	cg map_ngmlr -stack 1 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[lsort -dict [glob tmp/*.fq.gz]]
	# using samtools merge may result in differently (although still correctly) ordered bam each run
	# so first check if sorted correctly (only chromosome and start)
	cg sam2tsv tmp/ali.bam | cg select -f {qname chromosome begin e=$end strand mapquality ref2 begin2 strand2 tlen unmapped mateunmapped read secondary qcfail duplicate supplementary cigar seqlen seq quality} > tmp/ali.tsv
	cg checktsv tmp/ali.tsv
	# now check contents (compare sorted)
	cg sam2tsv tmp/ali.bam | cg select -s {chromosome begin end qname} -f {qname chromosome begin end strand mapquality ref2 begin2 strand2 tlen unmapped mateunmapped read secondary qcfail duplicate supplementary cigar seqlen seq quality} > tmp/ali.tsv
	cg sam2tsv data/ngmlr.bam | cg select -s {chromosome begin end qname} -f {qname chromosome begin end strand mapquality ref2 begin2 strand2 tlen unmapped mateunmapped read secondary qcfail duplicate supplementary cigar seqlen seq quality} > tmp/expected.tsv
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
	cg sam2tsv tmp/ali.bam | cg select -f {qname chromosome begin e=$end strand mapquality ref2 begin2 strand2 tlen unmapped mateunmapped read secondary qcfail duplicate supplementary cigar seqlen seq quality} > tmp/ali.tsv
	cg checktsv tmp/ali.tsv
	# now check contents (compare sorted)
	cg sam2tsv tmp/ali.bam | cg select -s {chromosome begin end qname} -f {qname chromosome begin end strand mapquality ref2 begin2 strand2 tlen unmapped mateunmapped read secondary qcfail duplicate supplementary cigar seqlen seq quality} | uniq > tmp/ali.tsv
	cg sam2tsv data/ngmlr.bam | cg select -s {chromosome begin end qname} -f {qname chromosome begin end strand mapquality ref2 begin2 strand2 tlen unmapped mateunmapped read secondary qcfail duplicate supplementary cigar seqlen seq quality} > tmp/expected.tsv
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
	cg sam2tsv tmp/ali.bam | cg select -f {qname chromosome begin e=$end strand mapquality ref2 begin2 strand2 tlen unmapped mateunmapped read secondary qcfail duplicate supplementary cigar seqlen seq quality} > tmp/ali.tsv
	cg checktsv tmp/ali.tsv
	# now check contents (compare sorted)
	cg sam2tsv tmp/ali.bam | cg select -s {chromosome begin end qname} -f {qname chromosome begin end strand mapquality ref2 begin2 strand2 tlen unmapped mateunmapped read secondary qcfail duplicate supplementary cigar seqlen seq quality} | uniq > tmp/ali.tsv
	cg sam2tsv data/ngmlr.bam | cg select -s {chromosome begin end qname} -f {qname chromosome begin end strand mapquality ref2 begin2 strand2 tlen unmapped mateunmapped read secondary qcfail duplicate supplementary cigar seqlen seq quality} > tmp/expected.tsv
	cg tsvdiff tmp/ali.tsv tmp/expected.tsv
	lindex [cg sam2tsv tmp/ali.bam | cg select -g all] end
} {400}

testsummarize
