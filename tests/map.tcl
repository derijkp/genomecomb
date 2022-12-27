#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test map_bwa {map_bwa basic} {
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	cg map_bwa -stack 1 -paired 1 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[bsort [glob tmp/*.fq.gz]]
	# chr21:42730799-42762826
	exec samtools view --no-PG -h tmp/ali.bam > tmp/ali.sam
	cg sam_sort -sort name data/bwa.sam tmp/expected.sam
	catch {
		exec diff -I {@PG	} -I {@PG	}  -I {@HD	} tmp/ali.sam tmp/expected.sam
	}
} 0

test map_bwa {map_bwa -paired 0} {
	test_cleantmp
	file copy data/seq_R1.fq.gz tmp
	cg map_bwa -stack 1 -paired 0 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[bsort [glob tmp/*.fq.gz]]
	set expectdfile data/bwa.sam
	set otherfields {AS XS MQ MC ms MD RG NM XA YS YT}
	set removefields {MC MQ YS YT XS read ms mapquality mateunmapped ref2	begin2	strand2	tlen	pair	properpair}
	exec samtools view --no-PG tmp/ali.bam | cg sam2tsv -fields $otherfields \
		| cg select -s {chromosome begin end} -rf $removefields > tmp/ali.tsv
	exec samtools view --no-PG $expectdfile | cg sam2tsv -fields $otherfields \
		| cg select -s {chromosome begin end} -q {$read == 1} -rf $removefields > tmp/expected.tsv
	cg tsvdiff tmp/ali.tsv tmp/expected.tsv
} {}

test map_bwa {map_bwa to stdout} {
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	cg map_bwa -stack 1 -paired 1 -compressionlevel 1 -.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[bsort [glob tmp/*.fq.gz]] > tmp/ali.bam
	# chr21:42730799-42762826
	exec samtools view --no-PG -h tmp/ali.bam > tmp/ali.sam
	cg sam_sort -sort name data/bwa.sam tmp/expected.sam
	catch {exec diff -I {@PG	} -I {@HD	} tmp/ali.sam tmp/expected.sam}
} 0

test map_bwa {map_bwa refseq error} {
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	mklink $::refseqdir/hg19/genome_hg19.ifas tmp/genome_hg19.ifas
	cg map_bwa -stack 1 -paired 1 tmp/ali.bam tmp/genome_hg19.ifas NA19240m {*}[bsort [glob tmp/*.fq.gz]]
} {The bwa version of the refseq does not exist (should be at */tmp/genome_hg19.ifas.bwa/genome_hg19.fa)
You can create it using:
cg refseq_bwa *genome_hg19.ifas*} match error

# This one takes to long to run every time
#test map_bwa {map_bwa cg refseq_bwa} {
#	test_cleantmp
#	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
#	mklink $::refseqdir/hg19/genome_hg19.ifas tmp/genome_hg19.ifas
#	cg refseq_bwa tmp/genome_hg19.ifas
#	cg map_bwa -stack 1 -paired 1 tmp/ali.bam tmp/genome_hg19.ifas NA19240m {*}[bsort [glob tmp/*.fq.gz]]
#	exec samtools view --no-PG -h tmp/ali.bam > tmp/ali.sam
#	cg sam_sort -sort name data/bwa.sam tmp/expected.sam
#	file delete -force [glob tmp/genome_hg19.ifas*]
#	catch {exec diff -I {@PG	ID:bwa	PN:bwa} -I {@HD	} tmp/ali.sam tmp/expected.sam}
#} 0

test map_bwa {map_bwa cram} {
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	cg map_bwa -stack 1 -paired 1 tmp/ali.cram $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[bsort [glob tmp/*.fq.gz]]
	# chr21:42730799-42762826
	dbdir $::refseqdir/hg19
	exec samtools sort --no-PG -O sam tmp/ali.cram > tmp/ali.sam
	cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA} tmp/ali.sam tmp/ali.sam.tsv
	cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA} data/bwa.sam tmp/bwa.sam.tsv
	catch {cg tsvdiff tmp/ali.sam.tsv tmp/bwa.sam.tsv}
} 0

test map {map -method bwa basic} {
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	cg map -stack 1 -method bwa -sort coordinate -paired 1 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[bsort [glob tmp/*.fq.gz]]
	exec samtools view --no-PG -h tmp/ali.bam > tmp/ali.sam
	cg sam_sort -sort coordinate data/bwa.sam tmp/expected.sam
	catch {exec diff -I {@PG	} -I {@HD	} tmp/ali.sam tmp/expected.sam}
} 0

test map {map -method bwa -paired 0} {
	test_cleantmp
	file copy data/seq_R1.fq.gz tmp
	cg map -method bwa -stack 1 -paired 0 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[bsort [glob tmp/*.fq.gz]]
	exec samtools view --no-PG tmp/ali.bam | cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA YS YT} \
		| cg select -s {chromosome begin end} \
		-rf {MC MQ YS YT XS read ms mapquality mateunmapped ref2	begin2	strand2	tlen	pair	properpair} > tmp/ali.tsv
	exec samtools view --no-PG data/bwa.sam | cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA YS YT} \
		| cg select -s {chromosome begin end} -q {$read == 1} \
		-rf {MC MQ YS YT XS read ms mapquality mateunmapped ref2	begin2	strand2	tlen	pair	properpair} > tmp/expected.tsv
	catch {cg tsvdiff tmp/ali.tsv tmp/expected.tsv}
} 0

test map {map -method bwa multiple} {
	test_cleantmp
	set temp [split [string trim [exec zcat data/seq_R1.fq.gz]] \n]
	file_write tmp/seq1_R1.fq [join [lrange $temp 0 199] \n]\n
	file_write tmp/seq2_R1.fq [join [lrange $temp 200 end] \n]\n
	set temp [split [string trim [exec zcat data/seq_R2.fq.gz]] \n]
	file_write tmp/seq1_R2.fq [join [lrange $temp 0 199] \n]\n
	file_write tmp/seq2_R2.fq [join [lrange $temp 200 end] \n]\n
	cg map -method bwa -stack 1 -paired 1 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[bsort [glob tmp/*.fq]]
	# chr21:42730799-42762826
	exec samtools view --no-PG -h tmp/ali.bam > tmp/ali.sam
	catch {exec diff -I {@PG	} tmp/ali.sam data/bwa.sam}
} 0

test map {map_bwa paired -.sam} {
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	cg map_bwa -paired 1 -.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[bsort [glob tmp/*.fq.gz]] > tmp/ali.bam
	cg sam_sort -sort coordinate tmp/ali.bam tmp/ali.sam
	cg sam_sort -sort coordinate data/bwa.sam tmp/expected.sam
	catch {exec diff -I {@PG	} -I {@HD	} tmp/ali.sam tmp/expected.sam}
} 0

test map_bowtie2 {map_bowtie2 basic} {
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	cg map_bowtie2 -stack 1 -paired 1 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[bsort [glob tmp/*.fq.gz]]
	exec samtools view --no-PG tmp/ali.bam | cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA} | cg select -s {chromosome begin end} -rf {MC MQ} > tmp/ali.tsv
	exec samtools view --no-PG data/bowtie2.sam | cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA} | cg select -s {chromosome begin end} -rf {MC MQ} > tmp/expected.tsv
	catch {cg tsvdiff tmp/ali.tsv tmp/expected.tsv}
} 0

test map_bowtie2 {map_bowtie2 -paired 0} {
	test_cleantmp
	file copy data/seq_R1.fq.gz tmp
	cg map_bowtie2 -stack 1 -paired 0 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[bsort [glob tmp/*.fq.gz]]
	set expectdfile data/bowtie2.sam
	set otherfields {AS XS MC MQ YS YT RG NM XA s1 s2 cm de rl ms}
	set removefields {AS XS MC MQ YS YT RG NM XA s1 s2 cm de rl read ms mapquality mateunmapped ref2 begin2 strand2 tlen pair properpair}
	exec samtools view --no-PG tmp/ali.bam | cg sam2tsv -fields $otherfields \
		| cg select -f {chromosome	begin	end	strand {qname="[string range $qname 0 end-2]"} *} \
		| cg select -s {chromosome begin end} -rf $removefields > tmp/ali.tsv
	exec samtools view --no-PG $expectdfile | cg sam2tsv -fields $otherfields \
		| cg select -s {chromosome begin end} -q {$read == 1} -rf $removefields > tmp/expected.tsv
	catch {cg tsvdiff tmp/ali.tsv tmp/expected.tsv}
} 0

#test map_minimap2 {map_minimap2 basic} {
#	test_cleantmp
#	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
#	cg map_minimap2 -stack 1 -paired 0 -preset sr tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[bsort [glob tmp/*.fq.gz]]
#	# chr21:42730799-42762826
#	exec samtools view --no-PG tmp/ali.bam > tmp/ali.sam
#	exec samtools view --no-PG data/minimap2.bam > tmp/expected.sam
#	exec diff tmp/ali.sam tmp/expected.sam
#} {}

test map_minimap2 {map_minimap2 paired} {
	if {![file exists $::refseqdir/hg19/genome_hg19.ifas.minimap2.sr]} {
		error "minimap2 sr index not made"
	}
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	cg map_minimap2 -stack 1 -paired 1 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[bsort [glob tmp/*.fq.gz]]
	# chr21:42730799-42762826
	exec samtools sort tmp/ali.bam | samtools view --no-PG -h > tmp/ali.sam
	cg sam2tsv -fields {RG NM AS nn tp cm s1 s2 MD MQ MC ms} tmp/ali.sam > tmp/alis.tsv
	cg sam2tsv -fields {RG NM AS nn tp cm s1 s2 MD MQ MC ms} data/minimap2-p.sam tmp/expected.tsv
	catch {cg tsvdiff tmp/alis.tsv tmp/expected.tsv}
} 0

test map_minimap2 {map_minimap2 paired -.sam} {
	if {![file exists $::refseqdir/hg19/genome_hg19.ifas.minimap2.sr]} {
		error "minimap2 sr index not made"
	}
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	cd tmp
	cg map_minimap2 -stack 1 -paired 1 -.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[bsort [glob *.fq.gz]] > ali.bam
	cd ..
	if {![catch {glob tmp/-.*} result]} {
		error "$result file was made"
	}
	# chr21:42730799-42762826
	exec samtools sort tmp/ali.bam | samtools view --no-PG -h > tmp/ali.sam
	cg sam2tsv -fields {RG NM AS nn tp cm s1 s2 MD MQ MC ms} tmp/ali.sam > tmp/alis.tsv
	cg sam2tsv -fields {RG NM AS nn tp cm s1 s2 MD MQ MC ms} data/minimap2-p.sam tmp/expected.tsv
	catch {cg tsvdiff tmp/alis.tsv tmp/expected.tsv}
} 0

cd $::testdir

test map_minimap2 {map_minimap2 -paired 0} {
	test_cleantmp
	file copy data/seq_R1.fq.gz tmp
	cg map_minimap2 -stack 1 -paired 0 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[bsort [glob tmp/*.fq.gz]]
	set otherfields {AS XS MC MQ YS YT s1 s2 cm de rl ms}
	set removefields {AS XS MC MQ YS YT s1 s2 cm de rl read ms mapquality mateunmapped ref2 begin2 strand2 tlen pair properpair}
	exec samtools view --no-PG tmp/ali.bam | cg sam2tsv -fields $otherfields \
		| cg select -f {chromosome	begin	end	strand {qname="[string range $qname 0 end-2]"} *} \
		| cg select -s {chromosome begin end} -rf $removefields > tmp/ali.tsv
	exec samtools view --no-PG data/minimap2-p.sam | cg sam2tsv -fields $otherfields \
		| cg select -s {chromosome begin end} -q {$read == 1} -rf $removefields > tmp/expected.tsv
	# we have one slight difference in alignment for unpaired alignement!
	catch {cg tsvdiff tmp/ali.tsv tmp/expected.tsv >& tmp/diff}
	file_write tmp/expectdiff [deindent {
		diff tmp/ali.tsv tmp/expected.tsv
		header
		  chromosome	begin	end	strand	qname	qstart	qend	unmapped	secondary	qcfail	duplicate	supplementary	cigar	seqlen	seq	quality	other
		25c25
		< chr21	42752084	42752180	-	SRR792091.1603286	0	95	0	0	0	0	0	91M1D4M5S	100	CCACGTCATTCTGAGGTTCGGATCTGGCAGCCGCTCCTCTCACTTCCTCGGTTCCTTCTCCTCTTCCTCAAGTCACCCCCACAGTGACCACCAGCACCAC	@<??@::CCCCC<(B@@@DCDDBACA<B@BBBAD@DCDDDCADBBDDFFHHA2HGGEDAHGGGIIIIIIIIHF@)IIIIHFIIIGHGHFFHHDBD:F@?@	RG:Z:NA19240m NM:i:1 nn:i:0 tp:A:P MD:Z:91^T4
		---
		> chr21	42752084	42752175	-	SRR792091.1603286	0	91	0	0	0	0	0	91M9S	100	CCACGTCATTCTGAGGTTCGGATCTGGCAGCCGCTCCTCTCACTTCCTCGGTTCCTTCTCCTCTTCCTCAAGTCACCCCCACAGTGACCACCAGCACCAC	@<??@::CCCCC<(B@@@DCDDBACA<B@BBBAD@DCDDDCADBBDDFFHHA2HGGEDAHGGGIIIIIIIIHF@)IIIIHFIIIGHGHFFHHDBD:F@?@	RG:Z:NA19240m NM:i:0 nn:i:0 tp:A:P MD:Z:91
	}]\n
	exec diff tmp/diff tmp/expectdiff
} {}

test map_minimap2 {error dir as refseq} {
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	file mkdir tmp/ref
	mklink $::refseqdir/hg19/genome_hg19.ifas tmp/ref/genome_hg19.ifas
	cg map_minimap2 -paired 1 tmp/ali.bam tmp/ref NA19240m {*}[bsort [glob tmp/*.fq.gz]]
} {The minimap2 version for preset sr of the refseq does not exist (should be at *tmp/ref/*.minimap2.sr)
You can create it using:
cg refseq_minimap2 '*tmp/ref*' sr} error match

test map_minimap2 {error missing fastq} {
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	mklink tmp/doesnotexists.fq tmp/bla.fq
	cg map_minimap2 -stack 1 -paired 0 -preset sr tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m tmp/bla.fq
} {*ERROR: failed to open file *tmp/bla.fq*} match error

test map_ngmlr {map_ngmlr basic} {
	test_cleantmp
	cg zcat data/seq_R1.fq.gz data/seq_R2.fq.gz | cg bgzip > tmp/seq.fq.gz
	cg map_ngmlr -stack 1 -paired 0 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m tmp/seq.fq.gz
	# using samtools merge may result in differently (although still correctly) ordered bam each run
	# so first check if sorted correctly (only chromosome and start)
	cg sam2tsv tmp/ali.bam | cg select -f {qname chromosome begin e=$end strand mapquality ref2 begin2 strand2 tlen unmapped mateunmapped read secondary qcfail duplicate supplementary cigar seqlen seq quality} > tmp/ali.tsv
	# now check contents (compare sorted)
	cg sam2tsv tmp/ali.bam | cg select -s {chromosome begin end qname} -f {qname chromosome begin end strand mapquality ref2 begin2 strand2 tlen unmapped mateunmapped read secondary qcfail duplicate supplementary cigar seqlen seq quality} > tmp/ali.tsv
	cg sam2tsv data/ngmlr.bam | cg select -s {chromosome begin end qname} -f {qname chromosome begin end strand mapquality ref2 begin2 strand2 tlen unmapped mateunmapped read secondary qcfail duplicate supplementary cigar seqlen seq quality} > tmp/expected.tsv
	cg tsvdiff tmp/ali.tsv tmp/expected.tsv
} {}

test map_ngmlr {map -method ngmlr 7 files -m 2} {
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	for {set i 3} {$i < 8} {incr i} {
		file copy data/seq_R1.fq.gz tmp/seq_R$i.fq.gz
	}
	cg map -stack 1 -method ngmlr -paired 0 -maxopenfiles 2 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[bsort [glob tmp/*.fq.gz]]
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
	cg map -method ngmlr -stack 1 -paired 0 -maxopenfiles 2 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[bsort [glob tmp/*.fq.gz]]
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

test map_hisat2 {map_hisat2 paired} {
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	cg map_hisat2 -stack 1 -paired 1 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[bsort [glob tmp/*.fq.gz]]
	# chr21:42730799-42762826
	exec samtools view --no-PG -h tmp/ali.bam > tmp/ali.sam
	cg sam2tsv -fields {RG NM AS MD MQ MC XN XM XO XG YS YT NH ms} tmp/ali.sam \
		| cg select -s - -rf {other de rl ROW RG NM AS MD MQ MC XN XM XO XG YS YT NH ms} > tmp/alis.tsv
	cg sam2tsv -fields {RG NM AS MD MQ MC XN XM XO XG YS YT NH ms} data/hisat2-p.sam \
		| cg select -s - -rf {other de rl ROW RG NM AS MD MQ MC XN XM XO XG YS YT NH ms} > tmp/expected.tsv
	catch {cg tsvdiff tmp/alis.tsv tmp/expected.tsv}
} 0

if 0 {
# do not run test by default, star uses too much memory for most test machines

test map_star {map_star paired} {
	if {![file exists $::refseqdir/hg19/genome_hg19.ifas.star]} {
		error "star index not made"
	}
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	cg map_star -paired 1 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[bsort [glob tmp/*.fq.gz]]
	# chr21:42730799-42762826
	exec samtools view --no-PG -h tmp/ali.bam > tmp/ali.sam
	cg sam2tsv -fields {RG NM AS MD MQ MC XN XM XO XG YS YT NH ms} tmp/ali.sam  > tmp/alis.tsv
	cg sam2tsv -fields {RG NM AS MD MQ MC XN XM XO XG YS YT NH ms} data/star-p.sam  > tmp/expected.tsv
	catch {cg tsvdiff tmp/alis.tsv tmp/expected.tsv}
} 0

test map_star {map_star paired 2p} {
	if {![file exists $::refseqdir/hg19/genome_hg19.ifas.star]} {
		error "star index not made"
	}
	test_cleantmp
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp
	cg map_star -paired 1 -preset 2p tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[bsort [glob tmp/*.fq.gz]]
	# chr21:42730799-42762826
	exec samtools view --no-PG -h tmp/ali.bam > tmp/ali.sam
	cg sam2tsv -fields {RG NM AS MD MQ MC XN XM XO XG YS YT NH ms} tmp/ali.sam  > tmp/alis.tsv
	cg sam2tsv -fields {RG NM AS MD MQ MC XN XM XO XG YS YT NH ms} data/star-p.sam  > tmp/expected.tsv
	catch {cg tsvdiff tmp/alis.tsv tmp/expected.tsv}
} 0

test map_star {map_star -paired 0} {
	if {![file exists $::refseqdir/hg19/genome_hg19.ifas.star]} {
		error "star index not made"
	}
	test_cleantmp
	file copy -force data/seq_R1.fq.gz tmp
	cg map_star -paired 0 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[bsort [glob tmp/*.fq.gz]]
	# chr21:42730799-42762826
	exec samtools view --no-PG -h tmp/ali.bam > tmp/ali.sam
	cg sam2tsv -fields {RG NM AS MD MQ MC XN XM XO XG YS YT NH ms} tmp/ali.sam  > tmp/alis.tsv
	cg sam2tsv -fields {RG NM AS MD MQ MC XN XM XO XG YS YT NH ms} data/star.sam  > tmp/expected.tsv
	catch {cg tsvdiff tmp/alis.tsv tmp/expected.tsv}
} 0

test map {map -method star basic} {
	if {![file exists $::refseqdir/hg19/genome_hg19.ifas.star]} {
		error "star index not made"
	}
	test_cleantmp
	file copy -force data/seq_R1.fq.gz tmp
	cg map -stack 1 -method star -paired 0 tmp/ali.bam $::refseqdir/hg19/genome_hg19.ifas NA19240m {*}[bsort [glob tmp/*.fq.gz]]
	# chr21:42730799-42762826
	exec samtools view --no-PG -h tmp/ali.bam > tmp/ali.sam
	cg sam2tsv -fields {RG NM AS MD MQ MC XN XM XO XG YS YT NH ms} tmp/ali.sam  > tmp/alis.tsv
	cg sam2tsv -fields {RG NM AS MD MQ MC XN XM XO XG YS YT NH ms} data/star.sam  > tmp/expected.tsv
	catch {cg tsvdiff tmp/alis.tsv tmp/expected.tsv}
} 0

}

testsummarize
