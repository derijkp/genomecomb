#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

proc splitfastqs {fastq1 pre1 post1 {perfile 20}} {
	set num 1
	catch {close $o1 ; close $f1}
	set f1 [gzopen $fastq1]
	set o1 [open $pre1${num}$post1 w]
	set lines 0
	while {[gets $f1 line1] != -1} {
		puts $o1 $line1
		incr lines
		if {$lines == $perfile} {
			close $o1
			if {[eof $f1]} break
			incr num
			set o1 [open $pre1${num}$post1 w]
			set lines 0
		}
	}
	close $o1; close $f1
	if {$lines == 0} {file delete $pre1${num}$post1}
}

test process_sample {bwa distrreg cram} {
	test_cleantmp
	file mkdir tmp/NA19240m/fastq
	file copy -force data/seq_R1.fq.gz data/seq_R2.fq.gz tmp/NA19240m/fastq
	cg process_sample -threads 1 -aligners bwa -aliformat cram -distrreg chr -dbdir $::refseqdir/hg19/genome_hg19.ifas tmp/NA19240m
	# chr21:42730799-42762826
	cg tsvdiff -x fastq -x info_analysis.tsv -x *.submitting -x fastqc* -x *.index -x log_jobs tmp/NA19240m data/NA19240m
} {}

test process_sample {bwa distrreg -removeduplicates 0 cram} {
	test_cleantmp
	file mkdir tmp/NA19240m/fastq
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp/NA19240m/fastq
	cg process_sample -aligners bwa -varcallers {} -removeduplicates 0 -aliformat cram -distrreg chr -dbdir $::refseqdir/hg19/genome_hg19.ifas tmp/NA19240m
	# chr21:42730799-42762826
	dbdir $::refseqdir/hg19
	exec samtools sort -O sam tmp/NA19240m/map-rsbwa-NA19240m.cram > tmp/ali.sam
	exec cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA} tmp/ali.sam | cg select -s {chromosome begin end qname} > tmp/ali.sam.tsv
	exec cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA} data/bwa.sam | cg select -s {chromosome begin end qname} > tmp/bwa.sam.tsv
	catch {cg tsvdiff tmp/ali.sam.tsv tmp/bwa.sam.tsv}
} 0

test process_sample {map bwa distrreg mutiple fastq} {
	test_cleantmp
	file mkdir tmp/NA19240m/fastq
	# make split up fastq
	splitfastqs data/seq_R1.fq.gz tmp/NA19240m/fastq/seq _R1.fq 20
	splitfastqs data/seq_R2.fq.gz tmp/NA19240m/fastq/seq _R2.fq 20
	file delete tmp/NA19240m/map-rdsbwa-NA19240m.bam
	cg process_sample -clip 0 -aligners bwa -varcallers {} -distrreg chr -varcallers {} -dbdir $::refseqdir/hg19/genome_hg19.ifas tmp/NA19240m >@ stdout 2>@ stderr
	# chr21:42730799-42762826
	exec samtools sort -O sam tmp/NA19240m/map-rdsbwa-NA19240m.bam > tmp/ali.sam
	exec cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA} tmp/ali.sam | cg select -rf {duplicate other properpair mapquality XS MQ} -s {chromosome begin end qname} > tmp/ali.sam.tsv
	exec cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA} data/bwa.sam | cg select -rf {duplicate other properpair mapquality XS MQ} -s {chromosome begin end qname} > tmp/bwa.sam.tsv
	catch {cg tsvdiff tmp/ali.sam.tsv tmp/bwa.sam.tsv}
} 0

test process_sample {map bwa distrreg mutiple fastq -maxfastqdistr 2} {
	test_cleantmp
	file mkdir tmp/NA19240m/fastq
	# make split up fastq
	splitfastqs data/seq_R1.fq.gz tmp/NA19240m/fastq/seq _R1.fq 20
	splitfastqs data/seq_R2.fq.gz tmp/NA19240m/fastq/seq _R2.fq 20
	file delete tmp/NA19240m/map-rdsbwa-NA19240m.bam
	cg process_sample -stack 1 -v 2 -clip 0 -maxfastqdistr 2 -aligners bwa -varcallers {} -distrreg chr -varcallers {} -dbdir $::refseqdir/hg19/genome_hg19.ifas tmp/NA19240m >@ stdout 2>@ stderr
	# chr21:42730799-42762826
	exec samtools sort -O sam tmp/NA19240m/map-rdsbwa-NA19240m.bam > tmp/ali.sam
	exec cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA} tmp/ali.sam | cg select -rf {duplicate other properpair mapquality XS MQ} -s {chromosome begin end qname} > tmp/ali.sam.tsv
	exec cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA} data/bwa.sam | cg select -rf {duplicate other properpair mapquality XS MQ} -s {chromosome begin end qname} > tmp/bwa.sam.tsv
	catch {cg tsvdiff tmp/ali.sam.tsv tmp/bwa.sam.tsv}
} 0

testsummarize

