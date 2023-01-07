#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test flames {basic SIRV test} {
	cd $::testdir
	file delete -force tmp/sirv_flames
	file mkdir tmp/sirv_flames/fastq
	foreach file [glob -nocomplain data/SIRV-flames/fastq/*] {
		mklink $file tmp/sirv_flames/fastq/[file tail $file]
	}
	foreach file [glob -nocomplain data/SIRV-flames/*] {
		if {$file eq "data/SIRV-flames/fastq"} continue
		mklink $file tmp/sirv_flames/[file tail $file]
	}
	# flames works directly form fastq (and we could give the fastq directory instead of the bam file
	# going via bam to fit into the the typical workflow
	cg refseq_minimap2 tmp/sirv_flames/SIRV_isoforms_multi-fasta_170612a.fasta splice
	cg map \
		-method minimap2 -preset splice -paired 0 \
		tmp/sirv_flames/map-sirv_flames.bam \
		tmp/sirv_flames/SIRV_isoforms_multi-fasta_170612a.fasta \
		tmp/sirv_flames/sirv_flames \
		tmp/sirv_flames/fastq/sample1.fastq.gz tmp/sirv_flames/fastq/sample2.fastq.gz
	exec samtools index tmp/sirv_flames/map-sirv_flames.bam
	cg iso_flames \
		-refseq tmp/sirv_flames/SIRV_isoforms_multi-fasta_170612a.fasta \
		-reftranscripts tmp/sirv_flames/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf \
		tmp/sirv_flames/map-sirv_flames.bam
	# check vs expected
	exec diff tmp/sirv_flames/isoform_counts-flames-sirv_flames.tsv data/isoform_counts-flames-sirv_flames.tsv
	exec diff tmp/sirv_flames/gene_counts-flames-sirv_flames.tsv data/gene_counts-flames-sirv_flames.tsv
} {}

testsummarize

