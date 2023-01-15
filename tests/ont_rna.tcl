#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test ont_rna {flames basic SIRV test} {
	file delete -force tmp/sirv
	file mkdir tmp/sirv/fastq
	foreach file [glob -nocomplain data/SIRV-flames/fastq/*] {
		mklink $file tmp/sirv/fastq/[file tail $file]
	}
	foreach file [glob -nocomplain data/SIRV-flames/*] {
		if {$file eq "data/SIRV-flames/fastq"} continue
		mklink $file tmp/sirv/[file tail $file]
	}
	# flames works directly form fastq (and we could give the fastq directory instead of the bam file
	# going via bam to fit into the the typical workflow
	cg refseq_minimap2 tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta splice
	cg map \
		-method minimap2 -preset splice -paired 0 \
		tmp/sirv/map-sirv.bam \
		tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		tmp/sirv/sirv \
		tmp/sirv/fastq/sample1.fastq.gz tmp/sirv/fastq/sample2.fastq.gz
	exec samtools index tmp/sirv/map-sirv.bam
	cg iso_flames \
		-refseq tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		-reftranscripts tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf \
		tmp/sirv/map-sirv.bam
	# check vs expected
	exec diff tmp/sirv/isoform_counts-flames-sirv.tsv data/isoform_counts-flames-sirv.tsv
	exec diff tmp/sirv/gene_counts-flames-sirv.tsv data/gene_counts-flames-sirv.tsv
} {}

test ont_rna {flair basic SIRV test} {
	file delete -force tmp/sirv_flair
	file mkdir tmp/sirv/fastq
	foreach file [glob -nocomplain data/SIRV-flames/fastq/*] {
		mklink $file tmp/sirv/fastq/[file tail $file]
	}
	foreach file [glob -nocomplain data/SIRV-flames/*] {
		if {$file eq "data/SIRV-flames/fastq"} continue
		mklink $file tmp/sirv/[file tail $file]
	}
	cg refseq_minimap2 tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta splice
	cg map \
		-method minimap2 -preset splice -paired 0 \
		tmp/sirv/map-sirv.bam \
		tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		tmp/sirv/sirv \
		tmp/sirv/fastq/sample1.fastq.gz tmp/sirv/fastq/sample2.fastq.gz
	exec samtools index tmp/sirv/map-sirv.bam
	cg flair -stack 1 \
		-refseq tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		-reftranscripts tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf \
		tmp/sirv/map-sirv.bam
	# check vs expected
	exec diff tmp/sirv/isoform_counts-sqanti3-flair-sirv.genepred.tsv data/isoform_counts-sqanti3-flair-sirv.genepred.tsv
	exec diff tmp/sirv/gene_counts-sqanti3-flair-sirv.tsv data/gene_counts-sqanti3-flair-sirv.tsv
} {}

test ont_rna {isoquant basic SIRV test} {
	file delete -force tmp/sirv
	file mkdir tmp/sirv/fastq
	foreach file [glob -nocomplain data/SIRV-flames/fastq/*] {
		mklink $file tmp/sirv/fastq/[file tail $file]
	}
	foreach file [glob -nocomplain data/SIRV-flames/*] {
		if {$file eq "data/SIRV-flames/fastq"} continue
		mklink $file tmp/sirv/[file tail $file]
	}
	exec samtools faidx tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta
	cg refseq_minimap2 tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta splice
	cg map \
		-method minimap2 -preset splice -paired 0 \
		tmp/sirv/map-sirv.bam \
		tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		tmp/sirv/sirv \
		tmp/sirv/fastq/sample1.fastq.gz tmp/sirv/fastq/sample2.fastq.gz
	exec samtools index tmp/sirv/map-sirv.bam
	cg iso_isoquant -stack 1 \
		-refseq tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		-reftranscripts tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf \
		tmp/sirv/map-sirv.bam
	# check vs expected
	exec diff tmp/sirv/isoform_counts-isoquant-sirv.tsv data/isoform_counts-isoquant-sirv.tsv
	exec diff tmp/sirv/gene_counts-isoquant-sirv.tsv data/gene_counts-isoquant-sirv.tsv
} {}

test ont_rna {flames SIRV test no ref} {
	file delete -force tmp/sirv
	file mkdir tmp/sirv/fastq
	foreach file [glob -nocomplain data/SIRV-flames/fastq/*] {
		mklink $file tmp/sirv/fastq/[file tail $file]
	}
	foreach file [glob -nocomplain data/SIRV-flames/*] {
		if {$file eq "data/SIRV-flames/fastq"} continue
		mklink $file tmp/sirv/[file tail $file]
	}
	exec samtools faidx tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta
	cg refseq_minimap2 tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta splice
	cg map \
		-method minimap2 -preset splice -paired 0 \
		tmp/sirv/map-sirv.bam \
		tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		tmp/sirv/sirv \
		tmp/sirv/fastq/sample1.fastq.gz tmp/sirv/fastq/sample2.fastq.gz
	exec samtools index tmp/sirv/map-sirv.bam
	# does not seem to find new genes, so add one of each
	exec grep SIRV101 tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf > tmp/sirv/part_SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf
	foreach id {SIRV201 SIRV301 SIRV403 SIRV501 SIRV601 SIRV701} {
		exec grep $id tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf >> tmp/sirv/part_SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf
	}
	cg iso_flames \
		-refseq tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		-reftranscripts tmp/sirv/part_SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf \
		tmp/sirv/map-sirv.bam
	cg gtf2tsv tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf | cg select -s - -f {* counts-ref=1} > tmp/sirv/ref.tsv
	file delete tmp/sirv/multitranscript.tsv
	cg multitranscript -match . tmp/sirv/multitranscript.tsv tmp/sirv/isoform_counts-flames-sirv.tsv tmp/sirv/ref.tsv 
	# check vs expected
	exec diff tmp/sirv/isoform_counts-flames-sirv.tsv data/isoform_counts-flames-noref_sirv.tsv
	exec diff tmp/sirv/gene_counts-flames-sirv.tsv data/gene_counts-flames-noref_sirv.tsv
} {}

test ont_rna {flair SIRV test no ref} {
	file delete -force tmp/sirv
	file mkdir tmp/sirv/fastq
	foreach file [glob -nocomplain data/SIRV-flames/fastq/*] {
		mklink $file tmp/sirv/fastq/[file tail $file]
	}
	foreach file [glob -nocomplain data/SIRV-flames/*] {
		if {$file eq "data/SIRV-flames/fastq"} continue
		mklink $file tmp/sirv/[file tail $file]
	}
	exec samtools faidx tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta
	cg refseq_minimap2 tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta splice
	cg map \
		-method minimap2 -preset splice -paired 0 \
		tmp/sirv/map-sirv.bam \
		tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		tmp/sirv/sirv \
		tmp/sirv/fastq/sample1.fastq.gz tmp/sirv/fastq/sample2.fastq.gz
	exec samtools index tmp/sirv/map-sirv.bam
	exec grep SIRV101 tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf > tmp/sirv/part_SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf
	foreach id {SIRV201 SIRV301 SIRV403 SIRV501 SIRV601 SIRV701} {
		exec grep $id tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf >> tmp/sirv/part_SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf
	}
	cg flair -stack 1 \
		-refseq tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		-reftranscripts tmp/sirv/part_SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf \
		tmp/sirv/map-sirv.bam
	cg gtf2tsv tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf | cg select -s - -f {* counts-ref=1} > tmp/sirv/ref.tsv
	file delete tmp/sirv/multitranscript.tsv
	cg multitranscript -match . tmp/sirv/multitranscript.tsv tmp/sirv/isoform_counts-sqanti3-flair-sirv.genepred.tsv tmp/sirv/ref.tsv 
	# check vs expected
	exec diff tmp/sirv/isoform_counts-sqanti3-flair-sirv.genepred.tsv data/isoform_counts-sqanti3-flair-noref_sirv.genepred.tsv
	exec diff tmp/sirv/gene_counts-sqanti3-flair-sirv.tsv data/gene_counts-sqanti3-flair-noref_sirv.tsv
} {}

test ont_rna {isoquant SIRV test no ref} {
	file delete -force tmp/sirv
	file mkdir tmp/sirv/fastq
	foreach file [glob -nocomplain data/SIRV-flames/fastq/*] {
		mklink $file tmp/sirv/fastq/[file tail $file]
	}
	foreach file [glob -nocomplain data/SIRV-flames/*] {
		if {$file eq "data/SIRV-flames/fastq"} continue
		mklink $file tmp/sirv/[file tail $file]
	}
	exec samtools faidx tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta
	cg refseq_minimap2 tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta splice
	cg map \
		-method minimap2 -preset splice -paired 0 \
		tmp/sirv/map-sirv.bam \
		tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		tmp/sirv/sirv \
		tmp/sirv/fastq/sample1.fastq.gz tmp/sirv/fastq/sample2.fastq.gz
	exec samtools index tmp/sirv/map-sirv.bam
	# isoquant will find the "novel" genes so only giving one (not all). It will find a few more isoforms with all genes given though
	cg gtf2tsv tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf | cg select -s - > tmp/sirv/ref.tsv
	cg select -q {$transcript in "SIRV101"} tmp/sirv/ref.tsv > tmp/sirv/part_ref.tsv
	cg iso_isoquant -stack 1 \
		-refseq tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		-reftranscripts tmp/sirv/part_ref.tsv \
		tmp/sirv/map-sirv.bam
	cg gtf2tsv tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf | cg select -s - -f {* counts-ref=1} > tmp/sirv/ref.tsv
	file delete tmp/sirv/multitranscript.tsv
	cg multitranscript -match . tmp/sirv/multitranscript.tsv tmp/sirv/isoform_counts-isoquant-sirv.tsv tmp/sirv/ref.tsv 
	# check vs expected
	exec diff tmp/sirv/isoform_counts-isoquant-sirv.tsv data/isoform_counts-isoquant-noref_sirv.tsv
	exec diff tmp/sirv/gene_counts-isoquant-sirv.tsv data/gene_counts-isoquant-noref_sirv.tsv
} {}

test ont_rna {process_project multi methods} {
	file delete -force tmp/sirv tmp/sirv2 tmp/ref
	file mkdir tmp/samples/sirv/fastq
	file mkdir tmp/samples/sirv2/fastq
	file mkdir tmp/ref/sirv
	foreach file [glob -nocomplain data/SIRV-flames/fastq/*] {
		mklink $file tmp/samples/sirv/fastq/[file tail $file]
		exec cg fastq2tsv $file | cg select -q {$ROW < 1000} | cg tsv2fastq | cg bgzip > tmp/samples/sirv2/fastq/[file tail $file]
	}
	mklink data/SIRV-flames/SIRV_isoforms_multi-fasta_170612a.fasta tmp/ref/sirv/genome_sirv.ifas
	mklink data/SIRV-flames/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf tmp/ref/sirv/gene_sirv.gtf
	exec samtools faidx tmp/ref/sirv/SIRV_isoforms_multi-fasta_170612a.ifas
	cg refseq_minimap2 tmp/ref/sirv/genome_sirv.ifas splice
	exec cg process_project -stack 1 -v 2 -split 1 \
		-threads 2 \
		-paired 0 -clip 0 \
		-maxfastqdistr 250 \
		-aligner {minimap2_splice} \
		-removeduplicates 0 \
		-realign 0 \
		-distrreg 0 \
		-svcallers {} \
		-varcallers {} \
		-isocallers {isoquant flair flames} \
		-reports {} \
		-dbdir tmp/ref/sirv \
		tmp >& tmp/ontrna.log


	cg map \
		-method minimap2 -preset splice -paired 0 \
		tmp/sirv/map-sirv.bam \
		tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		tmp/sirv/sirv \
		tmp/sirv/fastq/sample1.fastq.gz tmp/sirv/fastq/sample2.fastq.gz
	exec samtools index tmp/sirv/map-sirv.bam
	cg iso_isoquant -stack 1 \
		-refseq tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		-reftranscripts tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf \
		tmp/sirv/map-sirv.bam
	# check vs expected
	exec diff tmp/sirv/isoform_counts-isoquant-sirv.tsv data/isoform_counts-isoquant-sirv.tsv
	exec diff tmp/sirv/gene_counts-isoquant-sirv.tsv data/gene_counts-isoquant-sirv.tsv
} {}

testsummarize
