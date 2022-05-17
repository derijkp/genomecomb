#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test flames {basic SIRV test} {
	cd $::smalltestdir
	file delete -force tmp/sirv_flames
	file mkdir tmp/sirv_flames/fastq
	foreach file [glob -nocomplain ori/SIRV/data/fastq/*] {
		mklink $file tmp/sirv_flames/fastq/[file tail $file]
	}
	foreach file [glob -nocomplain ori/SIRV/data/*] {
		if {$file eq "ori/SIRV/data/fastq"} continue
		mklink $file tmp/sirv_flames/[file tail $file]
	}
	cd tmp/sirv_flames
#	puts time:[time {
#		exec bulk_long_pipeline.py \
#			--gff3 SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf \
#			--genomefa SIRV_isoforms_multi-fasta_170612a.fasta \
#			--outdir FLAMES_output \
#			--config_file SIRV_config.json \
#			--fq_dir fastq \
#			>@ stdout 2>@ stderr
#	}]
	cg refseq_minimap2 /home/peter/genomecomb.smalltestdata/tmp/sirv_flames/SIRV_isoforms_multi-fasta_170612a.fasta splice
	cg map \
		-method minimap2 -preset splice -paired 0 \
		map-sirv_flames.bam \
		SIRV_isoforms_multi-fasta_170612a.fasta \
		sirv_flames \
		fastq/sample1.fastq.gz fastq/sample2.fastq.gz
	exec samtools index map-sirv_flames.bam
	cg flames \
		-refseq SIRV_isoforms_multi-fasta_170612a.fasta \
		-reftranscripts SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf \
		map-sirv_flames.bam
				
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 -x *.bam -x *.bai \
		tmp/sirv_flames expected/sirv_flames]
	join [list_remove $result {}] \n
} {}

testsummarize

