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

	puts time:[time {
	exec bulk_long_pipeline.py \
		--gff3 SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf \
		--genomefa SIRV_isoforms_multi-fasta_170612a.fasta \
		--outdir FLAMES_output \
		--config_file SIRV_config.json \
		--fq_dir fastq \
		>@ stdout 2>@ stderr
	}]	

	puts time:[time {
		exec cg sirv_flames {*}$::dopts \
			-distrreg chr \
			-threads 6 \
			-refseq $::refseqdir/hg19 \
			tmp/sirv_flames/fast5 tmp/sirv_flames/fastq tmp/sirv_flames/map-hlongshot-sminimap2-methtest.bam \
			tmp/sirv_flames/meth-methtest.tsv.zst \
			>& tmp/sirv_flames/sirv_flames.log
	}]
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 -x *.bam -x *.bai -x fastqc_report.html \
		-x colinfo -x meth.html -x *.zsti -x *.lz4i -x *.finished -x info_analysis.tsv \
		-x *.analysisinfo -x *.png -x *.submitting \
		-x *log_jobs -x *.index -x *.log \
		tmp/sirv_flames expected/sirv_flames]
	join [list_remove $result {}] \n
} {}

testsummarize

