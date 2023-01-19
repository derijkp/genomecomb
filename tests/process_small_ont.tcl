#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl
set keepdir [pwd]

if {![info exists argv]} {set argv {}}

set test_cleantmp 0
set download 0
set run 0
set forcebenchmark 0
set distr 1
cg_options job_example argv {
	-run {
		set run $value
	}
	-distr - -d {
		set distr $value
	}
	-download {
		set download $value
	}
	-forcebenchmark {
		set forcebenchmark $value
	}
}

# check or run

proc benchmarkvars {args} {
	set target [lindex $args end]
	if {!$::forcebenchmark && [file exists $target]} {
		puts "skipping $target: already exists and no forcebenchmark"
		return $target
	}
	exec cg benchmarkvars {*}$args
}

# tests
# =====

# tests will run in ~/genomecomb.smalltestdata
# this directy should be present, and should have the following directories in it
# - ori with the original data to work from
# - expected with the expected results to compare to
# This directory can be made from scratch using make_testsdata_ont.tcl

test process_ontrna {ontrna} {
	cd $::smalltestdir
	if {$::run} {
		puts "Running ontrna"
		file delete -force tmp/ontrna
		file mkdir tmp/ontrna/samples
		foreach sample {
			HG001_NA12878_cDNA HG001_NA12878_directRNA HG001_NA12878_ivtRNA
		} {
			file mkdir tmp/ontrna/samples/$sample/fastq
			foreach file [glob ori/nanopore-wgs-consortium-rna/$sample/splitfastq/*fastq.gz] {
				mklink $file tmp/ontrna/samples/$sample/fastq/[file tail $file]
			}
		}
		# run
		exec cg process_project -stack 1 -v 2 -split 1 \
			-d $::distr \
			-threads 6 \
			-paired 0 -clip 0 \
			-maxfastqdistr 250 \
			-aligner {minimap2_splice} \
			-removeduplicates 0 \
			-realign 0 \
			-distrreg chr \
			-svcallers {} \
			-varcallers {} \
			-isocallers {isoquant flair} \
			-reports {-fastqc predictgender} \
			-dbdir /complgen/refseq/hg38 \
			tmp/ontrna >& tmp/ontrna.log
		puts "ontrna run started"
	} else {
		# check vs expected
		puts "Checking ontrna"
		set result {}
		# flair has some non-deterministic results -> have to exclude most of its results from comparison for now
		lappend result [tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x colinfo -x fastqc_report.html \
			-x *bam.dupmetrics -x info_analysis.tsv -x *.zsti -x *.lz4i -x *.finished -x *.index -x info_analysis.tsv \
			-x *.analysisinfo -x *.png -x *.submitting \
			-x GMST -x sqanti3-flair-* -x gene_counts-flair-* \
			-ignorefields {
				associated_gene
			} \
			tmp/ontrna expected/ontrna]
		join [list_remove $result {}] \n
	}
} {}

if {$run} {
	puts "All tests started"
} else {
	testsummarize
}
