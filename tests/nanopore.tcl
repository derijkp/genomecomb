#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test nano {basecaller_albacore basic} {
	file mkdir tmp/test/group/fast5/0
	file copy {*}[glob data/fast5/*.fast5] tmp/test//group/fast5/0
	cg basecaller_albacore -stack 1 -v 2 -subdirs 1 -numreads 4 -c r95_450bps_linear.cfg -threads 2 tmp/fastq tmp/test/group
	exec gunzip tmp/fastq/pass_group_0.fastq.gz
	exec gunzip -c data/expected-pass_group_0.fastq.gz > tmp/expected.fastq
	exec diff tmp/fastq/pass_group_0.fastq tmp/expected.fastq
} {}

testsummarize

