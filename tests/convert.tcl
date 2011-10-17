#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test select {vcf2sft} {
	exec cg vcf2sft data/test.vcf temp.tsv
	exec diff temp.tsv data/expected-test.vcf2sft
} {}

test select {vcf2sft} {
	exec cg vcf2sft data/test1000glow.vcf temp.tsv
	exec diff temp.tsv data/expected-test1000glow.vcf2sft
} {}

file delete temp.tsv

testsummarize
