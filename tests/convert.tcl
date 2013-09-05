#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test select {vcf2tsv} {
	exec cg vcf2tsv data/test.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test.vcf2tsv
} {}

test select {vcf2tsv ins and del} {
	exec cg vcf2tsv data/test2.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test2.vcf2tsv
} {}

test select {vcf2tsv} {
	exec cg vcf2tsv data/test1000glow.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test1000glow.vcf2tsv
} {}

test select {bed2sft} {
	exec cg bed2sft data/sample.bed
} {#browser position chr7:127471196-127495720
#browser hide all
#track name="ColorByStrandDemo" description="Color by strand demonstration" visibility=2 colorByStrand="255,0,0 0,0,255"
chromosome	begin	end	name	score	strand
chr7	127471196  127472363  Pos1  0  +
chr7	127472363  127473530  Pos2  0  +
chr7	127473530  127474697  Pos3  0  +
chr7	127474697  127475864  Pos4  0  +
chr7	127475864  127477031  Neg1  0  -
chr7	127477031  127478198  Neg2  0  -
chr7	127478198  127479365  Neg3  0  -
chr7	127479365  127480532  Pos5  0  +
chr7	127480532  127481699  Neg4  0  -}

file delete tmp/temp.tsv

test select {sft2gff} {
	exec cg vcf2tsv data/test1000glow.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test1000glow.vcf2tsv
} {}

testsummarize
