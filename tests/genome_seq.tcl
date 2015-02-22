#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test genome_seq {basic} {
	exec cg genome_seq -i name data/reg_genome_seq.tsv /complgen/refseq/hg19_test > tmp/temp.fas  2> /dev/null
	exec diff tmp/temp.fas data/expected-reg_genome_seq.fas
} {} 

test genome_seq {basic -f -1 -n -1} {
	exec cg genome_seq -i name -f -1 -n -1 data/reg_genome_seq.tsv /complgen/refseq/hg19_test > tmp/temp.fas  2> /dev/null
	exec diff tmp/temp.fas data/expected-reg_genome_seq_fullmasked.fas
} {} 

test genome_seq {basic -p {snp 1000g}} {
	exec cg genome_seq -i name -f -1 -n -1 -p {snp 1000g} data/reg_genome_seq.tsv /complgen/refseq/hg19_test > tmp/temp.fas  2> /dev/null
	exec diff tmp/temp.fas data/expected-reg_genome_seq_fullmasked.fas
} {} 

test genome_seq {basic -p {1000g}} {
	exec cg genome_seq -i name -f -1 -n -1 -p {1000g} data/reg_genome_seq2.tsv /complgen/refseq/hg19_test > tmp/temp.fas  2> /dev/null
	exec diff tmp/temp.fas data/expected-reg_genome_seq2.fas
} {} 

test genome_seq {outfile} {
	exec cg genome_seq -i name data/reg_genome_seq.tsv /complgen/refseq/hg19_test tmp/temp.fas  2> /dev/null
	exec diff tmp/temp.fas data/expected-reg_genome_seq.fas
} {} 

test genome_seq {concat and makemap} {
	exec cg genome_seq -c "" -m tmp/temp.map -i name data/reg_genome_seq.tsv /complgen/refseq/hg19_test > tmp/temp.fas  2> /dev/null
	exec diff tmp/temp.fas data/expected-reg_genome_seq.ifas
	exec diff tmp/temp.map data/expected-reg_genome_seq.map
} {} 

test genome_seq {gcsplit} {
	exec cg genome_seq -g 50 -gs 60 -i name data/reg_genome_seq.tsv /complgen/refseq/hg19_test tmp/temp.fas  2> /dev/null
	exec diff tmp/temp-lowgc.fas data/expected-reg_genome_seq-lowgc.fas
	exec diff tmp/temp-highgc.fas data/expected-reg_genome_seq-highgc.fas
} {} 

test genome_seq {-g 0} {
	exec cg genome_seq -g 0 data/reg_genome_seq.tsv /complgen/refseq/hg19_test tmp/temp.fas  2> /dev/null
	exec diff tmp/temp.fas data/expected-reg_genome_seq-noid-g0.fas
} {}

test genome_seq {-g 0 -gd 0} {
	exec cg genome_seq -g 0 -gd 0 data/reg_genome_seq.tsv /complgen/refseq/hg19_test tmp/temp.fas  2> /dev/null
	exec diff tmp/temp.fas data/expected-reg_genome_seq-noid.fas
} {}

test genome_seq {gcsplit and -l} {
	exec cg genome_seq -g 50 -gs 60 -l _ -i name data/reg_genome_seq.tsv /complgen/refseq/hg19_test tmp/temp.fas  2> /dev/null
	exec diff tmp/temp-highgc.fas data/expected-reg_genome_seq-highgc-l.fas
} {} 

test genome_seq {split} {
	file delete {*}[glob tmp/*.fas]
	exec cg genome_seq -g 50 -s 1 -i name data/reg_genome_seq.tsv /complgen/refseq/hg19_test tmp/  2> /dev/null
	exec diff tmp/r1.fas data/expected-reg_genome_seq-1.fas
	exec diff tmp/r9.fas data/expected-reg_genome_seq-9.fas
	lsort -dict [glob tmp/*.fas]
} {tmp/r1.fas tmp/r2.fas tmp/r3.fas tmp/r4.fas tmp/r5.fas tmp/r6.fas tmp/r7.fas tmp/r8.fas tmp/r9.fas tmp/r10.fas tmp/r11.fas tmp/r12.fas tmp/r13.fas tmp/r14.fas tmp/r15.fas tmp/r16.fas tmp/r17.fas} 

test genome_seq {gcsplit split} {
	file delete -force {*}[glob tmp/*]
	exec cg genome_seq -g 50 -s 1 -gs 60 -i name data/reg_genome_seq.tsv /complgen/refseq/hg19_test tmp/temp.fas  2> /dev/null
	exec diff tmp/temp-lowgc/r1.fas data/expected-reg_genome_seq-1.fas
	exec diff tmp/temp-highgc/r9.fas data/expected-reg_genome_seq-9.fas
	list [lsort -dict [glob tmp/temp-lowgc/*.fas]] [lsort -dict [glob tmp/temp-highgc/*.fas]]
} {{tmp/temp-lowgc/r1.fas tmp/temp-lowgc/r2.fas tmp/temp-lowgc/r4.fas tmp/temp-lowgc/r5.fas tmp/temp-lowgc/r6.fas tmp/temp-lowgc/r7.fas tmp/temp-lowgc/r8.fas tmp/temp-lowgc/r10.fas tmp/temp-lowgc/r11.fas tmp/temp-lowgc/r15.fas tmp/temp-lowgc/r16.fas} {tmp/temp-highgc/r3.fas tmp/temp-highgc/r9.fas tmp/temp-highgc/r12.fas tmp/temp-highgc/r13.fas tmp/temp-highgc/r14.fas tmp/temp-highgc/r17.fas}} 

test primercheck {basic} {
	exec cg primercheck data/primers.tsv /complgen/refseq/hg19_test tmp/primerinfo.tsv 2> /dev/null 2> /dev/null
	exec diff tmp/primerinfo.tsv data/primercheck-results.tsv
} {} 

file delete -force tmp/temp.sft

set ::env(PATH) $keeppath

testsummarize
