#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

# cg select -q 'oneof($locus,1224,1226,1598520,2816907,2816911,2818387,3038185,3042129,5019241,5054443,5211957,9605719,10679396,19364044,21233191,22632054,21542413)' cgNA19240/fannotvar-cgNA19240.tsv > testvars.tsv
# exit

test makeprimers {basic} {
	exec cg makeregions data/testvars.tsv 200 > tmp/valregs.tsv
	exec cg makeprimers tmp/valregs.tsv 600 500 /complgen/refseq/hg18 > tmp/primersvalregs.tsv 2> /dev/null
	exec diff tmp/primersvalregs.tsv data/makeprimers-results.tsv
} {} 

test makeprimers {basic with minfreq} {
	exec cg makeregions data/testvars.tsv 200 > tmp/valregs.tsv
	exec cg makeprimers tmp/valregs.tsv 600 500 /complgen/refseq/hg18 0.5 > tmp/primersvalregs.tsv 2> /dev/null
	exec diff tmp/primersvalregs.tsv data/makeprimers-results2.tsv
} {} 

test makesequenom {basic} {
	exec cg makesequenom data/testvars.tsv tmp/temp.sft /complgen/refseq/hg19_test 2> /dev/null
	exec diff tmp/temp.sft data/expected-makesequenom-testvars.tsv
} {} 

test makesequenom {basic} {
	exec cg makesequenom -f -1 -n -1 data/testvars.tsv tmp/temp.sft /complgen/refseq/hg19_test 2> /dev/null
	exec diff tmp/temp.sft data/expected-makesequenom-1-1-testvars.tsv
} {} 

test genome_seq {basic} {
	exec cg genome_seq -i name data/reg_genome_seq.tsv /complgen/refseq/hg19_test > tmp/temp.fas  2> /dev/null
	exec diff tmp/temp.fas data/expected-reg_genome_seq.fas
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

test primercheck {basic} {
	exec cg primercheck data/primers.tsv /complgen/refseq/hg19_test tmp/primerinfo.tsv 2> /dev/null 2> /dev/null
	exec diff tmp/primerinfo.tsv data/primercheck-results.tsv
} {} 

file delete -force tmp/temp.sft

set ::env(PATH) $keeppath

testsummarize
