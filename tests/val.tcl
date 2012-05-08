#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

# cg select -q 'oneof($locus,1224,1226,1598520,2816907,2816911,2818387,3038185,3042129,5019241,5054443,5211957,9605719,10679396,19364044,21233191,22632054,21542413)' cgNA19240/fannotvar-cgNA19240.tsv > testvars.tsv

puts "expected data no longer valid with new dbsnp"
# exit

test makesequenom {basic} {
	exec cg makesequenom data/testvars.tsv temp.sft /complgen/refseq/hg19 2> /dev/null
	exec diff temp.sft data/expected-makesequenom-testvars.tsv
} {} 

test makesequenom {basic} {
	exec cg makesequenom -f -1 -n -1 data/testvars.tsv temp.sft /complgen/refseq/hg19 2> /dev/null
	exec diff temp.sft data/expected-makesequenom-1-1-testvars.tsv
} {} 

test genome_seq {basic} {
	exec cg genome_seq -i name data/reg_genome_seq.tsv /complgen/refseq/hg19 > temp.sft  2> /dev/null
	exec diff temp.sft data/expected-reg_genome_seq.tsv
} {} 

test makesequenom {basic} {
	exec cg genome_seq -f -1 -n -1 data/reg_genome_seq.tsv /complgen/refseq/hg19 > temp.sft  2> /dev/null
	exec diff temp.sft data/expected-reg_genome_seq-1-1.tsv
} {} 

file delete -force temp.sft

set ::env(PATH) $keeppath

testsummarize
