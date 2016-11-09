#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

# cg select -q 'oneof($locus,1224,1226,1598520,2816907,2816911,2818387,3038185,3042129,5019241,5054443,5211957,9605719,10679396,19364044,21233191,22632054,21542413)' cgNA19240/fannotvar-cgNA19240.tsv > testvars.tsv
# exit

test makeprimers {basic} {
	exec cg makeregions data/testvars.tsv 200 > tmp/valregs.tsv
	exec cg makeprimers tmp/valregs.tsv 600 500 $::refseqdir/hg18_test > tmp/primersvalregs.tsv
	exec diff tmp/primersvalregs.tsv data/makeprimers-results.tsv
} {} 

test makeprimers {basic with minfreq} {
	exec cg makeregions data/testvars.tsv 200 > tmp/valregs.tsv
	exec cg makeprimers tmp/valregs.tsv 600 500 $::refseqdir/hg18_test 0.5 > tmp/primersvalregs.tsv
	exec diff tmp/primersvalregs.tsv data/makeprimers-results2.tsv
} {} 

test makesequenom {basic} {
	exec cg makesequenom data/testvars.tsv tmp/temp.sft $::refseqdir/hg19_test
	exec diff tmp/temp.sft data/expected-makesequenom-testvars.tsv
} {} 

test makesequenom {basic} {
	exec cg makesequenom -f -1 -n -1 data/testvars.tsv tmp/temp.sft $::refseqdir/hg19_test
	exec diff tmp/temp.sft data/expected-makesequenom-1-1-testvars.tsv
} {} 

file delete -force tmp/temp.sft

set ::env(PATH) $keeppath

testsummarize
