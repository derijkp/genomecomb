#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

# cg select -q 'oneof($locus,1224,1226,1598520,2816907,2816911,2818387,3038185,3042129,5019241,5054443,5211957,9605719,10679396,19364044,21233191,22632054,21542413)' cgNA19240/fannotvar-cgNA19240.tsv > testvars.tsv
# exit

test zyg {t} {
	zyg v A C A C
} t

test zyg {t} {
	zyg v C A A C
} t

test zyg {m} {
	zyg v C C A C
} m

test zyg {c} {
	zyg v C G A C
} c

test zyg {c} {
	zyg v G C A C
} c

test zyg {o} {
	zyg v G A A C
} o

test zyg {r} {
	zyg v A A A C
} r

file delete -force tmp/temp.sft

set ::env(PATH) $keeppath

testsummarize
