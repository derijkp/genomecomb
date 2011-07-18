#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test select {-f} {
	file delete -force temp.tsv
	exec cg updatevarfile -c data/updatavartest.tsv temp.tsv /complgen/refseq/hg18
	exec diff temp.tsv data/expected-updatavartest.tsv
} {}

test select {-f} {
	file delete -force temp2.tsv
	exec cg select -f {chromosome begin end type ref alt} data/updatavartest.tsv temp.tsv
	exec cg updatevarfile -c temp.tsv temp2.tsv /complgen/refseq/hg18
	exec diff temp2.tsv data/expected-updatavartest2.tsv
} {}

testsummarize
