#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test select {-f} {
	file delete -force tmp/temp.tsv
	exec cg updatevarfile -c data/updatavartest.tsv tmp/temp.tsv /complgen/refseq/hg18
	exec diff tmp/temp.tsv data/expected-updatavartest.tsv
} {}

test select {-f} {
	file delete -force tmp/temp2.tsv
	exec cg select -f {chromosome begin end type ref alt} data/updatavartest.tsv tmp/temp.tsv
	exec cg updatevarfile -c tmp/temp.tsv tmp/temp2.tsv /complgen/refseq/hg18
	exec diff tmp/temp2.tsv data/expected-updatavartest2.tsv
} {}

testsummarize
