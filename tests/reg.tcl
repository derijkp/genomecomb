#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl
set keeppath $::env(PATH)
append ::env(PATH) :[file dir [file dir [file normalize [info script]]]]/bin
putsvars ::env(PATH)

test regsubtract {basic} {
	exec cg regsubtract reg1.tsv reg2.tsv
} {chromosome	begin	end
1	10	15
1	55	60
2	100	150
2	160	170
2	180	200
3	2000	2100
1-10
2-100
3-1000} error

test covered {basic} {
	exec cg covered reg1.tsv
} {chromosome	bases
1	20
2	130
3	200

total	350
1-10
2-100
3-1000} error

test covered {basic} {
	exec cg covered reg2.tsv
} {chromosome	bases
1	20
2	220
3	100

total	340
1-15
2-150
3-1000} error

test getregions {above} {
	exec getregions < coverage.tsv chr1 0 1 10 1 0
} {chromosome	begin	end
chr1	25	27
chr1	28	31}

test getregions {below} {
	exec getregions < coverage.tsv chr1 0 1 10 0 0
} {chromosome	begin	end
chr1	20	24
chr1	27	28}

test getregions {shift} {
	exec getregions < coverage.tsv chr1 0 1 10 1 -1
} {chromosome	begin	end
chr1	24	26
chr1	27	30}

set ::env(PATH) $keeppath

testsummarize
