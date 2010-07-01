#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

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

testsummarize
