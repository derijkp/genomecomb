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

testsummarize
