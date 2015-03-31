#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test remap {basic} {
	cg remap data/vars1.sft data/remap.tsv tmp/temp.tsv
	list [cg select -f {chromosome begin end} tmp/temp.tsv] [cg select -f {chromosome begin end} tmp/temp.tsv.unmapped]
} {{chromosome	begin	end
chr2	1001	1002
chr2	2000	2010
chr3	4001	4002
chr3	4002	4003
chr3	4100	4101
chr3	5001	5001
chr3	5001	5002
chr3	5006	5007
chr3	5011	5012
chr3	5012	5013
chr3	8001	8002} {chromosome	begin	end
chr1	4000	4001
chr1	4099	4100
chr1	5020	5021}}

testsummarize
