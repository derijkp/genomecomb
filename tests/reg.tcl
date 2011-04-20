#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test multireg {basic} {
	file delete temp.tsv
	exec cg multireg temp.tsv data/reg1.tsv data/reg2.tsv 2> /dev/null
	exec diff temp.tsv data/expected-multireg-reg1-reg2.sft
} {}

test multireg {same} {
	file delete temp.tsv
	exec cg multireg temp.tsv data/reg1.tsv data/reg1b.tsv 2> /dev/null
	exec diff temp.tsv data/expected-multireg-reg1-reg1b.sft
} {}

test multireg {add empty} {
	file delete temp.tsv
	exec cg multireg temp.tsv data/reg1b.tsv data/regempty.tsv 2> /dev/null
	file_read temp.tsv
} {chromosome	begin	end	reg1b	regempty
1	10	20	1	0
1	50	60	1	0
}

test multireg {add empty first} {
	file delete temp.tsv
	exec cg multireg temp.tsv data/regempty.tsv data/reg1b.tsv 2> /dev/null
	file_read temp.tsv
} {chromosome	begin	end	regempty	reg1b
1	10	20	0	1
1	50	60	0	1
}

test multireg {add fully empty} {
	file delete temp.tsv
	exec cg multireg temp.tsv data/reg1b.tsv data/empty.tsv
} {not a region file} error regexp

test multireg {3 adds} {
	file delete temp.tsv
	exec cg multireg temp.tsv data/reg1.tsv data/reg1b.tsv data/reg2.tsv 2> /dev/null
	exec diff temp.tsv data/expected-multireg-reg1-reg1b-reg2.sft
} {}

test multireg {different chromosome naming} {
	file delete temp.tsv
	exec cg multireg temp.tsv data/reg1.tsv data/reg4.tsv 2> /dev/null
	exec diff temp.tsv data/expected-multireg-reg1-reg4.sft
} {}

test regsubtract {basic} {
	exec cg regsubtract data/reg1.tsv data/reg2.tsv
} {chromosome	begin	end
1	10	15
1	55	60
2	100	150
2	160	170
2	180	200
3	2000	2100
Y	1000	1010
Y	1900	2000
1-10
2-100
3-1000
M-10
X-100
Y-1000} error

test covered {basic} {
	exec cg covered data/reg1.tsv
} {chromosome	bases
1	20
2	130
3	200
M	10
X	100
Y	1000

total	1460
1-10
2-100
3-1000
M-10
X-100
Y-1000} error

test covered {basic} {
	exec cg covered data/reg2.tsv
} {chromosome	bases
1	20
2	220
3	100
M	15
X	110
Y	890

total	1355
1-15
2-150
3-1000
M-10
X-90
Y-1010} error

test getregions {above} {
	exec getregions < data/coverage-chr1.tsv chr1 0 1 10 1 0
} {chr1	25	27
chr1	28	31
chr1	40	42}

test getregions {below} {
	exec getregions < data/coverage-chr1.tsv chr1 0 1 10 0 0
} {chr1	20	24
chr1	27	28
chr1	42	43}

test getregions {shift} {
	exec getregions < data/coverage-chr1.tsv chr1 0 1 10 1 -1
} {chr1	24	26
chr1	27	30
chr1	39	41}

test cg_regextract {above} {
	exec cg regextract -qfields coverage -above 1 -shift 0 10 data/coverage-chr1.tsv
} {chromosome	begin	end
chr1	25	27
chr1	28	31
chr1	40	42
Processing data/coverage-chr1.tsv} error

test cg_regextract {below} {
	exec cg regextract -qfields coverage -above 0 -shift 0 10 data/coverage-chr1.tsv
} {chromosome	begin	end
chr1	20	24
chr1	27	28
chr1	42	43
Processing data/coverage-chr1.tsv} error

test cg_regextract {shift} {
	exec cg regextract -qfields coverage -above 1 -shift -1 10 data/coverage-chr1.tsv
} {chromosome	begin	end
chr1	24	26
chr1	27	30
chr1	39	41
Processing data/coverage-chr1.tsv} error

test regjoin {basic} {
	exec cg regjoin data/reg1.tsv data/reg2.tsv 2> /dev/null
} {chromosome	begin	end
1	10	25
1	45	60
2	100	200
2	300	500
3	1000	1100
3	2000	2100
M	10	25
X	90	200
Y	1000	2000}

test regjoin {basic} {
	exec cg regjoin data/reg1.tsv data/reg3.tsv 2> /dev/null
} {chromosome	begin	end
1	5	8
1	10	25
1	45	60
1	70	80
2	90	250
2	300	500
3	100	250
3	1000	2100
3	2500	2600
M	10	20
X	100	200
Y	1000	2000}

test regjoin {self} {
	exec cg regjoin data/reg3.tsv 2> /dev/null
} {chromosome	begin	end
1	5	8
1	15	18
1	19	25
1	45	50
1	70	80
2	90	100
2	200	250
2	300	500
3	100	250
3	1100	2000
3	2500	2600}

set ::env(PATH) $keeppath

file delete -force temp.tsv temp.tsv.old

testsummarize
