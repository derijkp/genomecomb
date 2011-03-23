#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test select {-f} {
	exec cg select -f "num text" data/table.tsv
} {num	text
4	a
10	b
2	c
100	d}

test select {-q} {
	exec cg select -q {$num <= 4} data/table.tsv
} {num	text	mixed	other
4	a	a4	aaaa
2	c	a2	cc}

test select {-s} {
	exec cg select -s num data/table.tsv
} {num	text	mixed	other
2	c	a2	cc
4	a	a4	aaaa
10	b	b2	bbbbbbbbbb
100	d	aa3	dd..}

test select {-s -f} {
	exec cg select -s num -f "num mixed" data/table.tsv
} {num	mixed
2	a2
4	a4
10	b2
100	aa3}

test select {-s -q} {
	exec cg select -s num -q {$num <= 4} data/table.tsv
} {num	text	mixed	other
2	c	a2	cc
4	a	a4	aaaa}

test select {-q -f} {
	exec cg select -f "num mixed" -q {$num <= 4} data/table.tsv
} {num	mixed
4	a4
2	a2}

test select {-s -f -q} {
	exec cg select -s num  -f "num mixed" -q {$num <= 4} data/table.tsv
} {num	mixed
2	a2
4	a4}

test select {-f *} {
	exec cg select -f {chromosome begin end alleleSeq1-*} -q {$begin == 4000} data/vars1.sft
} {chromosome	begin	end	alleleSeq1-sample1	alleleSeq1-sample2
chr1	4000	4001	A	A
chr2	4000	4001	G	G}

test select {-f calculated} {
	exec cg select -f {chromosome begin end {geno1=$alleleSeq1-sample1 "/" $alleleSeq2-sample1}} -q {$begin == 4000} data/vars1.sft
} {chromosome	begin	end	geno1
chr1	4000	4001	A/G
chr2	4000	4001	G/G}

test select {-f calculated functions} {
	exec cg select -f {chromosome begin end {countG=count($alleleSeq*, == "G")}} -q {$begin == 4000} data/vars1.sft
} {chromosome	begin	end	countG
chr1	4000	4001	2
chr2	4000	4001	4}



testsummarize
