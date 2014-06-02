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

test select {-rf} {
	exec cg select -rf "mixed other" data/table.tsv
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

test select {-q multiple lines, tabs} {
exec cg select -s num  -f "num mixed" -q {
	$num <= 4
	&& $text == "c"
} data/table.tsv
} {num	mixed
2	a2}

test select {default -s -} {
	exec cg select -s - data/expected-test1000glow.vcf2tsv tmp/result.tsv
	exec diff data/expected-test1000glow.vcf2tsv tmp/result.tsv
} {}

test select {-f *} {
	exec cg select -f {chromosome begin end alleleSeq1-*} -q {$begin == 4000} [gzfile data/vars1.sft]
} {chromosome	begin	end	alleleSeq1-sample1	alleleSeq1-sample2
chr1	4000	4001	A	A
chr2	4000	4001	G	G}

test select {-f calculated} {
	exec cg select -f {chromosome begin end {geno1="$alleleSeq1-sample1/$alleleSeq2-sample1"}} -q {$begin == 4000} [gzfile data/vars1.sft]
} {chromosome	begin	end	geno1
chr1	4000	4001	A/G
chr2	4000	4001	G/G}

test select {-f calculated functions} {
	exec cg select -f {chromosome begin end {countG=count($alleleSeq*, == "G")}} -q {$begin == 4000} [gzfile data/vars1.sft]
} {chromosome	begin	end	countG
chr1	4000	4001	2
chr2	4000	4001	4}

test select {-f wildcard error} {
	exec cg select -f {chromosome begin end {countG=count($alleleSeq-*, == "G")}} -q {$begin == 4000} [gzfile data/vars1.sft]
} {field "alleleSeq-*" not present in file and no sampledata file found} error

test select {-f calculated if} {
	exec cg select -f {chromosome begin end {countG=if(($alleleSeq1-sample1 == "A" || $alleleSeq2-sample1 == "A"),"hasA","noA")}} -q {$begin == 4000} [gzfile data/vars1.sft]
} {chromosome	begin	end	countG
chr1	4000	4001	hasA
chr2	4000	4001	noA}

test select {-f calculated if} {
	exec cg select -f {chromosome begin end {countG=if(count($alleleSeq*, == "G")<4,"<4",">=4")}} -q {$begin == 4000} [gzfile data/vars1.sft]
} {chromosome	begin	end	countG
chr1	4000	4001	<4
chr2	4000	4001	>=4}

test select {keep header info and format rtg: -hc} {
	exec cg select -hc 1 -s position data/rtgsnps.tsv tmp/temp.tsv
	file delete temp
	catch {exec diff tmp/temp.tsv data/rtgsnps.tsv > temp}
	file_read temp
} {3c3
< name	position	type	reference	prediction	posterior	coverage	correction	support_statistics
---
> #name	position	type	reference	prediction	posterior	coverage	correction	support_statistics
}

test select {keep header info and format rtg: -hc} {
	exec cg select -hc 2 -s position data/rtgsnps.tsv tmp/temp.tsv
	file delete temp
	catch {exec diff tmp/temp.tsv data/rtgsnps.tsv > temp}
	file_read temp
} {}

test select {-hc with chars that must be escaped} {
	write_tab tmp/temp.tsv {
		{# a comment with ;}
		#a {b c}
		1 {2 3}
	}
	write_tab tmp/expected.tsv {
		{# a comment with ;}
		{b c}
		{2 3}
	}
	exec cg select -hc 1 -f {{b c}} tmp/temp.tsv tmp/out.tsv
	exec diff tmp/out.tsv tmp/expected.tsv
} {}

test select {keep header info and format vcf} {
	exec cg select -s POS data/test.vcf tmp/temp.tsv
	file delete temp
	catch {exec diff tmp/temp.tsv data/test.vcf > temp}
	file_read temp
} {18c18
< #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003
---
> #CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
}

# cannot sort on calculated fields (yet)
#test select {-f calculated functions + sort} {
#	exec cg select -f {chromosome begin end {countG=count($alleleSeq*, == "G")}} -q {$begin == 4000} -s countG [gzfile data/vars1.sft]
#} {chromosome	begin	end	countG
#chr1	4000	4001	2
#chr2	4000	4001	4}

test groupby {groupby} {
	exec cg groupby s1 < data/table2.tsv
} {s1	pos	s2	s3
u	1	u	u
vx	2	u	u
v	3,4,5	u,v,v	v,u,v
u	6,7,4000000000	v,v,xx	u,v,xx}

test groupby {groupby -sumfields} {
	exec cg groupby -sumfields pos s2 < data/table2.tsv
} {s2	pos	s1	s3
u	6	u,vx,v	u,u,v
v	22	v,v,u,u	u,v,u,v
xx	4000000000	u	xx}

test groupby {groupby -sumfields} {
	exec cg groupby -sumfields pos s1 < data/table2.tsv
} {s1	pos	s2	s3
u	1	u	u
vx	2	u	u
v	12	u,v,v	v,u,v
u	4000000013	v,v,xx	u,v,xx}

test groupby {groupby -stats} {
	exec cg groupby -stats pos s2 < data/table2.tsv
} {s2	s1	s3	pos_min	pos_total	pos_count	pos_max
u	u,vx,v	u,u,v	1	6	3	3
v	v,v,u,u	u,v,u,v	4	22	4	7
xx	u	xx	4000000000	4000000000	1	4000000000}

test groupby {groupby -sumfields one count (non existing)} {
	exec cg groupby -f s1 -sumfields count s1 < data/table2.tsv
} {s1	count
u	1
vx	1
v	3
u	3}

test groupby {groupby -sumfields one count (non existing) -sorted 0} {
	exec cg groupby -sorted 0 -f s1 -sumfields count s1 < data/table2.tsv
} {s1	count
u	4
v	3
vx	1}

test groupby {groupby -sorted 0} {
	exec cg groupby -sorted 0 s1 < data/table2.tsv
} {s1	pos	s2	s3
u	1,6,7,4000000000	u,v,v,xx	u,u,v,xx
v	3,4,5	u,v,v	v,u,v
vx	2	u	u}

test groupby {groupby -sorted 0 -f pos} {
	exec cg groupby -sorted 0 -f pos s1 < data/table2.tsv
} {s1	pos
u	1,6,7,4000000000
v	3,4,5
vx	2}

test groupby {groupby -sorted 0 -sumfields} {
	exec cg groupby -sorted 0 -sumfields pos s2 < data/table2.tsv
} {s2	pos	s1	s3
u	6	u,vx,v	u,u,v
v	22	v,v,u,u	u,v,u,v
xx	4000000000	u	xx}

test groupby {groupby -sorted 0 -sumfields} {
	exec cg groupby -sorted 0 -sumfields pos s1 < data/table2.tsv
} {s1	pos	s2	s3
u	4000000014	u,v,v,xx	u,u,v,xx
v	12	u,v,v	v,u,v
vx	2	u	u}

test groupby {groupby -sorted 0 -stats} {
	exec cg groupby -sorted 0 -stats pos s2 < data/table2.tsv
} {s2	s1	s3	pos_min	pos_total	pos_count	pos_max
u	u,vx,v	u,u,v	1	6	3	3
v	v,v,u,u	u,v,u,v	4	22	4	7
xx	u	xx	4000000000	4000000000	1	4000000000}

test select {sm} {
	split [exec cg select -f {chromosome begin} -q {sm(sample1,sample2)} < data/vars1.sft] \n
} {{chromosome	begin} {chr1	4000} {chr1	4001} {chr1	4099} {chr2	4000} {chr2	4001} {chr2	4099} {chr2	5010}}

test select {df} {
	split [exec cg select -f {chromosome begin} -q {df(sample1,sample2)} < data/vars1.sft] \n
} {{chromosome	begin} {chr1	5020} {chr2	5000} {chr2	5011} {chr2	8000}}

test select {mm} {
	split [exec cg select -f {chromosome begin} -q {mm(sample1,sample2)} < data/vars1.sft] \n
} {{chromosome	begin} {chr2	5005}}

test select {-q count with wildcard} {
	exec cg select -f {chromosome begin end alleleSeq*} -q {count($alleleSeq*, == "G") > 3} [gzfile data/vars1.sft]
} {chromosome	begin	end	alleleSeq1-sample1	alleleSeq2-sample1	alleleSeq1-sample2	alleleSeq2-sample2
chr1	4001	4002	G	G	G	G
chr2	4000	4001	G	G	G	G
chr2	4001	4002	G	G	G	G}

test select {lmin column} {
	split [exec cg select -f {chromosome begin {lmin=lmin($list)}} < data/vars1.sft] \n
} {{chromosome	begin	lmin} {chr1	4000	4} {chr1	4001	1} {chr1	4099	1} {chr1	5000	1} {chr1	5020	3} {chr2	4000	2} {chr2	4001	2} {chr2	4099	2} {chr2	5000	2} {chr2	5000	3} {chr2	5005	4} {chr2	5010	20} {chr2	5011	NaN} {chr2	8000	NaN}}

test select {lmind column} {
	split [exec cg select -f {chromosome begin {lmin=lmind($list,10)}} < data/vars1.sft] \n
} {{chromosome	begin	lmin} {chr1	4000	4} {chr1	4001	1} {chr1	4099	1} {chr1	5000	1} {chr1	5020	3} {chr2	4000	2} {chr2	4001	2} {chr2	4099	2} {chr2	5000	2} {chr2	5000	3} {chr2	5005	4} {chr2	5010	10} {chr2	5011	10} {chr2	8000	10}}

test select {lmind select} {
	split [exec cg select -f {chromosome begin} -q {lmind($list,10) == 2} < data/vars1.sft] \n
} {{chromosome	begin} {chr2	4000} {chr2	4001} {chr2	4099} {chr2	5000}}

test select {lmax column} {
	split [exec cg select -f {chromosome begin {lmax=lmax($list)}} < data/vars1.sft] \n
} {{chromosome	begin	lmax} {chr1	4000	4} {chr1	4001	4} {chr1	4099	2} {chr1	5000	2} {chr1	5020	3} {chr2	4000	2} {chr2	4001	4} {chr2	4099	4} {chr2	5000	2} {chr2	5000	3} {chr2	5005	4} {chr2	5010	20} {chr2	5011	NaN} {chr2	8000	NaN}}

test select {lmaxd column} {
	split [exec cg select -f {chromosome begin {lmax=lmaxd($list,0)}} < data/vars1.sft] \n
} {{chromosome	begin	lmax} {chr1	4000	4} {chr1	4001	4} {chr1	4099	2} {chr1	5000	2} {chr1	5020	3} {chr2	4000	2} {chr2	4001	4} {chr2	4099	4} {chr2	5000	2} {chr2	5000	3} {chr2	5005	4} {chr2	5010	20} {chr2	5011	0} {chr2	8000	0}}

test select {lmaxd select} {
	split [exec cg select -f {chromosome begin} -q {lmaxd($list,10) == 2} < data/vars1.sft] \n
} {{chromosome	begin} {chr1	4099} {chr1	5000} {chr2	4000} {chr2	5000}}

test select {counthasone column} {
	split [exec cg select -f {chromosome begin {lmin=counthasone($list, ==2)}} < data/vars1.sft] \n
} {{chromosome	begin	lmin} {chr1	4000	0} {chr1	4001	1} {chr1	4099	1} {chr1	5000	1} {chr1	5020	0} {chr2	4000	1} {chr2	4001	1} {chr2	4099	1} {chr2	5000	1} {chr2	5000	0} {chr2	5005	0} {chr2	5010	0} {chr2	5011	0} {chr2	8000	0}}

test select {ROW} {
	split [exec cg select -f {chromosome begin end} -q {$ROW == 2} < data/vars1.sft] \n
} {{chromosome	begin	end} {chr1	4099	5000}}

test select {ROW} {
	split [exec cg select -f {chromosome begin end ROW} -q {$ROW >= 2 && $ROW <= 3} < data/vars1.sft] \n
} {{chromosome	begin	end	ROW} {chr1	4099	5000	2} {chr1	5000	5010	3}}

test select {shared objects bugcheck} {
	split [exec cg select -f {chromosome begin end {test=[set ::keep $begin]}} -q {$ROW between {2 3}} < data/vars1.sft] \n
} {{chromosome	begin	end	test} {chr1	4099	5000	4099} {chr1	5000	5010	5000}}

test select {shared objects bugcheck} {
	split [exec cg select -f {chromosome begin end {test="[get ::keep 1]-[set ::keep $begin]"}} -q {$ROW between {2 3}} < data/vars1.sft] \n
} {{chromosome	begin	end	test} {chr1	4099	5000	1-4099} {chr1	5000	5010	4099-5000}}

test select {start brace bugcheck} {
	split [exec cg select -f {chromosome begin end {type=($type == "snp")? "Snp" : (($type == "del")? "Deletion" : $type)}} -q {$ROW in {2 3 9 10}} < data/vars1.sft] \n
} {{chromosome	begin	end	type} {chr1	4099	5000	Snp} {chr1	5000	5010	Deletion} {chr2	5000	5010	ins} {chr2	5005	5006	Snp}}

test select {list @-} {
	split [exec cg select -f {begin {diff=vdef($freq-sample1,0) @- vdef($freq-sample2,0)}} < data/vars3.sft] \n
} {{begin	diff} {4000	0.09999999999999998} {4001	0.09999999999999998,-0.1} {4099	0.3,-0.5} {5020	-0.1,0.5,0.0} {4001	-0.9,0.8} {4099	0.2,0.4}}

test select {list @- @-} {
	split [exec cg select -f {begin {diff=vdef($freq-sample1,0) @- vdef($freq-sample2,0) @- vdef($freq-sample2,0)}} < data/vars3.sft] \n
} {{begin	diff} {4000	-0.30000000000000004} {4001	-0.30000000000000004,-0.30000000000000004} {4099	0.09999999999999998,-1.1} {5020	-0.30000000000000004,0.5,-0.2} {4001	-1.8,0.8} {4099	0.2,0.4}}

test select {list @- @<} {
	split [exec cg select \
		-f {begin {diff=vdef($freq-sample1,0) @- vdef($freq-sample2,0)}} \
		-q {lone((vdef($freq-sample1,0) @- vdef($freq-sample2,0)) @< 0)} \
		< data/vars3.sft] \n
} {{begin	diff} {4001	0.09999999999999998,-0.1} {4099	0.3,-0.5} {5020	-0.1,0.5,0.0} {4001	-0.9,0.8}}

test select {list @- NaN} {
	split [exec cg select -f {begin {diff=$freq-sample1 @- $freq-sample2}} < data/vars3.sft] \n
} {{begin	diff} {4000	0.09999999999999998} {4001	0.09999999999999998,-0.1} {4099	0.3,-0.5} {5020	-0.1,NaN,0.0} {4001	NaN,NaN} {4099	NaN,NaN}}

test select {list vabs} {
	split [exec cg select -f {begin {diff=vabs(vdef($freq-sample1,0) @- vdef($freq-sample2,0))}} < data/vars3.sft] \n
} {{begin	diff} {4000	0.09999999999999998} {4001	0.09999999999999998,0.1} {4099	0.3,0.5} {5020	0.1,0.5,0.0} {4001	0.9,0.8} {4099	0.2,0.4}}

test select {list vabs NaN} {
	split [exec cg select -f {begin {diff=vabs($freq-sample1 @- $freq-sample2)}} < data/vars3.sft] \n
} {{begin	diff} {4000	0.09999999999999998} {4001	0.09999999999999998,0.1} {4099	0.3,0.5} {5020	0.1,NaN,0.0} {4001	NaN,NaN} {4099	NaN,NaN}}

test select {list @+} {
	split [exec cg select -f {begin {diff=vabs(vdef($freq-sample1,0) @+ vdef($freq-sample2,0))}} < data/vars3.sft] \n
} {{begin	diff} {4000	0.9} {4001	0.9,0.30000000000000004} {4099	0.7,0.7} {5020	0.30000000000000004,0.5,0.4} {4001	0.9,0.8} {4099	0.2,0.4}}

test select {list @*} {
	split [exec cg select -f {begin {diff=vabs(vdef($freq-sample1,0) @* vdef($freq-sample2,0))}} < data/vars3.sft] \n
} {{begin	diff} {4000	0.2} {4001	0.2,0.020000000000000004} {4099	0.1,0.06} {5020	0.020000000000000004,0.0,0.04000000000000001} {4001	0.0,0.0} {4099	0.0,0.0}}

test select {list @/} {
	split [exec cg select -f {begin {diff=vabs(vdef($freq-sample1,0) @/ vdef($freq-sample2,0))}} < data/vars3.sft] \n
} {{begin	diff} {4000	1.25} {4001	1.25,0.5} {4099	2.5,0.16666666666666669} {5020	0.5,Inf,1.0} {4001	0.0,Inf} {4099	Inf,Inf}}

test tokenize {@newop precedence} {
	set code {$a @- 1 - $b * $b / 2 + 3}
	set header {a b}; set a 1; set b 1
	set tokens [tsv_select_tokenize $header $code n]
	set pquery [tsv_select_detokenize $tokens $header n]
} {((vminus(${a}, 1) - ((${b} * ${b}) / 2)) + 3)}

test tokenize {@newop precedence} {
	set code {$a @- 1 - $b * $b / 2 @+ 3}
	set header {a b}; set a 1; set b 1
	set tokens [tsv_select_tokenize $header $code n]
	set pquery [tsv_select_detokenize $tokens $header n]
} {vplus((vminus(${a}, 1) - ((${b} * ${b}) / 2)), 3)}

test select {lindex} {
	split [exec cg select \
		-f {chromosome begin} \
		-q {lindex($freq-sample1,1) == 0.1} \
		< data/vars3.sft] \n
} {{chromosome	begin} {chr1	4001} {chr1	4099}}

test select {and} {
	split [exec cg select \
		-f {chromosome begin} \
		-q {$freq-sample1 contains 0.5 and $freq-sample2 contains 0.2} \
		< data/vars3.sft] \n
} {{chromosome	begin} {chr1	4001} {chr1	4099} {chr1	5020}}

test select {or} {
	split [exec cg select \
		-f {chromosome begin} \
		-q {$freq-sample1 contains 0.5 or $freq-sample2 contains 0.2} \
		< data/vars3.sft] \n
} {{chromosome	begin} {chr1	4000} {chr1	4001} {chr1	4099} {chr1	5020}}

test select {vand} {
	split [exec cg select \
		-f {chromosome begin} \
		-q {lone($freq-sample1 @== 0.5 vand $freq-sample2 @== 0.2)} \
		< data/vars3.sft] \n
} {{chromosome	begin} {chr1	4099}}

test select {@&&} {
	split [exec cg select \
		-f {chromosome begin} \
		-q {lone($freq-sample1 @== 0.5 @&& $freq-sample2 @== 0.2)} \
		< data/vars3.sft] \n
} {{chromosome	begin} {chr1	4099}}

test select {vor} {
	split [exec cg select \
		-f {chromosome begin} \
		-q {lone($freq-sample1 @== 0.5 vor $freq-sample2 @== 0.2)} \
		< data/vars3.sft] \n
} {{chromosome	begin} {chr1	4000} {chr1	4001} {chr1	4099} {chr1	5020}}

test select {@||} {
	split [exec cg select \
		-f {chromosome begin} \
		-q {lone($freq-sample1 @== 0.5 @|| $freq-sample2 @== 0.2)} \
		< data/vars3.sft] \n
} {{chromosome	begin} {chr1	4000} {chr1	4001} {chr1	4099} {chr1	5020}}

test select {region} {
	split [exec cg select -f {chromosome begin} -q {region("chr2:4099-5020")} < data/vars1.sft] \n
} {{chromosome	begin} {chr2	4099} {chr2	5000} {chr2	5000} {chr2	5005} {chr2	5010} {chr2	5011}}

test select {region, only chr} {
	split [exec cg select -f {chromosome begin} -q {region("chr1")} < data/vars1.sft] \n
} {{chromosome	begin} {chr1	4000} {chr1	4001} {chr1	4099} {chr1	5000} {chr1	5020}}


test select {empty first field bug check} {
	cg select -f test data/select-bugcheck.tsv
} {test
0
100}

test select {-q error on non number <} {
	exec cg select -q {$mixed < 4} data/table.tsv
} {a4 is not a number} regexp error

test select {error missing quote} {
	exec cg select -q {$regtest != "aa} -f {chromosome begin end type ref alt alleleSeq1-sample1 alleleSeq2-sample1 coverage-sample1 sequenced-sample1 alleleSeq1-sample2 alleleSeq2-sample2 coverage-sample2 sequenced-sample2} ../tests/data/expected-vars1-reg_annot.sft tmp/tempexpected.tsv
} {error: incomplete quoted expression: "aa} error

test select {error missing quote empty} {
	exec cg select -q {$regtest != "} -f {chromosome begin end type ref alt alleleSeq1-sample1 alleleSeq2-sample1 coverage-sample1 sequenced-sample1 alleleSeq1-sample2 alleleSeq2-sample2 coverage-sample2 sequenced-sample2} ../tests/data/expected-vars1-reg_annot.sft tmp/tempexpected.tsv
} {error: incomplete quoted expression: "} error

test select {brokentable} {
	exec cg select -q {$other == "cc"} data/brokentable.tsv
} {wrong number of fields for line} regexp error

test select {-q use calculated column from -f} {
	exec cg select -f {{calc=$num+1} *} -q {$calc <= 4} data/table.tsv
} {calc	num	text	mixed	other
3	2	c	a2	cc}

test select {-f use calculated column in other calculated column} {
	exec cg select -f {{calc=$num+1} {calc2=$calc+1} *} -q {$calc <= 4} data/table.tsv
} {calc	calc2	num	text	mixed	other
3	4	2	c	a2	cc}

test select {-q use calculated column from -f without it being in the output (using -)} {
	exec cg select -f {{-calc=$num+1} *} -q {$calc <= 4} data/table.tsv
} {num	text	mixed	other
2	c	a2	cc}

test select {calculated column with wildcard} {
	exec cg select -f {test-*="$alleleSeq1-*/$alleleSeq2-*"} -q {$ROW <= 4} data/expected-multicompar_reannot-var_annotvar_annot2.sft
} {test-annot1	test-annot2
T/C	T/C
A/C	A/C
T/T	-/-
AGCGTGGCAA/	AGCGTGGCAA/
A/A	G/C}

test select {calculated column with wildcard also outside var} {
	exec cg select -f {{test-*="*: $alleleSeq1-*/$alleleSeq2-*"}} -q {$ROW < 2} data/expected-multicompar_reannot-var_annotvar_annot2.sft
} {test-annot1	test-annot2
annot1: T/C	annot2: T/C
annot1: A/C	annot2: A/C}

test select {calculated column with multiple wildcards} {
	exec cg select -f {a*-**="$alleleSeq*-**"} -q {$ROW < 2} data/expected-multicompar_reannot-var_annotvar_annot2.sft
} {a1-annot1	a1-annot2	a2-annot1	a2-annot2
T	T	C	C
A	A	C	C}

test select {sampledata in fields} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	freq-sample1	freq-sample2	freq-sample3
		1	0.4	0.8	1.0
	}
	write_tab tmp/temp.sampledata.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select -f {id gender-sample1 gender-sample2 gender-sample3} tmp/temp.tsv
} {id	gender-sample1	gender-sample2	gender-sample3
1	m	f	f}

test select {sampledata in fields, missing sample} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	freq-sample1	freq-sample2	freq-sample3
		1	0.4	0.8	1.0
	}
	write_tab tmp/temp.sampledata.tsv {
		id	gender
		sample1	m
		sample2	f
	}
	exec cg select -f {id gender-sample1 gender-sample2 gender-sample3} tmp/temp.tsv
} {field "gender-sample3" not present, also not in sampledata file tmp/temp.sampledata.tsv} error

test select {sampledata in fields, empty} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	freq-sample1	freq-sample2	freq-sample3
		1	0.4	0.8	1.0
	}
	write_tab tmp/temp.sampledata.tsv {
		id	gender
		sample1	m
		sample2	{}
		sample3	f
	}
	exec cg select -f {id gender-sample1 gender-sample2 gender-sample3} tmp/temp.tsv
} {id	gender-sample1	gender-sample2	gender-sample3
1	m		f}

test select {sampledata in fields using -sd} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	freq-sample1	freq-sample2	freq-sample3
		1	0.4	0.8	1.0
	}
	write_tab tmp/sampledata.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select -sd tmp/sampledata.tsv -f {id gender-sample1 gender-sample2 gender-sample3} tmp/temp.tsv
} {id	gender-sample1	gender-sample2	gender-sample3
1	m	f	f}

test select {sampledata in code of calculated column} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	freq-sample1	freq-sample2	freq-sample3
		1	0.4	0.8	1.0
		2	0.8	0.9	0.3
	}
	write_tab tmp/temp.sampledata.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select -f {id {count=scount($freq > 0.5 and $gender eq "f")}} tmp/temp.tsv
} {id	count
1	2
2	1}

test select {sampledata in code of query} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	freq-sample1	freq-sample2	freq-sample3
		1	0.4	0.8	1.0
		2	0.8	0.9	0.3
	}
	write_tab tmp/temp.sampledata.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select -f id -q {scount($freq > 0.5 and $gender eq "f") > 1} tmp/temp.tsv
} {id
1}

test select {sampledata in calc field wildcard} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	freq-sample1	freq-sample2	freq-sample3
		1	0.4	0.8	1.0
	}
	write_tab tmp/temp.sampledata.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select -f {id {g-*=$gender-*}} tmp/temp.tsv
} {id	g-sample1	g-sample2	g-sample3
1	m	f	f}

testsummarize
