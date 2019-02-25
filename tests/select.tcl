#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

set dbopt "-db test"
set dbopt ""
set dboptt ""

proc selecttests {} {
global dbopt
set dboptt " $::dbopt"

test select "-f$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -f "num text" data/table.tsv
} {num	text
4	a
10	b
2	c
100	d}

test select "-rf$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -rf "mixed other" data/table.tsv
} {num	text
4	a
10	b
2	c
100	d}

test select "-q$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -q {$num <= 4} data/table.tsv
} {num	text	mixed	other
4	a	a4	aaaa
2	c	a2	cc}

test select "-s$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -s num data/table.tsv
} {num	text	mixed	other
2	c	a2	cc
4	a	a4	aaaa
10	b	b2	bbbbbbbbbb
100	d	aa3	dd..}

test select "-s$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -sr num data/table.tsv
} {num	text	mixed	other
100	d	aa3	dd..
10	b	b2	bbbbbbbbbb
4	a	a4	aaaa
2	c	a2	cc}

test select "-s -f$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -s num -f "num mixed" data/table.tsv
} {num	mixed
2	a2
4	a4
10	b2
100	aa3}

test select "-s -q$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -s num -q {$num <= 4} data/table.tsv
} {num	text	mixed	other
2	c	a2	cc
4	a	a4	aaaa}

test select "-q -f$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -f "num mixed" -q {$num <= 4} data/table.tsv
} {num	mixed
4	a4
2	a2}

test select "-s -f -q$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -s num  -f "num mixed" -q {$num <= 4} data/table.tsv
} {num	mixed
2	a2
4	a4}

test select "-q multiple lines, tabs$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -s num  -f "num mixed" -q {
		$num <= 4
		&& $text == "c"
	} data/table.tsv
} {num	mixed
2	a2}

test select "default -s -$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -s - data/expected-test1000glow.vcf2tsv tmp/result.tsv
	exec diff data/expected-test1000glow.vcf2tsv tmp/result.tsv
} {}

test select "-f *$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -f {chromosome begin end alleleSeq1-*} -q {$begin == 4000} [gzfile data/vars1.sft]
} {chromosome	begin	end	alleleSeq1-sample1	alleleSeq1-sample2
chr1	4000	4001	A	A
chr2	4000	4001	G	G}

test select "-f calculated$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -f {chromosome begin end {geno1="$alleleSeq1-sample1/$alleleSeq2-sample1"}} -q {$begin == 4000} [gzfile data/vars1.sft]
} {chromosome	begin	end	geno1
chr1	4000	4001	A/G
chr2	4000	4001	G/G}

test select "-f calculated functions$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -f {chromosome begin end {countG=count($alleleSeq*, == "G")}} -q {$begin == 4000} [gzfile data/vars1.sft]
} {chromosome	begin	end	countG
chr1	4000	4001	2
chr2	4000	4001	4}

test select "-f calculated functions error unknown function$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -f {chromosome begin end {error=blabla($alt)}} [gzfile data/vars1.sft] tmp/temp.tsv
} {unknown function blabla} error

test select "-f calculated functions error catch$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -f {chromosome begin end {error=catch(blabla($alt),"error")}} -q {$ROW < 2} [gzfile data/vars1.sft]
} {chromosome	begin	end	error
chr1	4000	4001	error
chr1	4001	4002	error}

test select "-f calculated functions error not a number$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -f {chromosome begin end {error=($list<1)}} [gzfile data/vars1.sft] tmp/temp.tsv
} {1;2,3;4 is not a number} error

test select "-f wildcard error$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -f {chromosome begin end {countG=count($alleleSeq-*, == "G")}} -q {$begin == 4000} [gzfile data/vars1.sft]
} {field alleleSeq-* not present in file (or sampleinfo)} error

test select "-f calculated if$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -f {chromosome begin end {countG=if(($alleleSeq1-sample1 == "A" || $alleleSeq2-sample1 == "A"),"hasA","noA")}} -q {$begin == 4000} [gzfile data/vars1.sft]
} {chromosome	begin	end	countG
chr1	4000	4001	hasA
chr2	4000	4001	noA}

test select "-f calculated if$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -f {chromosome begin end {countG=if(count($alleleSeq*, == "G")<4,"<4",">=4")}} -q {$begin == 4000} [gzfile data/vars1.sft]
} {chromosome	begin	end	countG
chr1	4000	4001	<4
chr2	4000	4001	>=4}

test select "-hc 1$dboptt" {
	global dbopt
	file_write tmp/temp.tsv [deindent {
		#a	b
		1	2
	}]
	file delete tmp/result.tsv
	cg select {*}$dbopt -h tmp/temp.tsv >> tmp/result.tsv
	cg select {*}$dbopt -hc 1 -h tmp/temp.tsv >> tmp/result.tsv
	cg select {*}$dbopt -hc 1 -f a tmp/temp.tsv >> tmp/result.tsv
	file_read tmp/result.tsv
} {1
2
a
b
a
1
}

test select "keep header info and format rtg: -hc$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -hc 1 -s position data/rtgsnps.tsv tmp/temp.tsv
	file delete tmp/temp
	catch {exec diff tmp/temp.tsv data/rtgsnps.tsv > tmp/temp}
	file_read tmp/temp
} {3c3
< name	position	type	reference	prediction	posterior	coverage	correction	support_statistics
---
> #name	position	type	reference	prediction	posterior	coverage	correction	support_statistics
}

test select "keep header info and format rtg: -hc$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -hc 2 -s position data/rtgsnps.tsv tmp/temp.tsv
	file delete tmp/temp
	catch {exec diff tmp/temp.tsv data/rtgsnps.tsv > tmp/temp}
	file_read tmp/temp
} {}

test select "-hc with chars that must be escaped$dboptt" {
	global dbopt
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
	exec cg select {*}$dbopt -hc 1 -f {{b c}} tmp/temp.tsv tmp/out.tsv
	exec diff tmp/out.tsv tmp/expected.tsv
} {}

test select "-hf$dboptt" {
	global dbopt
	write_tab tmp/temp_header.tsv {
		key	value
	}
	write_tab tmp/temp.tsv {
		a	1
		b	2
	}
	exec cg select {*}$dbopt -hf tmp/temp_header.tsv -f {value} tmp/temp.tsv
} {value
1
2}

test select "-hp$dboptt" {
	global dbopt
	write_tab tmp/temp.tsv {
		a	1
		b	2
	}
	exec cg select {*}$dbopt -hp "key value" -f {value} tmp/temp.tsv
} {value
1
2}

test select "-hp with comments$dboptt" {
	global dbopt
	write_tab tmp/temp.tsv {
		#test	test2
		a	1
		b	2
	}
	exec cg select {*}$dbopt -hp "key value" -f {value} tmp/temp.tsv
} {value
1
2}

test select "keep header info and format vcf$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -s POS data/test.vcf tmp/temp.tsv
	file delete tmp/temp
	catch {exec diff tmp/temp.tsv data/test.vcf > tmp/temp}
	file_read tmp/temp
} {}

# cannot sort on calculated fields (yet)
if 0 {
test select "-f calculated functions + sort$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -f {chromosome begin end {countG=count($alleleSeq*, == "G")}} -q {$begin == 4000} -s countG [gzfile data/vars1.sft]
} {chromosome	begin	end	countG
chr1	4000	4001	2
chr2	4000	4001	4}
}

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

test select "sm$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt -f {chromosome begin} -q {sm(sample1,sample2)} < data/vars1.sft] \n
} {{chromosome	begin} {chr1	4000} {chr1	4001} {chr1	4099} {chr2	4000} {chr2	4001} {chr2	4099} {chr2	5010}}

test select "df$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt -f {chromosome begin} -q {df(sample1,sample2)} < data/vars1.sft] \n
} {{chromosome	begin} {chr1	5020} {chr2	5000} {chr2	5011} {chr2	8000}}

test select "mm$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt -f {chromosome begin} -q {mm(sample1,sample2)} < data/vars1.sft] \n
} {{chromosome	begin} {chr2	5005}}

test select "-q count with wildcard$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -f {chromosome begin end alleleSeq*} -q {count($alleleSeq*, == "G") > 3} [gzfile data/vars1.sft]
} {chromosome	begin	end	alleleSeq1-sample1	alleleSeq2-sample1	alleleSeq1-sample2	alleleSeq2-sample2
chr1	4001	4002	G	G	G	G
chr2	4000	4001	G	G	G	G
chr2	4001	4002	G	G	G	G}

test select "lmin column$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt -f {chromosome begin {lmin=lmin($list)}} < data/vars1.sft] \n
} {{chromosome	begin	lmin} {chr1	4000	4} {chr1	4001	1} {chr1	4099	1} {chr1	5000	1} {chr1	5020	3} {chr2	4000	2} {chr2	4001	2} {chr2	4099	2} {chr2	5000	3} {chr2	5000	2} {chr2	5005	4} {chr2	5010	20} {chr2	5011	NaN} {chr2	8000	NaN}}

test select "lmind column$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt -f {chromosome begin {lmin=lmind($list,10)}} < data/vars1.sft] \n
} {{chromosome	begin	lmin} {chr1	4000	4} {chr1	4001	1} {chr1	4099	1} {chr1	5000	1} {chr1	5020	3} {chr2	4000	2} {chr2	4001	2} {chr2	4099	2} {chr2	5000	3} {chr2	5000	2} {chr2	5005	4} {chr2	5010	10} {chr2	5011	10} {chr2	8000	10}}

test select "lmind select$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt -f {chromosome begin} -q {lmind($list,10) == 2} < data/vars1.sft] \n
} {{chromosome	begin} {chr2	4000} {chr2	4001} {chr2	4099} {chr2	5000}}

test select "lmax column$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt -f {chromosome begin {lmax=lmax($list)}} < data/vars1.sft] \n
} {{chromosome	begin	lmax} {chr1	4000	4} {chr1	4001	4} {chr1	4099	2} {chr1	5000	2} {chr1	5020	3} {chr2	4000	2} {chr2	4001	4} {chr2	4099	4} {chr2	5000	3} {chr2	5000	2} {chr2	5005	4} {chr2	5010	20} {chr2	5011	NaN} {chr2	8000	NaN}}

test select "lmaxd column$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt -f {chromosome begin {lmax=lmaxd($list,0)}} < data/vars1.sft] \n
} {{chromosome	begin	lmax} {chr1	4000	4} {chr1	4001	4} {chr1	4099	2} {chr1	5000	2} {chr1	5020	3} {chr2	4000	2} {chr2	4001	4} {chr2	4099	4} {chr2	5000	3} {chr2	5000	2} {chr2	5005	4} {chr2	5010	20} {chr2	5011	0} {chr2	8000	0}}

test select "lmaxd select$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt -f {chromosome begin} -q {lmaxd($list,10) == 2} < data/vars1.sft] \n
} {{chromosome	begin} {chr1	4099} {chr1	5000} {chr2	4000} {chr2	5000}}

test select "counthasone column$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt -f {chromosome begin {lmin=counthasone($list, ==2)}} < data/vars1.sft] \n
} {{chromosome	begin	lmin} {chr1	4000	0} {chr1	4001	1} {chr1	4099	1} {chr1	5000	1} {chr1	5020	0} {chr2	4000	1} {chr2	4001	1} {chr2	4099	1} {chr2	5000	0} {chr2	5000	1} {chr2	5005	0} {chr2	5010	0} {chr2	5011	0} {chr2	8000	0}}

test select "ROW$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt -f {chromosome begin end} -q {$ROW == 2} < data/vars1.sft] \n
} {{chromosome	begin	end} {chr1	4099	4100}}

test select "ROW$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt -f {chromosome begin end ROW} -q {$ROW >= 2 && $ROW <= 3} < data/vars1.sft] \n
} {{chromosome	begin	end	ROW} {chr1	4099	4100	2} {chr1	5000	5010	3}}

test select "shared objects bugcheck$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt -f {chromosome begin end {test=[set ::keep $begin]}} -q {$ROW between {2 3}} < data/vars1.sft] \n
} {{chromosome	begin	end	test} {chr1	4099	4100	4099} {chr1	5000	5010	5000}}

test select "shared objects bugcheck$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt -f {chromosome begin end {test="[get ::keep 1]-[set ::keep $begin]"}} -q {$ROW between {2 3}} < data/vars1.sft] \n
} {{chromosome	begin	end	test} {chr1	4099	4100	1-4099} {chr1	5000	5010	4099-5000}}

test select "start brace bugcheck$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt -f {chromosome begin end {type=($type == "snp")? "Snp" : (($type == "del")? "Deletion" : $type)}} -q {$ROW in {2 3 8 10}} < data/vars1.sft] \n
} {{chromosome	begin	end	type} {chr1	4099	4100	Snp} {chr1	5000	5010	Deletion} {chr2	5000	5000	ins} {chr2	5005	5006	Snp}}

test select "list @-$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt -f {begin {diff=vdef($freq-sample1,0) @- vdef($freq-sample2,0)}} < data/vars3.sft] \n
} {{begin	diff} {4000	0.09999999999999998} {4001	0.09999999999999998,-0.1} {4099	0.3,-0.5} {5020	-0.1,0.5,0.0} {4001	-0.9,0.8} {4099	0.2,0.4}}

test select "list @- @-$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt -f {begin {diff=vdef($freq-sample1,0) @- vdef($freq-sample2,0) @- vdef($freq-sample2,0)}} < data/vars3.sft] \n
} {{begin	diff} {4000	-0.30000000000000004} {4001	-0.30000000000000004,-0.30000000000000004} {4099	0.09999999999999998,-1.1} {5020	-0.30000000000000004,0.5,-0.2} {4001	-1.8,0.8} {4099	0.2,0.4}}

test select "list @- @<$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt \
		-f {begin {diff=vdef($freq-sample1,0) @- vdef($freq-sample2,0)}} \
		-q {lone((vdef($freq-sample1,0) @- vdef($freq-sample2,0)) @< 0)} \
		< data/vars3.sft] \n
} {{begin	diff} {4001	0.09999999999999998,-0.1} {4099	0.3,-0.5} {5020	-0.1,0.5,0.0} {4001	-0.9,0.8}}

test select "list @- NaN$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt -f {begin {diff=$freq-sample1 @- $freq-sample2}} < data/vars3.sft] \n
} {{begin	diff} {4000	0.09999999999999998} {4001	0.09999999999999998,-0.1} {4099	0.3,-0.5} {5020	-0.1,NaN,0.0} {4001	NaN,NaN} {4099	NaN,NaN}}

test select "list vabs$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt -f {begin {diff=vabs(vdef($freq-sample1,0) @- vdef($freq-sample2,0))}} < data/vars3.sft] \n
} {{begin	diff} {4000	0.09999999999999998} {4001	0.09999999999999998,0.1} {4099	0.3,0.5} {5020	0.1,0.5,0.0} {4001	0.9,0.8} {4099	0.2,0.4}}

test select "list vabs NaN$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt -f {begin {diff=vabs($freq-sample1 @- $freq-sample2)}} < data/vars3.sft] \n
} {{begin	diff} {4000	0.09999999999999998} {4001	0.09999999999999998,0.1} {4099	0.3,0.5} {5020	0.1,NaN,0.0} {4001	NaN,NaN} {4099	NaN,NaN}}

test select "list @+$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt -f {begin {diff=vabs(vdef($freq-sample1,0) @+ vdef($freq-sample2,0))}} < data/vars3.sft] \n
} {{begin	diff} {4000	0.9} {4001	0.9,0.30000000000000004} {4099	0.7,0.7} {5020	0.30000000000000004,0.5,0.4} {4001	0.9,0.8} {4099	0.2,0.4}}

test select "list @*$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt -f {begin {diff=vabs(vdef($freq-sample1,0) @* vdef($freq-sample2,0))}} < data/vars3.sft] \n
} {{begin	diff} {4000	0.2} {4001	0.2,0.020000000000000004} {4099	0.1,0.06} {5020	0.020000000000000004,0.0,0.04000000000000001} {4001	0.0,0.0} {4099	0.0,0.0}}

test select "list @/$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt -f {begin {diff=vabs(vdef($freq-sample1,0) @/ vdef($freq-sample2,0))}} < data/vars3.sft] \n
} {{begin	diff} {4000	1.25} {4001	1.25,0.5} {4099	2.5,0.16666666666666669} {5020	0.5,Inf,1.0} {4001	0.0,Inf} {4099	Inf,Inf}}

test tokenize {@newop precedence} {
	set code {$a @- 1 - $b * $b / 2 + 3}
	set header {a b}; set a 1; set b 1
	set tokens [tsv_select_tokenize $header $code n]
	set pquery [tsv_select_detokenize $tokens $header n]
} {((vminus(${a}, 1) - ((${b} * ${b}) *1.0/ 2)) + 3)}

test tokenize {@newop precedence} {
	set code {$a @- 1 - $b * $b / 2 @+ 3}
	set header {a b}; set a 1; set b 1
	set tokens [tsv_select_tokenize $header $code n]
	set pquery [tsv_select_detokenize $tokens $header n]
} {vplus((vminus(${a}, 1) - ((${b} * ${b}) *1.0/ 2)), 3)}

test select "lindex$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt \
		-f {chromosome begin} \
		-q {lindex($freq-sample1,1) == 0.1} \
		< data/vars3.sft] \n
} {{chromosome	begin} {chr1	4001} {chr1	4099}}

test select "and$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt \
		-f {chromosome begin} \
		-q {$freq-sample1 contains 0.5 and $freq-sample2 contains 0.2} \
		< data/vars3.sft] \n
} {{chromosome	begin} {chr1	4001} {chr1	4099} {chr1	5020}}

test select "or$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt \
		-f {chromosome begin} \
		-q {$freq-sample1 contains 0.5 or $freq-sample2 contains 0.2} \
		< data/vars3.sft] \n
} {{chromosome	begin} {chr1	4000} {chr1	4001} {chr1	4099} {chr1	5020}}

test select "vand$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt \
		-f {chromosome begin} \
		-q {lone($freq-sample1 @== 0.5 vand $freq-sample2 @== 0.2)} \
		< data/vars3.sft] \n
} {{chromosome	begin} {chr1	4099}}

test select "@&&$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt \
		-f {chromosome begin} \
		-q {lone($freq-sample1 @== 0.5 @&& $freq-sample2 @== 0.2)} \
		< data/vars3.sft] \n
} {{chromosome	begin} {chr1	4099}}

test select "vor$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt \
		-f {chromosome begin} \
		-q {lone($freq-sample1 @== 0.5 vor $freq-sample2 @== 0.2)} \
		< data/vars3.sft] \n
} {{chromosome	begin} {chr1	4000} {chr1	4001} {chr1	4099} {chr1	5020}}

test select "@||$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt \
		-f {chromosome begin} \
		-q {lone($freq-sample1 @== 0.5 @|| $freq-sample2 @== 0.2)} \
		< data/vars3.sft] \n
} {{chromosome	begin} {chr1	4000} {chr1	4001} {chr1	4099} {chr1	5020}}

test select "region$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt -f {chromosome begin} -q {region("chr2:4099-5020")} < data/vars1.sft] \n
} {{chromosome	begin} {chr2	4099} {chr2	5000} {chr2	5000} {chr2	5005} {chr2	5010} {chr2	5011}}

test select "region, only chr$dboptt" {
	global dbopt
	split [exec cg select {*}$dbopt -f {chromosome begin} -q {region("chr1")} < data/vars1.sft] \n
} {{chromosome	begin} {chr1	4000} {chr1	4001} {chr1	4099} {chr1	5000} {chr1	5020}}


test select "empty first field bug check$dboptt" {
	global dbopt
	cg select {*}$dbopt -f test data/select-bugcheck.tsv
} {test
0
100}

test select "-q error on non number <$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -q {$mixed < 4} data/table.tsv
} {a4 is not a number} regexp error

test select "error missing quote$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -q {$regtest != "aa} -f {chromosome begin end type ref alt alleleSeq1-sample1 alleleSeq2-sample1 coverage-sample1 sequenced-sample1 alleleSeq1-sample2 alleleSeq2-sample2 coverage-sample2 sequenced-sample2} ../tests/data/expected-vars1-reg_annot.sft tmp/tempexpected.tsv
} {error: incomplete quoted expression: "aa} error

test select "error missing quote empty$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -q {$regtest != "} -f {chromosome begin end type ref alt alleleSeq1-sample1 alleleSeq2-sample1 coverage-sample1 sequenced-sample1 alleleSeq1-sample2 alleleSeq2-sample2 coverage-sample2 sequenced-sample2} ../tests/data/expected-vars1-reg_annot.sft tmp/tempexpected.tsv
} {error: incomplete quoted expression: "} error

test select "brokentable$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -q {$other == "cc"} data/brokentable.tsv
} {wrong number of fields for line} regexp error

test select "-q use calculated column from -f$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -f {{calc=$num+1} *} -q {$calc <= 4} data/table.tsv
} {calc	num	text	mixed	other
3	2	c	a2	cc}

test select "-f use calculated column in other calculated column$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -f {{calc=$num+1} {calc2=$calc+1} *} -q {$calc <= 4} data/table.tsv
} {calc	calc2	num	text	mixed	other
3	4	2	c	a2	cc}

test select "-f use calculated column in other calculated column, using different vars$dboptt" {
	global dbopt
	write_tab tmp/temp.tsv {
		id	num1	num2
		1	1	2
		2	3	4
	}
	exec cg select {*}$dbopt -f {id {calc=$num1+1} {calc2=$calc+$num2}} -q {$calc2 > 5} tmp/temp.tsv
} {id	calc	calc2
2	4	8}

test select "-q use calculated column from -f without it being in the output (using -)$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -f {{-calc=$num+1} *} -q {$calc <= 4} data/table.tsv
} {num	text	mixed	other
2	c	a2	cc}

test select "calculated column with wildcard$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -f {test-*="$alleleSeq1-*/$alleleSeq2-*"} -q {$ROW <= 4} data/expected-multicompar_reannot-var_annotvar_annot2.sft
} {test-annot1	test-annot2
T/C	T/C
A/C	A/C
T/T	?/?
AGCGTGGCAA/	AGCGTGGCAA/
A/A	G/C}

test select "calculated column with wildcard also outside var$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -f {{test-*="*: $alleleSeq1-*/$alleleSeq2-*"}} -q {$ROW < 2} data/expected-multicompar_reannot-var_annotvar_annot2.sft
} {test-annot1	test-annot2
annot1: T/C	annot2: T/C
annot1: A/C	annot2: A/C}

test select "calculated column with multiple wildcards$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -f {a*-**="$alleleSeq*-**"} -q {$ROW < 2} data/expected-multicompar_reannot-var_annotvar_annot2.sft
} {a1-annot1	a1-annot2	a2-annot1	a2-annot2
T	T	C	C
A	A	C	C}

test select "sampleinfo in fields$dboptt" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	freq-sample1	freq-sample2	freq-sample3
		1	0.4	0.8	1.0
	}
	write_tab tmp/temp.tsv.sampleinfo {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select {*}$dbopt -f {id gender-sample1 gender-sample2 gender-sample3} tmp/temp.tsv
} {id	gender-sample1	gender-sample2	gender-sample3
1	m	f	f}

test select "sampleinfo in fields with full analysis$dboptt" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	0.4	0.8	1.0
	}
	write_tab tmp/temp.tsv.sampleinfo {
		id	gender
		sample1	m
		sample2	f
	}
	exec cg select {*}$dbopt -f {id gender-gatk-sample1	gender-sam-sample1 gender-gatk-sample2} tmp/temp.tsv
} {id	gender-gatk-sample1	gender-sam-sample1	gender-gatk-sample2
1	m	m	f}

test select "sampleinfo in fields with full analysis but not in fields asked$dboptt" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	0.4	0.8	1.0
	}
	write_tab tmp/temp.tsv.sampleinfo {
		id	gender
		sample1	m
		sample2	f
	}
	exec cg select {*}$dbopt -f {id gender-sample1 gender-sample2} tmp/temp.tsv
} {id	gender-sample1	gender-sample2
1	m	f}

test select "sampleinfo (other filename convention) in fields$dboptt" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	freq-sample1	freq-sample2	freq-sample3
		1	0.4	0.8	1.0
	}
	write_tab tmp/temp.tsv.sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select {*}$dbopt -f {id gender-sample1 gender-sample2 gender-sample3} tmp/temp.tsv
} {id	gender-sample1	gender-sample2	gender-sample3
1	m	f	f}

test select "sampleinfo in fields, missing sample$dboptt" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	freq-sample1	freq-sample2	freq-sample3
		1	0.4	0.8	1.0
	}
	write_tab tmp/temp.tsv.sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	f
	}
	exec cg select {*}$dbopt -f {id gender-sample1 gender-sample2 gender-sample3} tmp/temp.tsv
} {field "gender-sample3" not present, also not in sampleinfo file tmp/temp.tsv.sampleinfo.tsv} error

test select "sampleinfo in fields, empty$dboptt" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	freq-sample1	freq-sample2	freq-sample3
		1	0.4	0.8	1.0
	}
	write_tab tmp/temp.tsv.sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	{}
		sample3	f
	}
	exec cg select {*}$dbopt -f {id gender-sample1 gender-sample2 gender-sample3} tmp/temp.tsv
} {id	gender-sample1	gender-sample2	gender-sample3
1	m		f}

test select "sampleinfo in fields using -si$dboptt" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	freq-sample1	freq-sample2	freq-sample3
		1	0.4	0.8	1.0
	}
	write_tab tmp/sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select {*}$dbopt -si tmp/sampleinfo.tsv -f {id gender-sample1 gender-sample2 gender-sample3} tmp/temp.tsv
} {id	gender-sample1	gender-sample2	gender-sample3
1	m	f	f}

test select "sampleinfo in code of calculated column$dboptt" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	freq-sample1	freq-sample2	freq-sample3
		1	0.4	0.8	1.0
		2	0.8	0.9	0.3
	}
	write_tab tmp/temp.tsv.sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select {*}$dbopt -f {id {count=scount($freq > 0.5 and $gender eq "f")}} tmp/temp.tsv
} {id	count
1	2
2	1}

test select "sampleinfo in code of query$dboptt" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	freq-sample1	freq-sample2	freq-sample3
		1	0.4	0.8	1.0
		2	0.8	0.9	0.3
	}
	write_tab tmp/temp.tsv.sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select {*}$dbopt -f id -q {scount($freq > 0.5 and $gender eq "f") > 1} tmp/temp.tsv
} {id
1}

test select "sampleinfo in calc field wildcard$dboptt" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	freq-sample1	freq-sample2	freq-sample3
		1	0.4	0.8	1.0
	}
	write_tab tmp/temp.tsv.sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select {*}$dbopt -f {id {g-*=$gender-*}} tmp/temp.tsv
} {id	g-sample1	g-sample2	g-sample3
1	m	f	f}

test select "sampleinfo in in -f and -q saggregate$dboptt" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	freq-sample1	freq-sample2	freq-sample3
		1	0.4	0.8	1.0
	}
	write_tab tmp/temp.tsv.sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select {*}$dbopt -q {scount($gender eq "f" and $freq > 0.9) > 0} -f {id gender-sample1 gender-sample2 gender-sample3} tmp/temp.tsv
} {id	gender-sample1	gender-sample2	gender-sample3
1	m	f	f}

test select "sampleinfo ignoring prefix in in -f and -q saggregate$dboptt" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	freq-gatk-crsbwa-sample1	freq-gatk-crsbwa-sample2	freq-gatk-crsbwa-sample3
		1	0.4	0.8	1.0
	}
	write_tab tmp/temp.tsv.sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select {*}$dbopt -q {acount($gender eq "f" and $freq > 0.9) > 0} -f {id gender-gatk-crsbwa-sample1 gender-gatk-crsbwa-sample2 gender-gatk-crsbwa-sample3} tmp/temp.tsv
} {id	gender-gatk-crsbwa-sample1	gender-gatk-crsbwa-sample2	gender-gatk-crsbwa-sample3
1	m	f	f}

test select "wildcard calc column in query$dboptt" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	val-sample1	val-sample2	val-sample3
		1	2	0	1
		2	5	0	2
		3	0	1	3
		4	1	9	4
	}
	exec cg select {*}$dbopt -f {id {-valb-*=if($val-* == 0,0,1)}} -q {$valb-sample2} tmp/temp.tsv
} {id
3
4}

test select "do not make index on select only$dboptt" {
	test_cleantmp
	global dbopt
	file copy data/table.tsv tmp/table.tsv
	exec cg select {*}$dbopt -f "num text" tmp/table.tsv
	file exists tmp/table.tsv.index
} 0

test select "catch$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -f {chromosome begin end {countG=catch(someerror)}} -q {$begin == 4000} [gzfile data/vars1.sft]
} {chromosome	begin	end	countG
chr1	4000	4001	1
chr2	4000	4001	1}

test select "catch$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -f {chromosome begin end {countG=catch($begin)}} -q {$begin == 4000} [gzfile data/vars1.sft]
} {chromosome	begin	end	countG
chr1	4000	4001	0
chr2	4000	4001	0}

test select "catch with errorvalue $dboptt" {
	global dbopt
	exec cg select {*}$dbopt -f {chromosome begin end {countG=catch(someerror,"e")}} -q {$begin == 4000} [gzfile data/vars1.sft]
} {chromosome	begin	end	countG
chr1	4000	4001	e
chr2	4000	4001	e}

test select "catch with errorvalue $dboptt" {
	global dbopt
	exec cg select {*}$dbopt -f {chromosome begin end {countG=catch($begin,"e")}} -q {$begin == 4000} [gzfile data/vars1.sft]
} {chromosome	begin	end	countG
chr1	4000	4001	4000
chr2	4000	4001	4000}

test select "chr_clip -$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -q {chr_clip($chromosome) in {1 X Y}} data/allchr.tsv tmp/result.tsv
	file_read tmp/result.tsv
} {chromosome	pos	coverage
chr1	1	1
chrX	24	24
chrY	25	25
}

test select "toupper$dboptt" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	va
		1	a
		2	A
	}
	exec cg select {*}$dbopt -f {id {val=toupper($va)}} tmp/temp.tsv
} {id	val
1	A
2	A}

test select "long format with sampleinfo $dboptt" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		sample id	freq
		sample1	1	0.4
		sample2	1	0.8
		sample3	1	1.0
	}
	write_tab tmp/temp.tsv.sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select {*}$dbopt -q {$gender eq "f" and $freq > 0.9} -f {sample freq gender} tmp/temp.tsv
} {sample	freq	gender
sample3	1.0	f}

test select "long format with sampleinfo using analysis$dboptt" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		sample id	freq
		gatk-sample1	1	0.4
		gatk-sample2	1	0.8
		gatk-sample3	1	1.0
	}
	write_tab tmp/temp.tsv.sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select {*}$dbopt -q {$gender eq "f" and $freq > 0.9} -f {sample freq gender} tmp/temp.tsv
} {sample	freq	gender
gatk-sample3	1.0	f}

test select "long with * in -f" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		sample id	freq
		sample1	1	0.4
	}
	exec cg select {*}$dbopt -f {*} tmp/temp.tsv
} {sample	id	freq
sample1	1	0.4}


test select "long with sampleinfo in -f with *" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		sample id	freq
		sample1	1	0.4
		sample2	1	0.8
	}
	write_tab tmp/temp.tsv.sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select {*}$dbopt -f {*} tmp/temp.tsv
} {sample	id	freq
sample1	1	0.4
sample2	1	0.8}

test select "long with sampleinfo in -f with wildcards" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		sample id	freq
		sample1	1	0.4
		sample2	1	0.8
	}
	write_tab tmp/temp.tsv.sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select {*}$dbopt -f {* gender} tmp/temp.tsv
} {sample	id	freq	gender
sample1	1	0.4	m
sample2	1	0.8	f}

test select "long with sampleinfo in -f with **" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		sample id	freq
		sample1	1	0.4
		sample2	1	0.8
	}
	write_tab tmp/temp.tsv.sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select {*}$dbopt -f {**} tmp/temp.tsv
} {sample	id	freq	gender
sample1	1	0.4	m
sample2	1	0.8	f}

test select "long format with sampleinfo in calc $dboptt" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		sample id	freq
		sample1	1	0.4
		sample2	1	0.8
		sample3	1	1.0
	}
	write_tab tmp/temp.tsv.sampleinfo.tsv {
		id	gender
		sample1	m
		sample2	f
		sample3	f
	}
	exec cg select -stack 1 {*}$dbopt -q {$gender eq "f" and $freq > 0.9} -f {{t=concat($sample,$gender)}} tmp/temp.tsv
} {t
sample3f}

test select "old CGI gene format$dboptt" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		#BUILD  1.4.2.7
		#GENERATED_AT   2009-Oct-06 12:03:35.919080
		#GENERATED_BY   callannotate
		#GENE_SET       /REF/HUMAN-F_04-REF/gene_annotations.csv
		#TYPE   GENE-ANNOTATION
		#VARIATIONS     /ASM/GS000000078-ASM/callrpt/annotation/annotation_Variations20-40_ccf100.csv
		#VAR_ANN_SET    /REF/HUMAN-F_04-REF/dbSNP.csv
		#VAR_ANN_TYPE   dbSNP
		#VERSION        0.2
		
		
		>chromosome  begin   end
		0       1       2
		3       4       5
	}
	write_tab tmp/expected.tsv {
		#BUILD  1.4.2.7
		#GENERATED_AT   2009-Oct-06 12:03:35.919080
		#GENERATED_BY   callannotate
		#GENE_SET       /REF/HUMAN-F_04-REF/gene_annotations.csv
		#TYPE   GENE-ANNOTATION
		#VARIATIONS     /ASM/GS000000078-ASM/callrpt/annotation/annotation_Variations20-40_ccf100.csv
		#VAR_ANN_SET    /REF/HUMAN-F_04-REF/dbSNP.csv
		#VAR_ANN_TYPE   dbSNP
		#VERSION        0.2
		#
		#
		chromosome  begin   end
		0       1       2
		3       4       5
	}
	exec cg select {*}$dbopt -s {chromosome begin end} tmp/temp.tsv tmp/sorted.tsv
	exec diff tmp/sorted.tsv tmp/expected.tsv
} {}

test select "check bug in -s$dboptt" {
	global dbopt
	write_tab tmp/testsort.tsv {
		name
		test-100
		test-1-1
	}
	exec cg select {*}$dbopt -s name tmp/testsort.tsv
} {name
test-100
test-1-1}

test select "-samples$dboptt" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	freq-sample1	v-sample1	freq-sample2	v-sample2	freq-sample3	v-sample3	annot
		1	0.4	s1	0.8	s2	1.0	s3	a
	}
	exec cg select {*}$dbopt -samples {sample2 sample3} tmp/temp.tsv
} {id	freq-sample2	v-sample2	freq-sample3	v-sample3	annot
1	0.8	s2	1.0	s3	a}

test select "-ssamples$dboptt" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	freq-sample1	v-sample1	freq-sample2	v-sample2	freq-sample3	v-sample3	annot
		1	0.4	s1	0.8	s2	1.0	s3	a
	}
	exec cg select {*}$dbopt -ssamples {sample3 sample2} tmp/temp.tsv
} {id	freq-sample3	v-sample3	freq-sample2	v-sample2	annot
1	1.0	s3	0.8	s2	a}

test select "-ssamples$dboptt" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		chromosome	begin	end	freq-sample1	v-sample1	freq-sample2	v-sample2	freq-sample3	v-sample3	annot type
		chr1	1	2	0.4	s1	0.8	s2	1.0	s3	a	snp
	}
	exec cg select {*}$dbopt -ssamples {sample3 sample2} tmp/temp.tsv
} {chromosome	begin	end	type	freq-sample3	v-sample3	freq-sample2	v-sample2	annot
chr1	1	2	snp	1.0	s3	0.8	s2	a}

test select "-f {-*-*}$dboptt" {
	global dbopt
	test_cleantmp
	write_tab tmp/temp.tsv {
		id	freq-sample1	v-sample1	freq-sample2	v-sample2	freq-sample3	v-sample3	annot
		1	0.4	s1	0.8	s2	1.0	s3	a
	}
	exec cg select {*}$dbopt -f {-*-*} tmp/temp.tsv
} {id	annot
1	a}

test select "regextract $dboptt" {
	global dbopt
	write_tab tmp/testsort.tsv {
		id	num
		gatk-rdsbwa-test1	1
		gatk-rdsbwa-test1	3
	}
#cg select -f {{id=regextract($id,{^[^-]+-[^-]+(.+)})} *} annot_seqcap_multireg-cov20.tsv.sampleinfo.pre annot_seqcap_multireg-cov20.tsv.sampleinfo
#devcg select -f {{id=regsub($id,{^[^-]+-[^-]+-},"")} *} annot_seqcap_multireg-cov20.tsv.sampleinfo.pre annot_seqcap_multireg-cov20.tsv.sampleinfo
	exec cg select {*}$dbopt -f {{id=regextract($id,{^[^-]+-[^-]+-(.+)})} *} tmp/testsort.tsv
} {id	num
test1	1
test1	3}

test select "regsub $dboptt" {
	global dbopt
	write_tab tmp/testsort.tsv {
		id	num
		gatk-rdsbwa-test1	1
		gatk-rdsbwa-test1	3
	}
	exec cg select {*}$dbopt -f {{id=regsub($id,{^[^-]+-[^-]+-},"")} *} tmp/testsort.tsv
} {id	num
test1	1
test1	3}

test select "-f compressed$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -f "num text" data/table.tsv tmp/temp.lz4
	cg lz4cat tmp/temp.lz4
} {num	text
4	a
10	b
2	c
100	d}

test select "-q compressed$dboptt" {
	global dbopt
	exec cg select {*}$dbopt -q {$num <= 4} data/table.tsv tmp/temp.lz4
	cg lz4cat tmp/temp.lz4
} {num	text	mixed	other
4	a	a4	aaaa
2	c	a2	cc}

}

foreach dbopt {{}} {
	selecttests
} ;# dbopt

test select "maximpact" {
	write_tab tmp/testsort.tsv {
		id	impact
		test1	{}
		test2	CDSsilent,RNA
	}
	exec cg select -f {id max=maximpact($impact)} tmp/testsort.tsv
} {id	max
test1	
test2	CDSsilent}

test select "maximpact" {
	write_tab tmp/testsort.tsv {
		id	impact	impact2
		test1	{}	GENEDEL
		test2	CDSsilent,RNA	CDSMIS
	}
	exec cg select -f {id max=maximpact($impact,$impact2)} tmp/testsort.tsv
} {id	max
test1	GENEDEL
test2	CDSMIS}

test select "select -n" {
	global dbopt
	write_tab tmp/test.tsv {
		chromsome begin end zyg-gatk-rdsbwa-sample1 zyg-sam-rdsbwa-sample1 zyg-gatk-rdsbwa-sample2 zyg-sam-rdsbwa-sample2
		10	1000	1001	t	m	r	t
	}
	exec cg select -n tmp/test.tsv
} {sample1
sample2}

test select "select -a" {
	global dbopt
	write_tab tmp/test.tsv {
		chromsome begin end zyg-gatk-rdsbwa-sample1 zyg-sam-rdsbwa-sample1 zyg-gatk-rdsbwa-sample2 zyg-sam-rdsbwa-sample2
		10	1000	1001	t	m	r	t
	}
	exec cg select -a tmp/test.tsv
} {gatk-rdsbwa-sample1
sam-rdsbwa-sample1
gatk-rdsbwa-sample2
sam-rdsbwa-sample2}

test select "select allways make / to float (integer division is unexpected)" {
	global dbopt
	write_tab tmp/test.tsv {
		name value1 value2
		a	2	4
		b	2.0	4
	}
	exec cg select -f {name {r=$value1 / $value2}} tmp/test.tsv
} {name	r
a	0.5
b	0.5}

test select "-compressionlevel" {
	exec cg select -compressionlevel 1 -s - data/expected-test1000glow.vcf2tsv tmp/result.tsv.lz4
} {}

test select "same() on diff analyses" {
	write_tab tmp/test.tsv {
		chromsome begin end sequenced-gatk-rdsbwa-sample2	alleleSeq1-gatk-rdsbwa-sample2	alleleSeq2-gatk-rdsbwa-sample2	sequenced-gatk-rdsbwa-sample1	alleleSeq1-gatk-rdsbwa-sample1	alleleSeq2-gatk-rdsbwa-sample1
		10	1000	10001	v	A	T	v	A	T
		10	1001	10002	r	A	A	v	A	T
		10	1001	10002	r	A	A	v	A	A
	}
	exec cg select -g all -q {same("gatk-rdsbwa-sample1","gatk-rdsbwa-sample2")} tmp/test.tsv
} {all	count
all	2}

test select "sm() on diff analyses" {
	write_tab tmp/test.tsv {
		chromsome begin end sequenced-gatk-rdsbwa-sample2	alleleSeq1-gatk-rdsbwa-sample2	alleleSeq2-gatk-rdsbwa-sample2	sequenced-gatk-rdsbwa-sample1	alleleSeq1-gatk-rdsbwa-sample1	alleleSeq2-gatk-rdsbwa-sample1
		10	1000	10001	v	A	T	v	A	T
		10	1001	10002	r	A	A	v	A	T
		10	1001	10002	r	A	A	v	A	A
	}
	exec cg select -g all -q {sm("gatk-rdsbwa-sample1","gatk-rdsbwa-sample2")} tmp/test.tsv
} {all	count
all	1}

test select "compare() on diff analyses" {
	write_tab tmp/test.tsv {
		chromsome begin end sequenced-gatk-rdsbwa-sample2	alleleSeq1-gatk-rdsbwa-sample2	alleleSeq2-gatk-rdsbwa-sample2	sequenced-gatk-rdsbwa-sample1	alleleSeq1-gatk-rdsbwa-sample1	alleleSeq2-gatk-rdsbwa-sample1
		10	1000	10001	v	A	T	v	A	T
		10	1001	10002	r	A	A	v	A	T
		10	1001	10002	r	A	A	v	A	A
		10	1001	10002	u	A	A	v	A	A
		10	1000	10001	v	A	T	v	T	T
	}
	exec cg select -g comp -f {comp=compare("gatk-rdsbwa-sample1","gatk-rdsbwa-sample2")} tmp/test.tsv
} {comp	count
df	2
mm	1
sm	1
un	1}

test select "ROW with sort" {
	file_write tmp/test.tsv [deindent {
		key	value
		b	1
		a	2
		c	3
	}]\n
	exec cg select -sh /dev/null -f ROW -s key tmp/test.tsv
} {1
0
2}

test select "ROW with sort and query" {
	file_write tmp/test.tsv [deindent {
		key	value
		b	1
		a	2
		c	3
	}]\n
	exec cg select -sh /dev/null -f ROW -q {$key ne "b"} -s key tmp/test.tsv
} {1
2}

test select "ROW in calc field with sort and query" {
	file_write tmp/test.tsv [deindent {
		key	value
		b	1
		a	2
		c	3
	}]\n
	exec cg select -sh /dev/null -f {{test=$ROW/10.0}} -q {$key ne "b"} -s key tmp/test.tsv
} {0.1
0.2}

test select "ROW in query with sort and query" {
	file_write tmp/test.tsv [deindent {
		key	value
		b	1
		a	2
		c	3
	}]\n
	exec cg select -sh /dev/null -f ROW -q {$ROW != 0} -s key tmp/test.tsv
} {1
2}

test select "cg select -s {a -b}" {
	file_write tmp/test.tsv [deindent {
		a	b
		1	1
		1	2
		2	3
		2	4
	}]\n
	exec cg select -sh /dev/null -s {a -b} tmp/test.tsv
} {1	2
1	1
2	4
2	3}

test select "cg select -s {-a b}" {
	file_write tmp/test.tsv [deindent {
		a	b
		1	1
		1	2
		2	3
		2	4
	}]\n
	exec cg select -stack 1 -sh /dev/null -s {-a b} tmp/test.tsv
} {2	3
2	4
1	1
1	2}

test select "cg select -s {-a b}" {
	file_write tmp/test.tsv [deindent {
		a	b
		1	1
		1	2
		2	3
		2	4
	}]\n
	exec cg select -stack 1 -sh /dev/null -s {-a b} tmp/test.tsv
} {2	3
2	4
1	1
1	2}

test select "cg select -s bug fix reverse sort" {
	file_write tmp/test.tsv [deindent {
		a
		10
		99856
		9939162
	}]\n
	exec cg select -stack 1 -sh /dev/null -s {-a} tmp/test.tsv
} {9939162
99856
10}

testsummarize
