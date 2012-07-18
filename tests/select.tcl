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
} {{chromosome	begin	lmin} {chr1	4000	4} {chr1	4001	1} {chr1	4099	1} {chr1	5000	1} {chr1	5020	3} {chr2	4000	2} {chr2	4001	2} {chr2	4099	2} {chr2	5000	2} {chr2	5000	3} {chr2	5005	4} {chr2	5010	20} {chr2	5011	Inf} {chr2	8000	Inf}}

test select {lmind column} {
	split [exec cg select -f {chromosome begin {lmin=lmind($list,10)}} < data/vars1.sft] \n
} {{chromosome	begin	lmin} {chr1	4000	4} {chr1	4001	1} {chr1	4099	1} {chr1	5000	1} {chr1	5020	3} {chr2	4000	2} {chr2	4001	2} {chr2	4099	2} {chr2	5000	2} {chr2	5000	3} {chr2	5005	4} {chr2	5010	10} {chr2	5011	10} {chr2	8000	Inf}}

test select {lmind select} {
	split [exec cg select -f {chromosome begin} -q {lmind($list,10) == 2} < data/vars1.sft] \n
} {{chromosome	begin} {chr2	4000} {chr2	4001} {chr2	4099} {chr2	5000}}

test select {lmax column} {
	split [exec cg select -f {chromosome begin {lmax=lmax($list)}} < data/vars1.sft] \n
} {{chromosome	begin	lmax} {chr1	4000	4} {chr1	4001	4} {chr1	4099	2} {chr1	5000	2} {chr1	5020	3} {chr2	4000	2} {chr2	4001	4} {chr2	4099	4} {chr2	5000	2} {chr2	5000	3} {chr2	5005	4} {chr2	5010	20} {chr2	5011	-Inf} {chr2	8000	-Inf}}

test select {lmaxd column} {
	split [exec cg select -f {chromosome begin {lmax=lmaxd($list,0)}} < data/vars1.sft] \n
} {{chromosome	begin	lmax} {chr1	4000	4} {chr1	4001	4} {chr1	4099	2} {chr1	5000	2} {chr1	5020	3} {chr2	4000	2} {chr2	4001	4} {chr2	4099	4} {chr2	5000	2} {chr2	5000	3} {chr2	5005	4} {chr2	5010	20} {chr2	5011	0} {chr2	8000	-Inf}}

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
	split [exec cg select -f {chromosome begin end row=$ROW} -q {$ROW >= 2 && $ROW <= 3} < data/vars1.sft] \n
} {{chromosome	begin	end	row} {chr1	4099	5000	2} {chr1	5000	5010	3}}

test select {shared objects bugcheck} {
	split [exec cg select -f {chromosome begin end {test=[set ::keep $begin]}} -q {$ROW between {2 3}} < data/vars1.sft] \n
} {{chromosome	begin	end	test} {chr1	4099	5000	4099} {chr1	5000	5010	5000}}

test select {shared objects bugcheck} {
	split [exec cg select -f {chromosome begin end {test="[get ::keep 1]-[set ::keep $begin]"}} -q {$ROW between {2 3}} < data/vars1.sft] \n
} {{chromosome	begin	end	test} {chr1	4099	5000	1-4099} {chr1	5000	5010	4099-5000}}

test select {start brace bugcheck} {
	split [exec cg select -f {chromosome begin end {type=($type == "snp")? "Snp" : (($type == "del")? "Deletion" : $type)}} -q {$ROW in {2 3 9 10}} < data/vars1.sft] \n
} {{chromosome	begin	end	type} {chr1	4099	5000	Snp} {chr1	5000	5010	Deletion} {chr2	5000	5010	ins} {chr2	5005	5006	Snp}}

testsummarize
