#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

proc mselect_load {dbfile} {
	catch {exec cg monetdb sql test {drop table "test"}}
	exec cg tomonetdb test test $dbfile
}

test mselect {-f} {
	mselect_load data/table.tsv
	# exec cg monetdb sql test {select count(*) from "test"}
	exec cg mselect -f "num text" test test
} {num	text
4	a
10	b
2	c
100	d}

test mselect {-rf} {
	mselect_load data/table.tsv
	exec cg mselect -rf "mixed other" test test
} {num	text
4	a
10	b
2	c
100	d}

test mselect {-q} {
	mselect_load data/table.tsv
	exec cg mselect -q {$num <= 4} test test
} {num	text	mixed	other
4	a	a4	aaaa
2	c	a2	cc}

test mselect {-s} {
	mselect_load data/table.tsv
	exec cg mselect -s num test test
} {num	text	mixed	other
2	c	a2	cc
4	a	a4	aaaa
10	b	b2	bbbbbbbbbb
100	d	aa3	dd..}

test mselect {-s -f} {
	mselect_load data/table.tsv
	exec cg mselect -s num -f "num mixed" test test
} {num	mixed
2	a2
4	a4
10	b2
100	aa3}

test mselect {-s -q} {
	mselect_load data/table.tsv
	exec cg mselect -s num -q {$num <= 4} test test
} {num	text	mixed	other
2	c	a2	cc
4	a	a4	aaaa}

test mselect {-q -f} {
	mselect_load data/table.tsv
	exec cg mselect -f "num mixed" -q {$num <= 4} test test
} {num	mixed
4	a4
2	a2}

test mselect {-s -f -q} {
	mselect_load data/table.tsv
	exec cg mselect -s num  -f "num mixed" -q {$num <= 4} test test
} {num	mixed
2	a2
4	a4}

test mselect {-q multiple lines, tabs} {
mselect_load data/table.tsv
exec cg mselect -s num  -f "num mixed" -q {
	$num <= 4
	and $text == "c"
} test test
} {num	mixed
2	a2}

test mselect {-f *} {
	mselect_load data/vars1.sft
	exec cg mselect -f {chromosome begin end alleleSeq1_*} -q {$begin == 4000} test test
} {chromosome	begin	end	alleleSeq1_sample1	alleleSeq1_sample2
chr1	4000	4001	A	A
chr2	4000	4001	G	G}

test mselect {-f calculated} {
	mselect_load data/vars1.sft
	exec cg mselect -f {chromosome begin end {geno1="alleleSeq1_sample1" || '/' || "alleleSeq2_sample1"}} -q {"begin" = 4000} test test
} {chromosome	begin	end	geno1
chr1	4000	4001	A/G
chr2	4000	4001	G/G}

test mselect {-f calculated functions} {
	mselect_load data/vars1.sft
	exec cg mselect -f {chromosome begin end {countG=count($alleleSeq*, = "G")}} -q {$begin = 4000} test test
} {chromosome	begin	end	countG
chr1	4000	4001	2
chr2	4000	4001	4}

test select {-f calculated if} {
	mselect_load data/vars1.sft
	exec cg mselect -f {chromosome begin end {countG=if(count($alleleSeq*, == "G")<4,"<4",">=4")}} -q {$begin == 4000} test test
} {chromosome	begin	end	countG
chr1	4000	4001	<4
chr2	4000	4001	>=4}

#test mselect {keep header info and format rtg: -hc} {
#	exec cg mselect -hc 1 -s position data/rtgsnps.tsv tmp/temp.tsv
#	file delete temp
#	catch {exec diff tmp/temp.tsv data/rtgsnps.tsv > temp}
#	file_read temp
#} {3c3
#< name	position	type	reference	prediction	posterior	coverage	correction	support_statistics
#---
#> #name	position	type	reference	prediction	posterior	coverage	correction	support_statistics
#}

#test mselect {keep header info and format vcf} {
#	exec cg mselect -s POS data/test.vcf tmp/temp.tsv
#	file delete temp
#	catch {exec diff tmp/temp.tsv data/test.vcf > temp}
#	file_read temp
#} {18c18
#< #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003
#---
#> #CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
#}

# cannot sort on calculated fields (yet)
#test mselect {-f calculated functions + sort} {
#	exec cg mselect -f {chromosome begin end {countG=count($alleleSeq*, == "G")}} -q {$begin == 4000} -s countG [gzfile data/vars1.sft]
#} {chromosome	begin	end	countG
#chr1	4000	4001	2
#chr2	4000	4001	4}

#test groupby {groupby} {
#	exec cg groupby s1 < data/table2.tsv
#} {s1	pos	s2	s3
#u	1	u	u
#vx	2	u	u
#v	3,4,5	u,v,v	v,u,v
#u	6,7,4000000000	v,v,xx	u,v,xx}
#
#test groupby {groupby -sumfields} {
#	exec cg groupby -sumfields pos s2 < data/table2.tsv
#} {s2	pos	s1	s3
#u	6	u,vx,v	u,u,v
#v	22	v,v,u,u	u,v,u,v
#xx	4000000000	u	xx}
#
#test groupby {groupby -sumfields} {
#	exec cg groupby -sumfields pos s1 < data/table2.tsv
#} {s1	pos	s2	s3
#u	1	u	u
#vx	2	u	u
#v	12	u,v,v	v,u,v
#u	4000000013	v,v,xx	u,v,xx}
#
#test groupby {groupby -sumfields one count (non existing)} {
#	exec cg groupby -f s1 -sumfields count s1 < data/table2.tsv
#} {s1	count
#u	1
#vx	1
#v	3
#u	3}
#
#test groupby {groupby -sumfields one count (non existing) -sorted 0} {
#	exec cg groupby -sorted 0 -f s1 -sumfields count s1 < data/table2.tsv
#} {s1	count
#u	4
#v	3
#vx	1}
#
#test groupby {groupby -sorted 0} {
#	exec cg groupby -sorted 0 s1 < data/table2.tsv
#} {s1	pos	s2	s3
#u	1,6,7,4000000000	u,v,v,xx	u,u,v,xx
#v	3,4,5	u,v,v	v,u,v
#vx	2	u	u}
#
#test groupby {groupby -sorted 0 -f pos} {
#	exec cg groupby -sorted 0 -f pos s1 < data/table2.tsv
#} {s1	pos
#u	1,6,7,4000000000
#v	3,4,5
#vx	2}
#
#test groupby {groupby -sorted 0 -sumfields} {
#	exec cg groupby -sorted 0 -sumfields pos s2 < data/table2.tsv
#} {s2	pos	s1	s3
#u	6	u,vx,v	u,u,v
#v	22	v,v,u,u	u,v,u,v
#xx	4000000000	u	xx}
#
#test groupby {groupby -sorted 0 -sumfields} {
#	exec cg groupby -sorted 0 -sumfields pos s1 < data/table2.tsv
#} {s1	pos	s2	s3
#u	4000000014	u,v,v,xx	u,u,v,xx
#v	12	u,v,v	v,u,v
#vx	2	u	u}

test mselect {sm} {
	mselect_load data/vars1.sft
	split [exec cg mselect -f {chromosome begin} -q {sm("sample1","sample2")} test test] \n
} {{chromosome	begin} {chr1	4000} {chr1	4001} {chr1	4099} {chr2	4000} {chr2	4001} {chr2	4099} {chr2	5010}}

test mselect {df} {
	mselect_load data/vars1.sft
	split [exec cg mselect -f {chromosome begin} -q {df(sample1,sample2)} test test] \n
} {{chromosome	begin} {chr1	5020} {chr2	5000} {chr2	5011} {chr2	8000}}

test mselect {mm} {
	mselect_load data/vars1.sft
	split [exec cg mselect -f {chromosome begin} -q {mm(sample1,sample2)} test test] \n
} {{chromosome	begin} {chr2	5005}}

test mselect {-q count with wildcard} {
	mselect_load data/vars1.sft
	exec cg mselect -f {chromosome begin end alleleSeq*} -q {count(alleleSeq*, = 'G') > 3} test test
} {chromosome	begin	end	alleleSeq1_sample1	alleleSeq2_sample1	alleleSeq1_sample2	alleleSeq2_sample2
chr1	4001	4002	G	G	G	G
chr2	4000	4001	G	G	G	G
chr2	4001	4002	G	G	G	G}

test mselect {lmin column} {
	split [exec cg mselect -f {chromosome begin {lmin=lmin($list,10)}} < data/vars1.sft] \n
} {{chromosome	begin	lmin} {chr1	4000	4} {chr1	4001	1} {chr1	4099	1} {chr1	5000	1} {chr1	5020	3} {chr2	4000	2} {chr2	4001	2} {chr2	4099	2} {chr2	5000	3} {chr2	5000	2} {chr2	5005	4} {chr2	5010	} {chr2	5011	} {chr2	8000	}}

test mselect {lmin select} {
	split [exec cg mselect -f {chromosome begin} -q {lmin($list,10) == 2} < data/vars1.sft] \n
} {{chromosome	begin} {chr2	4000} {chr2	4001} {chr2	4099} {chr2	5000}}

test mselect {counthasone column} {
	split [exec cg mselect -f {chromosome begin {lmin=counthasone($list, ==2)}} < data/vars1.sft] \n
} {{chromosome	begin	lmin} {chr1	4000	0} {chr1	4001	1} {chr1	4099	1} {chr1	5000	1} {chr1	5020	0} {chr2	4000	1} {chr2	4001	1} {chr2	4099	1} {chr2	5000	1} {chr2	5000	0} {chr2	5005	0} {chr2	5010	0} {chr2	5011	0} {chr2	8000	0}}

testsummarize
