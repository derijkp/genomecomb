#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test select {scount basic} {
	exec cg select -f {chromosome begin 
		{test=scount($sequenced-gatk == "v")} sequenced-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	sequenced-gatk-sample1	sequenced-sam-sample1	sequenced-gatk-sample2
		1	259	2	v	v	v
		1	4001	1	v	v	r
		1	4050	1	v	u	u
		1	5000	1	v	v	r
		1	5020	1	v	r	r
		1	5020	1	r	v	v
		2	4000	2	v	v	v
		2	4001	1	v	r	r
		2	4001	1	v	r	r
		2	4010	1	u	v	v
		2	4010	1	u	v	v
		2	5010	2	v	v	v
		2	10000	1	v	r	r
		2	10000	1	r	v	v
		3	876	2	v	v	v
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {scount and} {
	exec cg select -f {chromosome begin
		{test=scount($sequenced-gatk == "v" and $freq-gatk > 0.5)} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0	0.1	0.11	0.2
		1	4001	0	0.2	0.1	0.2
		1	4050	0	0.3	?	?
		1	5000	0	0.4	0.6	0.6
		1	5020	0	0.5	?	?
		1	5020	0	?	0.4	0.5
		2	4000	2	0.6	0.6	0.6
		2	4001	1	0.8	?	?
		2	4001	1	0.7	0.01	?
		2	4010	1	?	0.8	0.8
		2	4010	1	?	0.7	0.7
		2	5010	2	0.9	0.9	0.9
		2	10000	1	0.9	?	?
		2	10000	1	?	0.9	0.9
		3	876	2	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {slist basic} {
	exec cg select -f {chromosome begin
		{test=slist($sequenced-gatk)} sequenced-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	sequenced-gatk-sample1	sequenced-sam-sample1	sequenced-gatk-sample2
		1	259	v,v	v	v	v
		1	4001	v,r	v	v	r
		1	4050	v,u	v	u	u
		1	5000	v,r	v	v	r
		1	5020	v,r	v	r	r
		1	5020	r,v	r	v	v
		2	4000	v,v	v	v	v
		2	4001	v,r	v	r	r
		2	4001	v,r	v	r	r
		2	4010	u,v	u	v	v
		2	4010	u,v	u	v	v
		2	5010	v,v	v	v	v
		2	10000	v,r	v	r	r
		2	10000	r,v	r	v	v
		3	876	v,v	v	v	v
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {slist nested} {
	exec cg select -f {chromosome begin
		{test=slist(if($sequenced-gatk != "v","u",$freq-gatk))} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0.1,0.2	0.1	0.11	0.2
		1	4001	0.2,u	0.2	0.1	0.2
		1	4050	0.3,u	0.3	?	?
		1	5000	0.4,u	0.4	0.6	0.6
		1	5020	0.5,u	0.5	?	?
		1	5020	u,0.5	?	0.4	0.5
		2	4000	0.6,0.6	0.6	0.6	0.6
		2	4001	0.8,u	0.8	?	?
		2	4001	0.7,u	0.7	0.01	?
		2	4010	u,0.8	?	0.8	0.8
		2	4010	u,0.7	?	0.7	0.7
		2	5010	0.9,0.9	0.9	0.9	0.9
		2	10000	0.9,u	0.9	?	?
		2	10000	u,0.9	?	0.9	0.9
		3	876	1,1	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {slist condition} {
	exec cg select -f {chromosome begin
		{test=slist($sequenced-gatk == "v" and $freq-gatk > 0.5,if($sequenced-gatk != "v","u",$freq-gatk))} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	{}	0.1	0.11	0.2
		1	4001	{}	0.2	0.1	0.2
		1	4050	{}	0.3	?	?
		1	5000	{}	0.4	0.6	0.6
		1	5020	{}	0.5	?	?
		1	5020	{}	?	0.4	0.5
		2	4000	0.6,0.6	0.6	0.6	0.6
		2	4001	0.8	0.8	?	?
		2	4001	0.7	0.7	0.01	?
		2	4010	0.8	?	0.8	0.8
		2	4010	0.7	?	0.7	0.7
		2	5010	0.9,0.9	0.9	0.9	0.9
		2	10000	0.9	0.9	?	?
		2	10000	0.9	?	0.9	0.9
		3	876	1,1	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {sdistinct} {
	exec cg select -f {chromosome begin 
		{test=sdistinct(if($sequenced-gatk != "v","u",$freq-gatk))} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0.1,0.2	0.1	0.11	0.2
		1	4001	0.2,u	0.2	0.1	0.2
		1	4050	0.3,u	0.3	?	?
		1	5000	0.4,u	0.4	0.6	0.6
		1	5020	0.5,u	0.5	?	?
		1	5020	u,0.5	?	0.4	0.5
		2	4000	0.6	0.6	0.6	0.6
		2	4001	0.8,u	0.8	?	?
		2	4001	0.7,u	0.7	0.01	?
		2	4010	u,0.8	?	0.8	0.8
		2	4010	u,0.7	?	0.7	0.7
		2	5010	0.9	0.9	0.9	0.9
		2	10000	0.9,u	0.9	?	?
		2	10000	u,0.9	?	0.9	0.9
		3	876	1	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {sucount} {
	exec cg select -f {chromosome begin 
		{test=sucount(if($sequenced-gatk != "v","u",$freq-gatk))} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	2	0.1	0.11	0.2
		1	4001	2	0.2	0.1	0.2
		1	4050	2	0.3	?	?
		1	5000	2	0.4	0.6	0.6
		1	5020	2	0.5	?	?
		1	5020	2	?	0.4	0.5
		2	4000	1	0.6	0.6	0.6
		2	4001	2	0.8	?	?
		2	4001	2	0.7	0.01	?
		2	4010	2	?	0.8	0.8
		2	4010	2	?	0.7	0.7
		2	5010	1	0.9	0.9	0.9
		2	10000	2	0.9	?	?
		2	10000	2	?	0.9	0.9
		3	876	1	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {sdistinct condition} {
	exec cg select -f {chromosome begin 
		{test=sdistinct($sequenced-gatk == "v" and $freq-gatk > 0.5,if($sequenced-gatk != "v","u",$freq-gatk))} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	{}	0.1	0.11	0.2
		1	4001	{}	0.2	0.1	0.2
		1	4050	{}	0.3	?	?
		1	5000	{}	0.4	0.6	0.6
		1	5020	{}	0.5	?	?
		1	5020	{}	?	0.4	0.5
		2	4000	0.6	0.6	0.6	0.6
		2	4001	0.8	0.8	?	?
		2	4001	0.7	0.7	0.01	?
		2	4010	0.8	?	0.8	0.8
		2	4010	0.7	?	0.7	0.7
		2	5010	0.9	0.9	0.9	0.9
		2	10000	0.9	0.9	?	?
		2	10000	0.9	?	0.9	0.9
		3	876	1	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {sdistinct condition} {
	exec cg select -f {chromosome begin 
		{test=sdistinct($alleleSeq2-gatk == "N",$sequenced-gatk)} alleleSeq2-* sequenced-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	alleleSeq2-gatk-sample1	alleleSeq2-sam-sample1	alleleSeq2-gatk-sample2	sequenced-gatk-sample1	sequenced-sam-sample1	sequenced-gatk-sample2
		1	259	{}	C	C	C	v	v	v
		1	4001	{}	C	C	C	v	v	r
		1	4050	{}	T	-	-	v	u	u
		1	5000	{}	{}	{}	{}	v	v	r
		1	5020	{}	A	G	G	v	r	r
		1	5020	{}	G	C	C	r	v	v
		2	4000	{}	A	A	A	v	v	v
		2	4001	{}	C	A	A	v	r	r
		2	4001	{}	C	A	A	v	r	r
		2	4010	{}	-	C	C	u	v	v
		2	4010	{}	-	C	C	u	v	v
		2	5010	v	{}	N	N	v	v	v
		2	10000	r,v	N	N	N	v	r	r
		2	10000	v,r	N	N	N	r	v	v
		3	876	{}	A	A	A	v	v	v
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {smin} {
	exec cg select -f {chromosome begin 
		{test=smin($freq-gatk)} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0.1	0.1	0.11	0.2
		1	4001	0.2	0.2	0.1	0.2
		1	4050	0.3	0.3	?	?
		1	5000	0.4	0.4	0.6	0.6
		1	5020	0.5	0.5	?	?
		1	5020	0.5	?	0.4	0.5
		2	4000	0.6	0.6	0.6	0.6
		2	4001	0.8	0.8	?	?
		2	4001	0.7	0.7	0.01	?
		2	4010	0.8	?	0.8	0.8
		2	4010	0.7	?	0.7	0.7
		2	5010	0.9	0.9	0.9	0.9
		2	10000	0.9	0.9	?	?
		2	10000	0.9	?	0.9	0.9
		3	876	1	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {smin condition} {
	exec cg select -f {chromosome begin 
		{test=smin($sequenced-gatk == "v", $freq-gatk)} sequenced-* freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	sequenced-gatk-sample1	sequenced-sam-sample1	sequenced-gatk-sample2	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0.1	v	v	v	0.1	0.11	0.2
		1	4001	0.2	v	v	r	0.2	0.1	0.2
		1	4050	0.3	v	u	u	0.3	?	?
		1	5000	0.4	v	v	r	0.4	0.6	0.6
		1	5020	0.5	v	r	r	0.5	?	?
		1	5020	0.5	r	v	v	?	0.4	0.5
		2	4000	0.6	v	v	v	0.6	0.6	0.6
		2	4001	0.8	v	r	r	0.8	?	?
		2	4001	0.7	v	r	r	0.7	0.01	?
		2	4010	0.8	u	v	v	?	0.8	0.8
		2	4010	0.7	u	v	v	?	0.7	0.7
		2	5010	0.9	v	v	v	0.9	0.9	0.9
		2	10000	0.9	v	r	r	0.9	?	?
		2	10000	0.9	r	v	v	?	0.9	0.9
		3	876	1	v	v	v	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {smax} {
	exec cg select -f {chromosome begin 
		{test=smax($freq-gatk)} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0.2	0.1	0.11	0.2
		1	4001	0.2	0.2	0.1	0.2
		1	4050	0.3	0.3	?	?
		1	5000	0.6	0.4	0.6	0.6
		1	5020	0.5	0.5	?	?
		1	5020	0.5	?	0.4	0.5
		2	4000	0.6	0.6	0.6	0.6
		2	4001	0.8	0.8	?	?
		2	4001	0.7	0.7	0.01	?
		2	4010	0.8	?	0.8	0.8
		2	4010	0.7	?	0.7	0.7
		2	5010	0.9	0.9	0.9	0.9
		2	10000	0.9	0.9	?	?
		2	10000	0.9	?	0.9	0.9
		3	876	1	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {smax condition} {
	exec cg select -f {chromosome begin 
		{test=smax(lmin($freq-gatk) < 0.2, $freq-gatk)} sequenced-* freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	sequenced-gatk-sample1	sequenced-sam-sample1	sequenced-gatk-sample2	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0.1	v	v	v	0.1	0.11	0.2
		1	4001	NaN	v	v	r	0.2	0.1	0.2
		1	4050	NaN	v	u	u	0.3	?	?
		1	5000	NaN	v	v	r	0.4	0.6	0.6
		1	5020	NaN	v	r	r	0.5	?	?
		1	5020	NaN	r	v	v	?	0.4	0.5
		2	4000	NaN	v	v	v	0.6	0.6	0.6
		2	4001	NaN	v	r	r	0.8	?	?
		2	4001	NaN	v	r	r	0.7	0.01	?
		2	4010	NaN	u	v	v	?	0.8	0.8
		2	4010	NaN	u	v	v	?	0.7	0.7
		2	5010	NaN	v	v	v	0.9	0.9	0.9
		2	10000	NaN	v	r	r	0.9	?	?
		2	10000	NaN	r	v	v	?	0.9	0.9
		3	876	NaN	v	v	v	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {ssum} {
	exec cg select -f {chromosome begin 
		{test=ssum($freq-gatk)} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0.30000000000000004	0.1	0.11	0.2
		1	4001	0.4	0.2	0.1	0.2
		1	4050	0.3	0.3	?	?
		1	5000	1.0	0.4	0.6	0.6
		1	5020	0.5	0.5	?	?
		1	5020	0.5	?	0.4	0.5
		2	4000	1.2	0.6	0.6	0.6
		2	4001	0.8	0.8	?	?
		2	4001	0.7	0.7	0.01	?
		2	4010	0.8	?	0.8	0.8
		2	4010	0.7	?	0.7	0.7
		2	5010	1.8	0.9	0.9	0.9
		2	10000	0.9	0.9	?	?
		2	10000	0.9	?	0.9	0.9
		3	876	2.0	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {ssum condition} {
	exec cg select -f {chromosome begin 
		{test=ssum($sequenced-gatk == "v", $freq-gatk)} sequenced-* freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	sequenced-gatk-sample1	sequenced-sam-sample1	sequenced-gatk-sample2	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0.30000000000000004	v	v	v	0.1	0.11	0.2
		1	4001	0.2	v	v	r	0.2	0.1	0.2
		1	4050	0.3	v	u	u	0.3	?	?
		1	5000	0.4	v	v	r	0.4	0.6	0.6
		1	5020	0.5	v	r	r	0.5	?	?
		1	5020	0.5	r	v	v	?	0.4	0.5
		2	4000	1.2	v	v	v	0.6	0.6	0.6
		2	4001	0.8	v	r	r	0.8	?	?
		2	4001	0.7	v	r	r	0.7	0.01	?
		2	4010	0.8	u	v	v	?	0.8	0.8
		2	4010	0.7	u	v	v	?	0.7	0.7
		2	5010	1.8	v	v	v	0.9	0.9	0.9
		2	10000	0.9	v	r	r	0.9	?	?
		2	10000	0.9	r	v	v	?	0.9	0.9
		3	876	2.0	v	v	v	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {savg} {
	exec cg select -f {chromosome begin 
		{test=format("%.2f",savg($freq-gatk))} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome begin test freq-gatk-sample1 freq-sam-sample1 freq-gatk-sample2
		1	259	0.15	0.1	0.11	0.2
		1	4001	0.2	0.2	0.1	0.2
		1	4050	0.3	0.3	?	?
		1	5000	0.5	0.4	0.6	0.6
		1	5020	0.5	0.5	?	?
		1	5020	0.5	?	0.4	0.5
		2	4000	0.6	0.6	0.6	0.6
		2	4001	0.8	0.8	?	?
		2	4001	0.7	0.7	0.01	?
		2	4010	0.8	?	0.8	0.8
		2	4010	0.7	?	0.7	0.7
		2	5010	0.9	0.9	0.9	0.9
		2	10000	0.9	0.9	?	?
		2	10000	0.9	?	0.9	0.9
		3	876	1.0	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {savg cond allways} {
	exec cg select -f {chromosome begin 
		{test=format("%.2f",savg(1,$freq-gatk))} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome begin test freq-gatk-sample1 freq-sam-sample1 freq-gatk-sample2
		1	259	0.15	0.1	0.11	0.2
		1	4001	0.2	0.2	0.1	0.2
		1	4050	0.3	0.3	?	?
		1	5000	0.5	0.4	0.6	0.6
		1	5020	0.5	0.5	?	?
		1	5020	0.5	?	0.4	0.5
		2	4000	0.6	0.6	0.6	0.6
		2	4001	0.8	0.8	?	?
		2	4001	0.7	0.7	0.01	?
		2	4010	0.8	?	0.8	0.8
		2	4010	0.7	?	0.7	0.7
		2	5010	0.9	0.9	0.9	0.9
		2	10000	0.9	0.9	?	?
		2	10000	0.9	?	0.9	0.9
		3	876	1.0	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {sstdev} {
	exec cg select -overwrite 1 -f {chromosome begin 
		{test=sstdev($freq-gatk)} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0.05	0.1	0.11	0.2
		1	4001	0.0	0.2	0.1	0.2
		1	4050	0.0	0.3	?	?
		1	5000	0.09999999999999998	0.4	0.6	0.6
		1	5020	0.0	0.5	?	?
		1	5020	0.0	?	0.4	0.5
		2	4000	0.0	0.6	0.6	0.6
		2	4001	0.0	0.8	?	?
		2	4001	0.0	0.7	0.01	?
		2	4010	0.0	?	0.8	0.8
		2	4010	0.0	?	0.7	0.7
		2	5010	0.0	0.9	0.9	0.9
		2	10000	0.0	0.9	?	?
		2	10000	0.0	?	0.9	0.9
		3	876	0.0	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {sstdev cond allways} {
	exec cg select -overwrite 1 -f {chromosome begin 
		{test=sstdev(1,$freq-gatk)} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0.05	0.1	0.11	0.2
		1	4001	0.0	0.2	0.1	0.2
		1	4050	0.0	0.3	?	?
		1	5000	0.09999999999999998	0.4	0.6	0.6
		1	5020	0.0	0.5	?	?
		1	5020	0.0	?	0.4	0.5
		2	4000	0.0	0.6	0.6	0.6
		2	4001	0.0	0.8	?	?
		2	4001	0.0	0.7	0.01	?
		2	4010	0.0	?	0.8	0.8
		2	4010	0.0	?	0.7	0.7
		2	5010	0.0	0.9	0.9	0.9
		2	10000	0.0	0.9	?	?
		2	10000	0.0	?	0.9	0.9
		3	876	0.0	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

# todo replace

test select {smedian} {
	exec cg select -f {chromosome begin 
		{test=smedian($freq-gatk)} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0.15000000000000002	0.1	0.11	0.2
		1	4001	0.2	0.2	0.1	0.2
		1	4050	0.3	0.3	?	?
		1	5000	0.5	0.4	0.6	0.6
		1	5020	0.5	0.5	?	?
		1	5020	0.5	?	0.4	0.5
		2	4000	0.6	0.6	0.6	0.6
		2	4001	0.8	0.8	?	?
		2	4001	0.7	0.7	0.01	?
		2	4010	0.8	?	0.8	0.8
		2	4010	0.7	?	0.7	0.7
		2	5010	0.9	0.9	0.9	0.9
		2	10000	0.9	0.9	?	?
		2	10000	0.9	?	0.9	0.9
		3	876	1.0	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}


test select {smedian cond allways} {
	exec cg select -f {chromosome begin 
		{test=smedian(1,$freq-gatk)} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0.15000000000000002	0.1	0.11	0.2
		1	4001	0.2	0.2	0.1	0.2
		1	4050	0.3	0.3	?	?
		1	5000	0.5	0.4	0.6	0.6
		1	5020	0.5	0.5	?	?
		1	5020	0.5	?	0.4	0.5
		2	4000	0.6	0.6	0.6	0.6
		2	4001	0.8	0.8	?	?
		2	4001	0.7	0.7	0.01	?
		2	4010	0.8	?	0.8	0.8
		2	4010	0.7	?	0.7	0.7
		2	5010	0.9	0.9	0.9	0.9
		2	10000	0.9	0.9	?	?
		2	10000	0.9	?	0.9	0.9
		3	876	1.0	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {smode} {
	exec cg select -f {chromosome begin 
		{test=smode($freq-gatk)} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0.2,0.1	0.1	0.11	0.2
		1	4001	0.2	0.2	0.1	0.2
		1	4050	0.3,?	0.3	?	?
		1	5000	0.6,0.4	0.4	0.6	0.6
		1	5020	?,0.5	0.5	?	?
		1	5020	0.5,?	?	0.4	0.5
		2	4000	0.6	0.6	0.6	0.6
		2	4001	0.8,?	0.8	?	?
		2	4001	0.7,?	0.7	0.01	?
		2	4010	0.8,?	?	0.8	0.8
		2	4010	0.7,?	?	0.7	0.7
		2	5010	0.9	0.9	0.9	0.9
		2	10000	?,0.9	0.9	?	?
		2	10000	0.9,?	?	0.9	0.9
		3	876	1	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {smode cond allways} {
	exec cg select -f {chromosome begin 
		{test=smode(1,$freq-gatk)} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0.2,0.1	0.1	0.11	0.2
		1	4001	0.2	0.2	0.1	0.2
		1	4050	0.3,?	0.3	?	?
		1	5000	0.6,0.4	0.4	0.6	0.6
		1	5020	?,0.5	0.5	?	?
		1	5020	0.5,?	?	0.4	0.5
		2	4000	0.6	0.6	0.6	0.6
		2	4001	0.8,?	0.8	?	?
		2	4001	0.7,?	0.7	0.01	?
		2	4010	0.8,?	?	0.8	0.8
		2	4010	0.7,?	?	0.7	0.7
		2	5010	0.9	0.9	0.9	0.9
		2	10000	?,0.9	0.9	?	?
		2	10000	0.9,?	?	0.9	0.9
		3	876	1	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {spercent} {
	exec cg select -f {chromosome begin 
		{test=spercent($sequenced-gatk == "v", lmax($freq-gatk) > 0.1)} sequenced-* freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	sequenced-gatk-sample1	sequenced-sam-sample1	sequenced-gatk-sample2	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	50.0	v	v	v	0.1	0.11	0.2
		1	4001	100.0	v	v	r	0.2	0.1	0.2
		1	4050	100.0	v	u	u	0.3	?	?
		1	5000	100.0	v	v	r	0.4	0.6	0.6
		1	5020	100.0	v	r	r	0.5	?	?
		1	5020	100.0	r	v	v	?	0.4	0.5
		2	4000	100.0	v	v	v	0.6	0.6	0.6
		2	4001	100.0	v	r	r	0.8	?	?
		2	4001	100.0	v	r	r	0.7	0.01	?
		2	4010	100.0	u	v	v	?	0.8	0.8
		2	4010	100.0	u	v	v	?	0.7	0.7
		2	5010	100.0	v	v	v	0.9	0.9	0.9
		2	10000	100.0	v	r	r	0.9	?	?
		2	10000	100.0	r	v	v	?	0.9	0.9
		3	876	100.0	v	v	v	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {scount sample} {
	exec cg select -f {chromosome begin 
		{test=scount(($sample eq "sample1") and $sequenced-gatk == "v")} sequenced-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	sequenced-gatk-sample1	sequenced-sam-sample1	sequenced-gatk-sample2
		1	259	1	v	v	v
		1	4001	1	v	v	r
		1	4050	1	v	u	u
		1	5000	1	v	v	r
		1	5020	1	v	r	r
		1	5020	0	r	v	v
		2	4000	1	v	v	v
		2	4001	1	v	r	r
		2	4001	1	v	r	r
		2	4010	0	u	v	v
		2	4010	0	u	v	v
		2	5010	1	v	v	v
		2	10000	1	v	r	r
		2	10000	0	r	v	v
		3	876	1	v	v	v
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {scount sample in query} {
	exec cg select \
		-q {scount(($sample eq "sample1") and $sequenced-gatk == "v") > 0} \
		-f {chromosome begin sequenced-*} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	sequenced-gatk-sample1	sequenced-sam-sample1	sequenced-gatk-sample2
		1	259	v	v	v
		1	4001	v	v	r
		1	4050	v	u	u
		1	5000	v	v	r
		1	5020	v	r	r
		2	4000	v	v	v
		2	4001	v	r	r
		2	4001	v	r	r
		2	5010	v	v	v
		2	10000	v	r	r
		3	876	v	v	v
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {scount sample in query and fields} {
	exec cg select \
		-q {scount(($sample eq "sample1") and $sequenced-gatk == "v") > 0} \
		-f {chromosome begin {test=scount(($sample eq "sample1") and $sequenced-gatk == "v")} sequenced-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	sequenced-gatk-sample1	sequenced-sam-sample1	sequenced-gatk-sample2
		1	259	1	v	v	v
		1	4001	1	v	v	r
		1	4050	1	v	u	u
		1	5000	1	v	v	r
		1	5020	1	v	r	r
		2	4000	1	v	v	v
		2	4001	1	v	r	r
		2	4001	1	v	r	r
		2	5010	1	v	v	v
		2	10000	1	v	r	r
		3	876	1	v	v	v
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {incorrect use of sample outside sample aggregates} {
	exec cg select \
		-q {$sample eq "gatk-sample1"} \
		-f {chromosome begin sequenced-*
	} data/vars-saggr.tsv
} {field sample not present in file (or sampleinfo)} error

test select {scount field not present in all samples} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1	sequenced-sam-sample1 zyg-gatk-sample2
		chr1 4001 4002 snp A G,C	v	v	t
	}
	exec cg select -f {chromosome begin {test=scount($sequenced-gatk eq "v")}} tmp/vars.tsv
} {chromosome	begin	test
chr1	4001	1}

test select {scount field not present at all} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1	sequenced-sam-sample1 zyg-gatk-sample2
		chr1 4001 4002 snp A G,C	v	v	t
	}
	exec cg select -f {chromosome begin {test=scount($test-gatk eq "v")}} tmp/vars.tsv
} {error in scount: all samples are missing one or more needed fields (test-gatk-sample1,test-gatk-sample2)} error

test select {scount used in other function field not present in all samples} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1	sequenced-sam-sample1 zyg-gatk-sample2
		chr1 4001 4002 snp A G,C	v	v	t
	}
	exec cg select -f {chromosome begin {test=if(scount($sequenced-gatk eq "v") == 1, scount($sequenced-gatk eq "v")+20, scount($sequenced-gatk eq "v")+10)}} tmp/vars.tsv
} {chromosome	begin	test
chr1	4001	21}

test select {scount used in other function field not present in all samples, but in calccols} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1	sequenced-sam-sample1 zyg-gatk-sample2
		chr1 4001 4002 snp A G,C	v	v	t
	}
	exec cg select -f {chromosome begin -sequenced-gatk-sample2="v" {test=if(scount($sequenced-gatk eq "v") == 2, scount($sequenced-gatk eq "v")+20, scount($sequenced-gatk eq "v")+10)}} tmp/vars.tsv
} {chromosome	begin	test
chr1	4001	22}

test select {scount field not present in all samples, but in calccol} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1	sequenced-sam-sample1 zyg-gatk-sample2
		chr1 4001 4002 snp A G,C	v	v	t
	}
	exec cg select -stack 1 -f {chromosome begin -sequenced-gatk-sample2="v" {test=scount($sequenced-gatk eq "v")}} tmp/vars.tsv
} {chromosome	begin	test
chr1	4001	2}

test select {scount field not present in all samples, but in sampleinfo} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1	sequenced-sam-sample1 zyg-gatk-sample2
		chr1 4001 4002 snp A G,C	v	v	t
	}
	write_tab tmp/vars.tsv.sampleinfo.tsv {
		sample	sequenced
		gatk-sample2	v
	}
	exec cg select -stack 1 -f {chromosome begin {test=scount($sequenced-gatk eq "v")}} tmp/vars.tsv
} {chromosome	begin	test
chr1	4001	2}

test select {spercent field not present in all samples} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1	sequenced-sam-sample1 zyg-gatk-sample2 sequenced-gatk-sample3
		chr1 4001 4002 snp A G,C	v	r	t r
	}
	exec cg select -f {chromosome begin {test=spercent($sequenced-gatk eq "v")}} tmp/vars.tsv
} {chromosome	begin	test
chr1	4001	50.0}

test select {spercent field not present in all samples, 2 cond} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1 num-gatk-sample1	sequenced-sam-sample1 num-sam-sample1	zyg-gatk-sample2
		chr1 4001 4002 snp A G,C	v 1 t 0 t
	}
	exec cg select -f {chromosome begin {test=spercent($sequenced-gatk eq "v",$num-gatk == 1)}} tmp/vars.tsv
} {chromosome	begin	test
chr1	4001	100.0}

test select {spercent field not present in all samples, 2 cond different missing} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1 num-gatk-sample1	sequenced-sam-sample1	zyg-gatk-sample2 num-gatk-sample2 sequenced-gatk-sample3
		chr1 4001 4002 snp A G,C	v 1 t t 0 v
	}
	exec cg select -f {chromosome begin {test=spercent($sequenced-gatk eq "v",$num-gatk == 1)}} tmp/vars.tsv
} {chromosome	begin	test
chr1	4001	100.0}

test select {spercent field not present at all} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1	sequenced-sam-sample1 zyg-gatk-sample2
		chr1 4001 4002 snp A G,C	v	t	t
	}
	exec cg select -f {chromosome begin {test=spercent($test-gatk eq "v")}} tmp/vars.tsv
} {error in spercent: all samples are missing one or more needed fields (test-gatk-sample1,test-gatk-sample2)} error

test select {ssum field not present in all samples} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1	sequenced-sam-sample1 zyg-gatk-sample2
		chr1 4001 4002 snp A G,C	v	v	t
	}
	exec cg select -f {chromosome begin {test=ssum($sequenced-gatk eq "v")}} tmp/vars.tsv
} {chromosome	begin	test
chr1	4001	1.0}

test select {ssum field not present in all samples, 2 cond} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1 num-gatk-sample1	sequenced-sam-sample1 num-sam-sample1	zyg-gatk-sample2 sequenced-gatk-sample3 num-gatk-sample3
		chr1 4001 4002 snp A G,C	v 1 v 2 t v 3
	}
	exec cg select -f {chromosome begin {test=ssum($sequenced-gatk eq "v",$num-gatk)}} tmp/vars.tsv
} {chromosome	begin	test
chr1	4001	4.0}

test select {ssum field not present in all samples, 2 cond different missing} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1 num-gatk-sample1	sequenced-sam-sample1	zyg-gatk-sample2 num-gatk-sample2
		chr1 4001 4002 snp A G,C	v 1 v t 0
	}
	exec cg select -f {chromosome begin {test=ssum($sequenced-gatk eq "v",$num-gatk)}} tmp/vars.tsv
} {chromosome	begin	test
chr1	4001	1.0}

test select {ssum field not present in all samples, 2 cond different missing but in calc and sampleinfo} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1 num-gatk-sample1	sequenced-sam-sample1	zyg-gatk-sample2 num-gatk-sample2 sequenced-gatk-sample3
		chr1 4001 4002 snp A G,C	v 1 v t 3 v
	}
	write_tab tmp/vars.tsv.sampleinfo.tsv {
		sample	num
		gatk-sample3	2
	}
	exec cg select -f {chromosome begin -sequenced-gatk-sample2="v" {test=ssum($sequenced-gatk eq "v",$num-gatk)}} tmp/vars.tsv
} {chromosome	begin	test
chr1	4001	6.0}

test select {ssum field not present at all} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1 num-gatk-sample1	sequenced-sam-sample1 num-sam-sample1	zyg-gatk-sample2
		chr1 4001 4002 snp A G,C	v 1 t 2 t
	}
	exec cg select -f {chromosome begin {test=ssum($test-gatk eq "v",$num)}} tmp/vars.tsv
} {error in ssum: all samples are missing one or more needed fields (test-gatk-sample1,test-gatk-sample2)} error

test select {ssum field not present at all} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1 num-gatk-sample1	sequenced-sam-sample1 num-sam-sample1	zyg-gatk-sample2
		chr1 4001 4002 snp A G,C	v 1 t 2 t
	}
	exec cg select -f {chromosome begin {test=ssum($sequenced-gatk eq "v",$test-gatk)}} tmp/vars.tsv
} {error in ssum: all samples are missing one or more needed fields (test-gatk-sample1,sequenced-gatk-sample2)} error

test select {calc field in -f used in query scount} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	zyg-gatk-sample1 zyg-sam-sample1	zyg-gatk-sample2 zyg-sam-sample2
		chr1 4001 4002 snp A G,C	m m t m
		chr1 4001 4002 snp A G,C	m m t t
	}
	exec cg select -f {chromosome begin {samzyg-gatk-*=$zyg-sam-*}} -q {scount($samzyg-gatk eq "m") == 2} tmp/vars.tsv
} {chromosome	begin	samzyg-gatk-sample1	samzyg-gatk-sample2
chr1	4001	m	m}

test select {query between analyses} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	zyg-gatk-sample1 zyg-sam-sample1	zyg-gatk-sample2 zyg-sam-sample2
		chr1 4001 4002 snp A G,C	m m t m
		chr1 4002 4003 snp A G,C	m m m m
	}
	exec cg select -f {chromosome begin {scount=scount($zyg-gatk eq $zyg-sam)}} tmp/vars.tsv
} {chromosome	begin	scount
chr1	4001	1
chr1	4002	2}

test select {use sample in scount condition} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	zyg-gatk-samplea1	zyg-gatk-samplea2	zyg-sam-samplea2	zyg-gatk-sampleb3
		chr1	4001	4002	snp	A	G	t	r	t	t
	}
	exec cg select -f {chromosome begin {count=scount($zyg-gatk == "t" and $sample matches "samplea*")}} tmp/vars.tsv
} {chromosome	begin	count
chr1	4001	1}

test select {trying to use analysis in scount condition error} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	zyg-gatk-samplea1	zyg-gatk-samplea2	zyg-sam-samplea2	zyg-gatk-sampleb3
		chr1	4001	4002	snp	A	G	t	r	t	t
	}
	exec cg select -f {chromosome begin {count=scount($zyg-gatk == "t" and $analysis matches "gatk-samplea*")}} tmp/vars.tsv
} {error in scount: all samples are missing one or more needed fields (analysis-samplea1,analysis-samplea2,analysis-sampleb3)} error

test select {use sample in slist condition} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	zyg-gatk-samplea1	zyg-gatk-samplea2	zyg-sam-samplea2	zyg-gatk-sampleb3
		chr1	4001	4002	snp	A	G	t	r	t	t
	}
	exec cg select -f {chromosome begin {list=slist($sample matches "samplea*",$zyg-gatk)}} tmp/vars.tsv
} {chromosome	begin	list
chr1	4001	t,r}

testsummarize
