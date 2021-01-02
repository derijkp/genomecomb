#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test select {acount basic} {
	exec cg select -f {chromosome begin 
		{test=acount($sequenced == "v")} sequenced-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	sequenced-gatk-sample1	sequenced-sam-sample1	sequenced-gatk-sample2
		1	259	3	v	v	v
		1	4001	2	v	v	r
		1	4050	1	v	u	u
		1	5000	2	v	v	r
		1	5020	1	v	r	r
		1	5020	2	r	v	v
		2	4000	3	v	v	v
		2	4001	1	v	r	r
		2	4001	1	v	r	r
		2	4010	2	u	v	v
		2	4010	2	u	v	v
		2	5010	3	v	v	v
		2	10000	1	v	r	r
		2	10000	2	r	v	v
		3	876	3	v	v	v
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {acount and} {
	exec cg select -f {chromosome begin
		{test=acount($sequenced == "v" and $freq > 0.5)} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0	0.1	0.11	0.2
		1	4001	0	0.2	0.1	0.2
		1	4050	0	0.3	?	?
		1	5000	1	0.4	0.6	0.6
		1	5020	0	0.5	?	?
		1	5020	0	?	0.4	0.5
		2	4000	3	0.6	0.6	0.6
		2	4001	1	0.8	?	?
		2	4001	1	0.7	0.01	?
		2	4010	2	?	0.8	0.8
		2	4010	2	?	0.7	0.7
		2	5010	3	0.9	0.9	0.9
		2	10000	1	0.9	?	?
		2	10000	2	?	0.9	0.9
		3	876	3	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {alist basic} {
	exec cg select -f {chromosome begin
		{test=alist($sequenced)} sequenced-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	sequenced-gatk-sample1	sequenced-sam-sample1	sequenced-gatk-sample2
		1	259	v,v,v	v	v	v
		1	4001	v,v,r	v	v	r
		1	4050	v,u,u	v	u	u
		1	5000	v,v,r	v	v	r
		1	5020	v,r,r	v	r	r
		1	5020	r,v,v	r	v	v
		2	4000	v,v,v	v	v	v
		2	4001	v,r,r	v	r	r
		2	4001	v,r,r	v	r	r
		2	4010	u,v,v	u	v	v
		2	4010	u,v,v	u	v	v
		2	5010	v,v,v	v	v	v
		2	10000	v,r,r	v	r	r
		2	10000	r,v,v	r	v	v
		3	876	v,v,v	v	v	v
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {alist nested} {
	exec cg select -f {chromosome begin
		{test=alist(if($sequenced != "v","u",$freq))} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0.1,0.11,0.2	0.1	0.11	0.2
		1	4001	0.2,0.1,u	0.2	0.1	0.2
		1	4050	0.3,u,u	0.3	?	?
		1	5000	0.4,0.6,u	0.4	0.6	0.6
		1	5020	0.5,u,u	0.5	?	?
		1	5020	u,0.4,0.5	?	0.4	0.5
		2	4000	0.6,0.6,0.6	0.6	0.6	0.6
		2	4001	0.8,u,u	0.8	?	?
		2	4001	0.7,u,u	0.7	0.01	?
		2	4010	u,0.8,0.8	?	0.8	0.8
		2	4010	u,0.7,0.7	?	0.7	0.7
		2	5010	0.9,0.9,0.9	0.9	0.9	0.9
		2	10000	0.9,u,u	0.9	?	?
		2	10000	u,0.9,0.9	?	0.9	0.9
		3	876	1,1,1	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {alist condition} {
	exec cg select -f {chromosome begin
		{test=alist($sequenced == "v" and $freq > 0.5,if($sequenced != "v","u",$freq))} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	{}	0.1	0.11	0.2
		1	4001	{}	0.2	0.1	0.2
		1	4050	{}	0.3	?	?
		1	5000	0.6	0.4	0.6	0.6
		1	5020	{}	0.5	?	?
		1	5020	{}	?	0.4	0.5
		2	4000	0.6,0.6,0.6	0.6	0.6	0.6
		2	4001	0.8	0.8	?	?
		2	4001	0.7	0.7	0.01	?
		2	4010	0.8,0.8	?	0.8	0.8
		2	4010	0.7,0.7	?	0.7	0.7
		2	5010	0.9,0.9,0.9	0.9	0.9	0.9
		2	10000	0.9	0.9	?	?
		2	10000	0.9,0.9	?	0.9	0.9
		3	876	1,1,1	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {adistinct} {
	exec cg select -f {chromosome begin 
		{test=adistinct(if($sequenced != "v","u",$freq))} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0.1,0.11,0.2	0.1	0.11	0.2
		1	4001	0.2,0.1,u	0.2	0.1	0.2
		1	4050	0.3,u	0.3	?	?
		1	5000	0.4,0.6,u	0.4	0.6	0.6
		1	5020	0.5,u	0.5	?	?
		1	5020	u,0.4,0.5	?	0.4	0.5
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

test select {aucount} {
	exec cg select -f {chromosome begin 
		{test=aucount(if($sequenced != "v","u",$freq))} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	3	0.1	0.11	0.2
		1	4001	3	0.2	0.1	0.2
		1	4050	2	0.3	?	?
		1	5000	3	0.4	0.6	0.6
		1	5020	2	0.5	?	?
		1	5020	3	?	0.4	0.5
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

test select {adistinct condition} {
	exec cg select -f {chromosome begin 
		{test=adistinct($sequenced == "v" and $freq > 0.5,if($sequenced != "v","u",$freq))} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	{}	0.1	0.11	0.2
		1	4001	{}	0.2	0.1	0.2
		1	4050	{}	0.3	?	?
		1	5000	0.6	0.4	0.6	0.6
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

test select {adistinct condition} {
	exec cg select -f {chromosome begin 
		{test=adistinct($alleleSeq2 == "N",$sequenced)} alleleSeq2-* sequenced-*
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

test select {amin} {
	exec cg select -f {chromosome begin 
		{test=amin($freq)} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0.1	0.1	0.11	0.2
		1	4001	0.1	0.2	0.1	0.2
		1	4050	0.3	0.3	?	?
		1	5000	0.4	0.4	0.6	0.6
		1	5020	0.5	0.5	?	?
		1	5020	0.4	?	0.4	0.5
		2	4000	0.6	0.6	0.6	0.6
		2	4001	0.8	0.8	?	?
		2	4001	0.01	0.7	0.01	?
		2	4010	0.8	?	0.8	0.8
		2	4010	0.7	?	0.7	0.7
		2	5010	0.9	0.9	0.9	0.9
		2	10000	0.9	0.9	?	?
		2	10000	0.9	?	0.9	0.9
		3	876	1	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {amin condition} {
	exec cg select -f {chromosome begin 
		{test=amin($sequenced == "v", $freq)} sequenced-* freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	sequenced-gatk-sample1	sequenced-sam-sample1	sequenced-gatk-sample2	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0.1	v	v	v	0.1	0.11	0.2
		1	4001	0.1	v	v	r	0.2	0.1	0.2
		1	4050	0.3	v	u	u	0.3	?	?
		1	5000	0.4	v	v	r	0.4	0.6	0.6
		1	5020	0.5	v	r	r	0.5	?	?
		1	5020	0.4	r	v	v	?	0.4	0.5
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

test select {amax} {
	exec cg select -f {chromosome begin 
		{test=amax($freq)} freq-*
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

test select {amax condition} {
	exec cg select -f {chromosome begin 
		{test=amax(lmin($freq) < 0.2, $freq)} sequenced-* freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	sequenced-gatk-sample1	sequenced-sam-sample1	sequenced-gatk-sample2	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0.11	v	v	v	0.1	0.11	0.2
		1	4001	0.1	v	v	r	0.2	0.1	0.2
		1	4050	NaN	v	u	u	0.3	?	?
		1	5000	NaN	v	v	r	0.4	0.6	0.6
		1	5020	NaN	v	r	r	0.5	?	?
		1	5020	NaN	r	v	v	?	0.4	0.5
		2	4000	NaN	v	v	v	0.6	0.6	0.6
		2	4001	NaN	v	r	r	0.8	?	?
		2	4001	0.01	v	r	r	0.7	0.01	?
		2	4010	NaN	u	v	v	?	0.8	0.8
		2	4010	NaN	u	v	v	?	0.7	0.7
		2	5010	NaN	v	v	v	0.9	0.9	0.9
		2	10000	NaN	v	r	r	0.9	?	?
		2	10000	NaN	r	v	v	?	0.9	0.9
		3	876	NaN	v	v	v	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {asum} {
	exec cg select -f {chromosome begin 
		{test=asum($freq)} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0.41000000000000003	0.1	0.11	0.2
		1	4001	0.5	0.2	0.1	0.2
		1	4050	0.3	0.3	?	?
		1	5000	1.6	0.4	0.6	0.6
		1	5020	0.5	0.5	?	?
		1	5020	0.9	?	0.4	0.5
		2	4000	1.7999999999999998	0.6	0.6	0.6
		2	4001	0.8	0.8	?	?
		2	4001	0.71	0.7	0.01	?
		2	4010	1.6	?	0.8	0.8
		2	4010	1.4	?	0.7	0.7
		2	5010	2.7	0.9	0.9	0.9
		2	10000	0.9	0.9	?	?
		2	10000	1.8	?	0.9	0.9
		3	876	3.0	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {asum condition} {
	exec cg select -f {chromosome begin 
		{test=asum($sequenced == "v", $freq)} sequenced-* freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	sequenced-gatk-sample1	sequenced-sam-sample1	sequenced-gatk-sample2	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0.41000000000000003	v	v	v	0.1	0.11	0.2
		1	4001	0.30000000000000004	v	v	r	0.2	0.1	0.2
		1	4050	0.3	v	u	u	0.3	?	?
		1	5000	1.0	v	v	r	0.4	0.6	0.6
		1	5020	0.5	v	r	r	0.5	?	?
		1	5020	0.9	r	v	v	?	0.4	0.5
		2	4000	1.7999999999999998	v	v	v	0.6	0.6	0.6
		2	4001	0.8	v	r	r	0.8	?	?
		2	4001	0.7	v	r	r	0.7	0.01	?
		2	4010	1.6	u	v	v	?	0.8	0.8
		2	4010	1.4	u	v	v	?	0.7	0.7
		2	5010	2.7	v	v	v	0.9	0.9	0.9
		2	10000	0.9	v	r	r	0.9	?	?
		2	10000	1.8	r	v	v	?	0.9	0.9
		3	876	3.0	v	v	v	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {aavg} {
	exec cg select -f {chromosome begin 
		{test=format("%.2f",aavg($freq))} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0.14	0.1	0.11	0.2
		1	4001	0.17	0.2	0.1	0.2
		1	4050	0.3	0.3	?	?
		1	5000	0.53	0.4	0.6	0.6
		1	5020	0.5	0.5	?	?
		1	5020	0.45	?	0.4	0.5
		2	4000	0.6	0.6	0.6	0.6
		2	4001	0.8	0.8	?	?
		2	4001	0.35	0.7	0.01	?
		2	4010	0.8	?	0.8	0.8
		2	4010	0.7	?	0.7	0.7
		2	5010	0.9	0.9	0.9	0.9
		2	10000	0.9	0.9	?	?
		2	10000	0.9	?	0.9	0.9
		3	876	1.0	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {aavg cond allways} {
	exec cg select -f {chromosome begin 
		{test=format("%.2f",aavg(1,$freq))} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0.14	0.1	0.11	0.2
		1	4001	0.17	0.2	0.1	0.2
		1	4050	0.3	0.3	?	?
		1	5000	0.53	0.4	0.6	0.6
		1	5020	0.5	0.5	?	?
		1	5020	0.45	?	0.4	0.5
		2	4000	0.6	0.6	0.6	0.6
		2	4001	0.8	0.8	?	?
		2	4001	0.35	0.7	0.01	?
		2	4010	0.8	?	0.8	0.8
		2	4010	0.7	?	0.7	0.7
		2	5010	0.9	0.9	0.9	0.9
		2	10000	0.9	0.9	?	?
		2	10000	0.9	?	0.9	0.9
		3	876	1.0	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {astdev} {
	exec cg select -overwrite 1 -f {chromosome begin 
		{test=astdev($freq)} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0.044969125210773474	0.1	0.11	0.2
		1	4001	0.04714045207910317	0.2	0.1	0.2
		1	4050	0.0	0.3	?	?
		1	5000	0.09428090415820632	0.4	0.6	0.6
		1	5020	0.0	0.5	?	?
		1	5020	0.04999999999999999	?	0.4	0.5
		2	4000	0.0	0.6	0.6	0.6
		2	4001	0.0	0.8	?	?
		2	4001	0.345	0.7	0.01	?
		2	4010	0.0	?	0.8	0.8
		2	4010	0.0	?	0.7	0.7
		2	5010	0.0	0.9	0.9	0.9
		2	10000	0.0	0.9	?	?
		2	10000	0.0	?	0.9	0.9
		3	876	0.0	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {astdev cond allways} {
	exec cg select -overwrite 1 -f {chromosome begin 
		{test=astdev(1,$freq)} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome begin test freq-gatk-sample1 freq-sam-sample1 freq-gatk-sample2
		1	259	0.044969125210773474	0.1	0.11	0.2
		1	4001	0.04714045207910317	0.2	0.1	0.2
		1	4050	0.0	0.3	?	?
		1	5000	0.09428090415820632	0.4	0.6	0.6
		1	5020	0.0	0.5	?	?
		1	5020	0.04999999999999999	?	0.4	0.5
		2	4000	0.0	0.6	0.6	0.6
		2	4001	0.0	0.8	?	?
		2	4001	0.345	0.7	0.01	?
		2	4010	0.0	?	0.8	0.8
		2	4010	0.0	?	0.7	0.7
		2	5010	0.0	0.9	0.9	0.9
		2	10000	0.0	0.9	?	?
		2	10000	0.0	?	0.9	0.9
		3	876	0.0	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {amedian} {
	exec cg select -f {chromosome begin 
		{test=amedian($freq)} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0.11	0.1	0.11	0.2
		1	4001	0.2	0.2	0.1	0.2
		1	4050	0.3	0.3	?	?
		1	5000	0.6	0.4	0.6	0.6
		1	5020	0.5	0.5	?	?
		1	5020	0.45	?	0.4	0.5
		2	4000	0.6	0.6	0.6	0.6
		2	4001	0.8	0.8	?	?
		2	4001	0.355	0.7	0.01	?
		2	4010	0.8	?	0.8	0.8
		2	4010	0.7	?	0.7	0.7
		2	5010	0.9	0.9	0.9	0.9
		2	10000	0.9	0.9	?	?
		2	10000	0.9	?	0.9	0.9
		3	876	1	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}


test select {amedian cond allways} {
	exec cg select -f {chromosome begin 
		{test=amedian(1,$freq)} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0.11	0.1	0.11	0.2
		1	4001	0.2	0.2	0.1	0.2
		1	4050	0.3	0.3	?	?
		1	5000	0.6	0.4	0.6	0.6
		1	5020	0.5	0.5	?	?
		1	5020	0.45	?	0.4	0.5
		2	4000	0.6	0.6	0.6	0.6
		2	4001	0.8	0.8	?	?
		2	4001	0.355	0.7	0.01	?
		2	4010	0.8	?	0.8	0.8
		2	4010	0.7	?	0.7	0.7
		2	5010	0.9	0.9	0.9	0.9
		2	10000	0.9	0.9	?	?
		2	10000	0.9	?	0.9	0.9
		3	876	1	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {amode} {
	exec cg select -f {chromosome begin 
		{test=amode($freq)} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0.2,0.11,0.1	0.1	0.11	0.2
		1	4001	0.2	0.2	0.1	0.2
		1	4050	?	0.3	?	?
		1	5000	0.6	0.4	0.6	0.6
		1	5020	?	0.5	?	?
		1	5020	0.4,0.5,?	?	0.4	0.5
		2	4000	0.6	0.6	0.6	0.6
		2	4001	?	0.8	?	?
		2	4001	0.7,?,0.01	0.7	0.01	?
		2	4010	0.8	?	0.8	0.8
		2	4010	0.7	?	0.7	0.7
		2	5010	0.9	0.9	0.9	0.9
		2	10000	?	0.9	?	?
		2	10000	0.9	?	0.9	0.9
		3	876	1	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {amode cond allways} {
	exec cg select -f {chromosome begin 
		{test=amode(1,$freq)} freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	0.2,0.11,0.1	0.1	0.11	0.2
		1	4001	0.2	0.2	0.1	0.2
		1	4050	?	0.3	?	?
		1	5000	0.6	0.4	0.6	0.6
		1	5020	?	0.5	?	?
		1	5020	0.4,0.5,?	?	0.4	0.5
		2	4000	0.6	0.6	0.6	0.6
		2	4001	?	0.8	?	?
		2	4001	0.7,?,0.01	0.7	0.01	?
		2	4010	0.8	?	0.8	0.8
		2	4010	0.7	?	0.7	0.7
		2	5010	0.9	0.9	0.9	0.9
		2	10000	?	0.9	?	?
		2	10000	0.9	?	0.9	0.9
		3	876	1	1	1	1
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {apercent} {
	exec cg select -f {chromosome begin 
		{test=apercent($sequenced == "v", lmax($freq) > 0.1)} sequenced-* freq-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	sequenced-gatk-sample1	sequenced-sam-sample1	sequenced-gatk-sample2	freq-gatk-sample1	freq-sam-sample1	freq-gatk-sample2
		1	259	66.66666666666667	v	v	v	0.1	0.11	0.2
		1	4001	50.0	v	v	r	0.2	0.1	0.2
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

test select {acount analysis} {
	exec cg select -overwrite 1 -f {chromosome begin 
		{test=acount(($analysis eq "gatk-sample1" or $analysis eq "sam-sample1") and $sequenced == "v")} sequenced-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	sequenced-gatk-sample1	sequenced-sam-sample1	sequenced-gatk-sample2
		1	259	2	v	v	v
		1	4001	2	v	v	r
		1	4050	1	v	u	u
		1	5000	2	v	v	r
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

test select {acount analysis in query} {
	exec cg select \
		-q {acount(($analysis eq "gatk-sample1" or $analysis eq "sam-sample1") and $sequenced == "v") > 1} \
		-f {chromosome begin sequenced-*} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	sequenced-gatk-sample1	sequenced-sam-sample1	sequenced-gatk-sample2
		1	259	v	v	v
		1	4001	v	v	r
		1	5000	v	v	r
		2	4000	v	v	v
		2	5010	v	v	v
		3	876	v	v	v
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {acount analysis in query and fields} {
	exec cg select -overwrite 1 \
		-q {acount(($analysis eq "gatk-sample1" or $analysis eq "sam-sample1") and $sequenced == "v") > 1} \
		-f {chromosome begin {test=acount(($analysis eq "gatk-sample1" or $analysis eq "sam-sample1") and $sequenced == "v")} sequenced-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	sequenced-gatk-sample1	sequenced-sam-sample1	sequenced-gatk-sample2
		1	259	2	v	v	v
		1	4001	2	v	v	r
		1	5000	2	v	v	r
		2	4000	2	v	v	v
		2	5010	2	v	v	v
		3	876	2	v	v	v
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {acount sample in query and fields} {
	exec cg select -overwrite 1 \
		-q {acount(($sample eq "sample1") and $sequenced == "v") > 1} \
		-f {chromosome begin {test=acount(($sample eq "sample1") and $sequenced == "v")} sequenced-*
	} data/vars-saggr.tsv tmp/results.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	test	sequenced-gatk-sample1	sequenced-sam-sample1	sequenced-gatk-sample2
		1	259	2	v	v	v
		1	4001	2	v	v	r
		1	5000	2	v	v	r
		2	4000	2	v	v	v
		2	5010	2	v	v	v
		3	876	2	v	v	v
	}
	exec diff tmp/results.tsv tmp/expected.tsv
} {}

test select {incorrect use of sample outside sample aggregates} {
	exec cg select \
		-q {$sample eq "gatk-sample1"} \
		-f {chromosome begin sequenced-*
	} data/vars-saggr.tsv
} {field sample not present in file (or sampleinfo)} error

test select {acount field not present in all analyses} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1	sequenced-sam-sample1 zyg-gatk-sample2
		chr1 4001 4002 snp A G,C	v	v	t
	}
	exec cg select -f {chromosome begin {test=acount($sequenced eq "v")}} tmp/vars.tsv
} {chromosome	begin	test
chr1	4001	2}

test select {acount field not present at all} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1	sequenced-sam-sample1 zyg-gatk-sample2
		chr1 4001 4002 snp A G,C	v	v	t
	}
	exec cg select -f {chromosome begin {test=acount($test eq "v")}} tmp/vars.tsv
} {error in acount: all analyses are missing one or more needed fields (test-gatk-sample1,test-sam-sample1,test-gatk-sample2)} error

test select {acount used in other function field not present in all analyses} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1	sequenced-sam-sample1 zyg-gatk-sample2
		chr1 4001 4002 snp A G,C	v	v	t
	}
	exec cg select -f {chromosome begin {test=if(acount($sequenced eq "v") == 2, acount($sequenced eq "v")+20, acount($sequenced eq "v")+10)}} tmp/vars.tsv
} {chromosome	begin	test
chr1	4001	22}

test select {acount used in other function field not present in all analyses, but in calccols} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1	sequenced-sam-sample1 zyg-gatk-sample2
		chr1 4001 4002 snp A G,C	v	v	t
	}
	exec cg select -f {chromosome begin -sequenced-gatk-sample2="v" {test=if(acount($sequenced eq "v") == 3, acount($sequenced eq "v")+20, acount($sequenced eq "v")+10)}} tmp/vars.tsv
} {chromosome	begin	test
chr1	4001	23}

test select {acount field not present in all analyses, but in calccol} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1	sequenced-sam-sample1 zyg-gatk-sample2
		chr1 4001 4002 snp A G,C	v	v	t
	}
	exec cg select -stack 1 -f {chromosome begin -sequenced-gatk-sample2="v" {test=acount($sequenced eq "v")}} tmp/vars.tsv
} {chromosome	begin	test
chr1	4001	3}

test select {acount field not present in all analyses, but in sampleinfo} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1	sequenced-sam-sample1 zyg-gatk-sample2
		chr1 4001 4002 snp A G,C	v	v	t
	}
	write_tab tmp/vars.tsv.sampleinfo.tsv {
		sample	sequenced
		gatk-sample2	v
	}
	exec cg select -stack 1 -f {chromosome begin {test=acount($sequenced eq "v")}} tmp/vars.tsv
} {chromosome	begin	test
chr1	4001	3}

test select {apercent field not present in all analyses} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1	sequenced-sam-sample1 zyg-gatk-sample2
		chr1 4001 4002 snp A G,C	v	t	t
	}
	exec cg select -f {chromosome begin {test=apercent($sequenced eq "v")}} tmp/vars.tsv
} {chromosome	begin	test
chr1	4001	50.0}

test select {apercent field not present in all analyses, 2 cond} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1 num-gatk-sample1	sequenced-sam-sample1 num-sam-sample1	zyg-gatk-sample2
		chr1 4001 4002 snp A G,C	v 1 t 0 t
	}
	exec cg select -f {chromosome begin {test=apercent($sequenced eq "v",$num == 1)}} tmp/vars.tsv
} {chromosome	begin	test
chr1	4001	100.0}

test select {apercent field not present in all analyses, 2 cond different missing} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1 num-gatk-sample1	sequenced-sam-sample1	zyg-gatk-sample2 num-gatk-sample2
		chr1 4001 4002 snp A G,C	v 1 t t 0
	}
	exec cg select -f {chromosome begin {test=apercent($sequenced eq "v",$num == 1)}} tmp/vars.tsv
} {chromosome	begin	test
chr1	4001	100.0}

test select {apercent field not present at all} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1	sequenced-sam-sample1 zyg-gatk-sample2
		chr1 4001 4002 snp A G,C	v	t	t
	}
	exec cg select -f {chromosome begin {test=apercent($test eq "v")}} tmp/vars.tsv
} {error in apercent: all analyses are missing one or more needed fields (test-gatk-sample1,test-sam-sample1,test-gatk-sample2)} error

test select {asum field not present in all analyses} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1	sequenced-sam-sample1 zyg-gatk-sample2
		chr1 4001 4002 snp A G,C	v	v	t
	}
	exec cg select -f {chromosome begin {test=asum($sequenced eq "v")}} tmp/vars.tsv
} {chromosome	begin	test
chr1	4001	2.0}

test select {asum field not present in all analyses, 2 cond} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1 num-gatk-sample1	sequenced-sam-sample1 num-sam-sample1	zyg-gatk-sample2
		chr1 4001 4002 snp A G,C	v 1 v 2 t
	}
	exec cg select -f {chromosome begin {test=asum($sequenced eq "v",$num)}} tmp/vars.tsv
} {chromosome	begin	test
chr1	4001	3.0}

test select {asum field not present in all analyses, 2 cond different missing} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1 num-gatk-sample1	sequenced-sam-sample1	zyg-gatk-sample2 num-gatk-sample2
		chr1 4001 4002 snp A G,C	v 1 v t 0
	}
	exec cg select -f {chromosome begin {test=asum($sequenced eq "v",$num)}} tmp/vars.tsv
} {chromosome	begin	test
chr1	4001	1.0}

test select {asum field not present in all analyses, 2 cond different missing but in calc and sampleinfo} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1 num-gatk-sample1	sequenced-sam-sample1	zyg-gatk-sample2 num-gatk-sample2
		chr1 4001 4002 snp A G,C	v 1 v t 3
	}
	write_tab tmp/vars.tsv.sampleinfo.tsv {
		sample	num
		sam-sample1	2
	}
	exec cg select -f {chromosome begin -sequenced-gatk-sample2="v" {test=asum($sequenced eq "v",$num)}} tmp/vars.tsv
} {chromosome	begin	test
chr1	4001	6.0}

test select {asum field not present at all} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1 num-gatk-sample1	sequenced-sam-sample1 num-sam-sample1	zyg-gatk-sample2
		chr1 4001 4002 snp A G,C	v 1 t 2 t
	}
	exec cg select -f {chromosome begin {test=asum($test eq "v",$num)}} tmp/vars.tsv
} {error in asum: all analyses are missing one or more needed fields (test-gatk-sample1,test-sam-sample1,test-gatk-sample2)} error

test select {asum field not present at all} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1 num-gatk-sample1	sequenced-sam-sample1 num-sam-sample1	zyg-gatk-sample2
		chr1 4001 4002 snp A G,C	v 1 t 2 t
	}
	exec cg select -f {chromosome begin {test=asum($sequenced eq "v",$test)}} tmp/vars.tsv
} {error in asum: all analyses are missing one or more needed fields (test-gatk-sample1,test-sam-sample1,sequenced-gatk-sample2)} error

test select {use analysis in acount condition} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	zyg-gatk-samplea1	zyg-gatk-samplea2	zyg-sam-samplea2	zyg-gatk-sampleb3
		chr1	4001	4002	snp	A	G	t	r	t	t
	}
	exec cg select -f {chromosome begin {count=acount($zyg == "t" and $analysis matches "gatk-samplea*")}} tmp/vars.tsv
} {chromosome	begin	count
chr1	4001	1}

test select {use sample in acount condition} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	zyg-gatk-samplea1	zyg-gatk-samplea2	zyg-sam-samplea2	zyg-gatk-sampleb3
		chr1	4001	4002	snp	A	G	t	r	t	t
	}
	exec cg select -f {chromosome begin {count=acount($zyg == "t" and $sample matches "samplea*")}} tmp/vars.tsv
} {chromosome	begin	count
chr1	4001	2}

test select {use analysis in alist condition} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	zyg-gatk-samplea1	zyg-gatk-samplea2	zyg-sam-samplea2	zyg-gatk-sampleb3
		chr1	4001	4002	snp	A	G	t	r	t	t
	}
	exec cg select -f {chromosome begin {list=alist($analysis matches "gatk-samplea*",$zyg)}} tmp/vars.tsv
} {chromosome	begin	list
chr1	4001	t,r}

test select {use analysis in acount in query} {
	write_tab tmp/vars.tsv {
		chromosome begin end type ref alt	sequenced-gatk-sample1	coverage-gatk-sample1	sequenced-sam-sample1	coverage-sam-sample1 sequenced-gatk-sample2 coverage-gatk-sample2
		chr1	4001	4002	snp	A	G	v	25	r	25	v	10
		chr1	4002	4003	snp	A	G	v	10	v	25	v	10
	}
	cg select -stack 1 -q {
	    $type=="snp" and acount($sequenced == "v" 
	    and $coverage >= 20 and $analysis matches "gatk-*") >= 1
	} -g all tmp/vars.tsv
} {all	count
all	1}

testsummarize
