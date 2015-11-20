#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test select {scount basic} {
	csv_parse [exec cg select -f {chromosome begin 
		{test=scount($sequenced == "v")} sequenced-*
	} data/vars-saggr.tsv] \t
} {{chromosome begin test sequenced-sample1 sequenced-sample2 sequenced-sample3} {1 259 3 v v v} {1 4001 2 v v r} {1 4050 1 v u u} {1 5000 2 v v r} {1 5020 1 v r r} {1 5020 2 r v v} {2 4000 3 v v v} {2 4001 1 v r r} {2 4001 1 v r r} {2 4010 2 u v v} {2 4010 2 u v v} {2 5010 3 v v v} {2 10000 1 v r r} {2 10000 2 r v v} {3 876 3 v v v}}

test select {scount and} {
	csv_parse [exec cg select -f {chromosome begin
		{test=scount($sequenced == "v" and $freq > 0.5)} freq-*
	} data/vars-saggr.tsv] \t
} {{chromosome begin test freq-sample1 freq-sample2 freq-sample3} {1 259 0 0.1 0.11 0.2} {1 4001 0 0.2 0.1 0.2} {1 4050 0 0.3 ? ?} {1 5000 1 0.4 0.6 0.6} {1 5020 0 0.5 ? ?} {1 5020 0 ? 0.4 0.5} {2 4000 3 0.6 0.6 0.6} {2 4001 1 0.8 ? ?} {2 4001 1 0.7 0.01 ?} {2 4010 2 ? 0.8 0.8} {2 4010 2 ? 0.7 0.7} {2 5010 3 0.9 0.9 0.9} {2 10000 1 0.9 ? ?} {2 10000 2 ? 0.9 0.9} {3 876 3 1 1 1}}

test select {slist basic} {
	csv_parse [exec cg select -f {chromosome begin
		{test=slist($sequenced)} sequenced-*
	} data/vars-saggr.tsv] \t
} {{chromosome begin test sequenced-sample1 sequenced-sample2 sequenced-sample3} {1 259 v,v,v v v v} {1 4001 v,v,r v v r} {1 4050 v,u,u v u u} {1 5000 v,v,r v v r} {1 5020 v,r,r v r r} {1 5020 r,v,v r v v} {2 4000 v,v,v v v v} {2 4001 v,r,r v r r} {2 4001 v,r,r v r r} {2 4010 u,v,v u v v} {2 4010 u,v,v u v v} {2 5010 v,v,v v v v} {2 10000 v,r,r v r r} {2 10000 r,v,v r v v} {3 876 v,v,v v v v}}

test select {slist nested} {
	csv_parse [exec cg select -f {chromosome begin
		{test=slist(if($sequenced != "v","u",$freq))} freq-*
	} data/vars-saggr.tsv] \t
} {{chromosome begin test freq-sample1 freq-sample2 freq-sample3} {1 259 0.1,0.11,0.2 0.1 0.11 0.2} {1 4001 0.2,0.1,u 0.2 0.1 0.2} {1 4050 0.3,u,u 0.3 ? ?} {1 5000 0.4,0.6,u 0.4 0.6 0.6} {1 5020 0.5,u,u 0.5 ? ?} {1 5020 u,0.4,0.5 ? 0.4 0.5} {2 4000 0.6,0.6,0.6 0.6 0.6 0.6} {2 4001 0.8,u,u 0.8 ? ?} {2 4001 0.7,u,u 0.7 0.01 ?} {2 4010 u,0.8,0.8 ? 0.8 0.8} {2 4010 u,0.7,0.7 ? 0.7 0.7} {2 5010 0.9,0.9,0.9 0.9 0.9 0.9} {2 10000 0.9,u,u 0.9 ? ?} {2 10000 u,0.9,0.9 ? 0.9 0.9} {3 876 1,1,1 1 1 1}}

test select {slist condition} {
	csv_parse [exec cg select -f {chromosome begin
		{test=slist($sequenced == "v" and $freq > 0.5,if($sequenced != "v","u",$freq))} freq-*
	} data/vars-saggr.tsv] \t
} {{chromosome begin test freq-sample1 freq-sample2 freq-sample3} {1 259 {} 0.1 0.11 0.2} {1 4001 {} 0.2 0.1 0.2} {1 4050 {} 0.3 ? ?} {1 5000 0.6 0.4 0.6 0.6} {1 5020 {} 0.5 ? ?} {1 5020 {} ? 0.4 0.5} {2 4000 0.6,0.6,0.6 0.6 0.6 0.6} {2 4001 0.8 0.8 ? ?} {2 4001 0.7 0.7 0.01 ?} {2 4010 0.8,0.8 ? 0.8 0.8} {2 4010 0.7,0.7 ? 0.7 0.7} {2 5010 0.9,0.9,0.9 0.9 0.9 0.9} {2 10000 0.9 0.9 ? ?} {2 10000 0.9,0.9 ? 0.9 0.9} {3 876 1,1,1 1 1 1}}

test select {sdistinct} {
	csv_parse [exec cg select -f {chromosome begin 
		{test=sdistinct(if($sequenced != "v","u",$freq))} freq-*
	} data/vars-saggr.tsv] \t
} {{chromosome begin test freq-sample1 freq-sample2 freq-sample3} {1 259 0.1,0.11,0.2 0.1 0.11 0.2} {1 4001 0.2,0.1,u 0.2 0.1 0.2} {1 4050 0.3,u 0.3 ? ?} {1 5000 0.4,0.6,u 0.4 0.6 0.6} {1 5020 0.5,u 0.5 ? ?} {1 5020 u,0.4,0.5 ? 0.4 0.5} {2 4000 0.6 0.6 0.6 0.6} {2 4001 0.8,u 0.8 ? ?} {2 4001 0.7,u 0.7 0.01 ?} {2 4010 u,0.8 ? 0.8 0.8} {2 4010 u,0.7 ? 0.7 0.7} {2 5010 0.9 0.9 0.9 0.9} {2 10000 0.9,u 0.9 ? ?} {2 10000 u,0.9 ? 0.9 0.9} {3 876 1 1 1 1}}

test select {sucount} {
	csv_parse [exec cg select -f {chromosome begin 
		{test=sucount(if($sequenced != "v","u",$freq))} freq-*
	} data/vars-saggr.tsv] \t
} {{chromosome begin test freq-sample1 freq-sample2 freq-sample3} {1 259 3 0.1 0.11 0.2} {1 4001 3 0.2 0.1 0.2} {1 4050 2 0.3 ? ?} {1 5000 3 0.4 0.6 0.6} {1 5020 2 0.5 ? ?} {1 5020 3 ? 0.4 0.5} {2 4000 1 0.6 0.6 0.6} {2 4001 2 0.8 ? ?} {2 4001 2 0.7 0.01 ?} {2 4010 2 ? 0.8 0.8} {2 4010 2 ? 0.7 0.7} {2 5010 1 0.9 0.9 0.9} {2 10000 2 0.9 ? ?} {2 10000 2 ? 0.9 0.9} {3 876 1 1 1 1}}

test select {sdistinct condition} {
	csv_parse [exec cg select -f {chromosome begin 
		{test=sdistinct($sequenced == "v" and $freq > 0.5,if($sequenced != "v","u",$freq))} freq-*
	} data/vars-saggr.tsv] \t
} {{chromosome begin test freq-sample1 freq-sample2 freq-sample3} {1 259 {} 0.1 0.11 0.2} {1 4001 {} 0.2 0.1 0.2} {1 4050 {} 0.3 ? ?} {1 5000 0.6 0.4 0.6 0.6} {1 5020 {} 0.5 ? ?} {1 5020 {} ? 0.4 0.5} {2 4000 0.6 0.6 0.6 0.6} {2 4001 0.8 0.8 ? ?} {2 4001 0.7 0.7 0.01 ?} {2 4010 0.8 ? 0.8 0.8} {2 4010 0.7 ? 0.7 0.7} {2 5010 0.9 0.9 0.9 0.9} {2 10000 0.9 0.9 ? ?} {2 10000 0.9 ? 0.9 0.9} {3 876 1 1 1 1}}

test select {sdistinct condition} {
	csv_parse [exec cg select -f {chromosome begin 
		{test=sdistinct($alleleSeq2 == "N",$sequenced)} alleleSeq2-* sequenced-*
	} data/vars-saggr.tsv] \t
} {{chromosome begin test alleleSeq2-sample1 alleleSeq2-sample2 alleleSeq2-sample3 sequenced-sample1 sequenced-sample2 sequenced-sample3} {1 259 {} C C C v v v} {1 4001 {} C C C v v r} {1 4050 {} T - - v u u} {1 5000 {} {} {} {} v v r} {1 5020 {} A G G v r r} {1 5020 {} G C C r v v} {2 4000 {} A A A v v v} {2 4001 {} C A A v r r} {2 4001 {} C A A v r r} {2 4010 {} - C C u v v} {2 4010 {} - C C u v v} {2 5010 v {} N N v v v} {2 10000 r,v N N N v r r} {2 10000 v,r N N N r v v} {3 876 {} A A A v v v}}

test select {smin} {
	csv_parse [exec cg select -f {chromosome begin 
		{test=smin($freq)} freq-*
	} data/vars-saggr.tsv] \t
} {{chromosome begin test freq-sample1 freq-sample2 freq-sample3} {1 259 0.1 0.1 0.11 0.2} {1 4001 0.1 0.2 0.1 0.2} {1 4050 0.3 0.3 ? ?} {1 5000 0.4 0.4 0.6 0.6} {1 5020 0.5 0.5 ? ?} {1 5020 0.4 ? 0.4 0.5} {2 4000 0.6 0.6 0.6 0.6} {2 4001 0.8 0.8 ? ?} {2 4001 0.01 0.7 0.01 ?} {2 4010 0.8 ? 0.8 0.8} {2 4010 0.7 ? 0.7 0.7} {2 5010 0.9 0.9 0.9 0.9} {2 10000 0.9 0.9 ? ?} {2 10000 0.9 ? 0.9 0.9} {3 876 1 1 1 1}}

test select {smin condition} {
	csv_parse [exec cg select -f {chromosome begin 
		{test=smin($sequenced == "v", $freq)} sequenced-* freq-*
	} data/vars-saggr.tsv] \t
} {{chromosome begin test sequenced-sample1 sequenced-sample2 sequenced-sample3 freq-sample1 freq-sample2 freq-sample3} {1 259 0.1 v v v 0.1 0.11 0.2} {1 4001 0.1 v v r 0.2 0.1 0.2} {1 4050 0.3 v u u 0.3 ? ?} {1 5000 0.4 v v r 0.4 0.6 0.6} {1 5020 0.5 v r r 0.5 ? ?} {1 5020 0.4 r v v ? 0.4 0.5} {2 4000 0.6 v v v 0.6 0.6 0.6} {2 4001 0.8 v r r 0.8 ? ?} {2 4001 0.7 v r r 0.7 0.01 ?} {2 4010 0.8 u v v ? 0.8 0.8} {2 4010 0.7 u v v ? 0.7 0.7} {2 5010 0.9 v v v 0.9 0.9 0.9} {2 10000 0.9 v r r 0.9 ? ?} {2 10000 0.9 r v v ? 0.9 0.9} {3 876 1 v v v 1 1 1}}

test select {smax} {
	csv_parse [exec cg select -f {chromosome begin 
		{test=smax($freq)} freq-*
	} data/vars-saggr.tsv] \t
} {{chromosome begin test freq-sample1 freq-sample2 freq-sample3} {1 259 0.2 0.1 0.11 0.2} {1 4001 0.2 0.2 0.1 0.2} {1 4050 0.3 0.3 ? ?} {1 5000 0.6 0.4 0.6 0.6} {1 5020 0.5 0.5 ? ?} {1 5020 0.5 ? 0.4 0.5} {2 4000 0.6 0.6 0.6 0.6} {2 4001 0.8 0.8 ? ?} {2 4001 0.7 0.7 0.01 ?} {2 4010 0.8 ? 0.8 0.8} {2 4010 0.7 ? 0.7 0.7} {2 5010 0.9 0.9 0.9 0.9} {2 10000 0.9 0.9 ? ?} {2 10000 0.9 ? 0.9 0.9} {3 876 1 1 1 1}}

test select {smax condition} {
	csv_parse [exec cg select -f {chromosome begin 
		{test=smax(lmin($freq) < 0.2, $freq)} sequenced-* freq-*
	} data/vars-saggr.tsv] \t
} {{chromosome begin test sequenced-sample1 sequenced-sample2 sequenced-sample3 freq-sample1 freq-sample2 freq-sample3} {1 259 0.11 v v v 0.1 0.11 0.2} {1 4001 0.1 v v r 0.2 0.1 0.2} {1 4050 NaN v u u 0.3 ? ?} {1 5000 NaN v v r 0.4 0.6 0.6} {1 5020 NaN v r r 0.5 ? ?} {1 5020 NaN r v v ? 0.4 0.5} {2 4000 NaN v v v 0.6 0.6 0.6} {2 4001 NaN v r r 0.8 ? ?} {2 4001 0.01 v r r 0.7 0.01 ?} {2 4010 NaN u v v ? 0.8 0.8} {2 4010 NaN u v v ? 0.7 0.7} {2 5010 NaN v v v 0.9 0.9 0.9} {2 10000 NaN v r r 0.9 ? ?} {2 10000 NaN r v v ? 0.9 0.9} {3 876 NaN v v v 1 1 1}}

test select {ssum} {
	csv_parse [exec cg select -f {chromosome begin 
		{test=ssum($freq)} freq-*
	} data/vars-saggr.tsv] \t
} {{chromosome begin test freq-sample1 freq-sample2 freq-sample3} {1 259 0.41000000000000003 0.1 0.11 0.2} {1 4001 0.5 0.2 0.1 0.2} {1 4050 0.3 0.3 ? ?} {1 5000 1.6 0.4 0.6 0.6} {1 5020 0.5 0.5 ? ?} {1 5020 0.9 ? 0.4 0.5} {2 4000 1.7999999999999998 0.6 0.6 0.6} {2 4001 0.8 0.8 ? ?} {2 4001 0.71 0.7 0.01 ?} {2 4010 1.6 ? 0.8 0.8} {2 4010 1.4 ? 0.7 0.7} {2 5010 2.7 0.9 0.9 0.9} {2 10000 0.9 0.9 ? ?} {2 10000 1.8 ? 0.9 0.9} {3 876 3.0 1 1 1}}

test select {ssum condition} {
	csv_parse [exec cg select -f {chromosome begin 
		{test=ssum($sequenced == "v", $freq)} sequenced-* freq-*
	} data/vars-saggr.tsv] \t
} {{chromosome begin test sequenced-sample1 sequenced-sample2 sequenced-sample3 freq-sample1 freq-sample2 freq-sample3} {1 259 0.41000000000000003 v v v 0.1 0.11 0.2} {1 4001 0.30000000000000004 v v r 0.2 0.1 0.2} {1 4050 0.3 v u u 0.3 ? ?} {1 5000 1.0 v v r 0.4 0.6 0.6} {1 5020 0.5 v r r 0.5 ? ?} {1 5020 0.9 r v v ? 0.4 0.5} {2 4000 1.7999999999999998 v v v 0.6 0.6 0.6} {2 4001 0.8 v r r 0.8 ? ?} {2 4001 0.7 v r r 0.7 0.01 ?} {2 4010 1.6 u v v ? 0.8 0.8} {2 4010 1.4 u v v ? 0.7 0.7} {2 5010 2.7 v v v 0.9 0.9 0.9} {2 10000 0.9 v r r 0.9 ? ?} {2 10000 1.8 r v v ? 0.9 0.9} {3 876 3.0 v v v 1 1 1}}

test select {savg} {
	file_write tmp/temp1 [join [csv_parse [exec cg select -f {chromosome begin 
		{test=format("%.2f",savg($freq))} freq-*
	} data/vars-saggr.tsv] \t] \n]
	file_write tmp/temp2 [join {{chromosome begin test freq-sample1 freq-sample2 freq-sample3} {1 259 0.14 0.1 0.11 0.2} {1 4001 0.17 0.2 0.1 0.2} {1 4050 0.3 0.3 ? ?} {1 5000 0.53 0.4 0.6 0.6} {1 5020 0.5 0.5 ? ?} {1 5020 0.45 ? 0.4 0.5} {2 4000 0.6 0.6 0.6 0.6} {2 4001 0.8 0.8 ? ?} {2 4001 0.35 0.7 0.01 ?} {2 4010 0.8 ? 0.8 0.8} {2 4010 0.7 ? 0.7 0.7} {2 5010 0.9 0.9 0.9 0.9} {2 10000 0.9 0.9 ? ?} {2 10000 0.9 ? 0.9 0.9} {3 876 1.0 1 1 1}} \n]
	exec diff tmp/temp1 tmp/temp2
} {}

test select {savg cond allways} {
	file_write tmp/temp1 [join [csv_parse [exec cg select -f {chromosome begin 
		{test=format("%.2f",savg(1,$freq))} freq-*
	} data/vars-saggr.tsv] \t] \n]
	file_write tmp/temp2 [join {{chromosome begin test freq-sample1 freq-sample2 freq-sample3} {1 259 0.14 0.1 0.11 0.2} {1 4001 0.17 0.2 0.1 0.2} {1 4050 0.3 0.3 ? ?} {1 5000 0.53 0.4 0.6 0.6} {1 5020 0.5 0.5 ? ?} {1 5020 0.45 ? 0.4 0.5} {2 4000 0.6 0.6 0.6 0.6} {2 4001 0.8 0.8 ? ?} {2 4001 0.35 0.7 0.01 ?} {2 4010 0.8 ? 0.8 0.8} {2 4010 0.7 ? 0.7 0.7} {2 5010 0.9 0.9 0.9 0.9} {2 10000 0.9 0.9 ? ?} {2 10000 0.9 ? 0.9 0.9} {3 876 1.0 1 1 1}} \n]
	exec diff tmp/temp1 tmp/temp2
} {}

test select {sstdev} {
	csv_parse [exec cg select -f {chromosome begin 
		{test=sstdev($freq)} freq-*
	} data/vars-saggr.tsv] \t
} {{chromosome begin test freq-sample1 freq-sample2 freq-sample3} {1 259 0.044969125210773474 0.1 0.11 0.2} {1 4001 0.04714045207910317 0.2 0.1 0.2} {1 4050 0.0 0.3 ? ?} {1 5000 0.09428090415820632 0.4 0.6 0.6} {1 5020 0.0 0.5 ? ?} {1 5020 0.04999999999999999 ? 0.4 0.5} {2 4000 0.0 0.6 0.6 0.6} {2 4001 0.0 0.8 ? ?} {2 4001 0.345 0.7 0.01 ?} {2 4010 0.0 ? 0.8 0.8} {2 4010 0.0 ? 0.7 0.7} {2 5010 0.0 0.9 0.9 0.9} {2 10000 0.0 0.9 ? ?} {2 10000 0.0 ? 0.9 0.9} {3 876 0.0 1 1 1}}

test select {sstdev cond allways} {
	csv_parse [exec cg select -f {chromosome begin 
		{test=sstdev(1,$freq)} freq-*
	} data/vars-saggr.tsv] \t
} {{chromosome begin test freq-sample1 freq-sample2 freq-sample3} {1 259 0.044969125210773474 0.1 0.11 0.2} {1 4001 0.04714045207910317 0.2 0.1 0.2} {1 4050 0.0 0.3 ? ?} {1 5000 0.09428090415820632 0.4 0.6 0.6} {1 5020 0.0 0.5 ? ?} {1 5020 0.04999999999999999 ? 0.4 0.5} {2 4000 0.0 0.6 0.6 0.6} {2 4001 0.0 0.8 ? ?} {2 4001 0.345 0.7 0.01 ?} {2 4010 0.0 ? 0.8 0.8} {2 4010 0.0 ? 0.7 0.7} {2 5010 0.0 0.9 0.9 0.9} {2 10000 0.0 0.9 ? ?} {2 10000 0.0 ? 0.9 0.9} {3 876 0.0 1 1 1}}

test select {smedian} {
	csv_parse [exec cg select -f {chromosome begin 
		{test=smedian($freq)} freq-*
	} data/vars-saggr.tsv] \t
} {{chromosome begin test freq-sample1 freq-sample2 freq-sample3} {1 259 0.11 0.1 0.11 0.2} {1 4001 0.2 0.2 0.1 0.2} {1 4050 0.3 0.3 ? ?} {1 5000 0.6 0.4 0.6 0.6} {1 5020 0.5 0.5 ? ?} {1 5020 0.45 ? 0.4 0.5} {2 4000 0.6 0.6 0.6 0.6} {2 4001 0.8 0.8 ? ?} {2 4001 0.355 0.7 0.01 ?} {2 4010 0.8 ? 0.8 0.8} {2 4010 0.7 ? 0.7 0.7} {2 5010 0.9 0.9 0.9 0.9} {2 10000 0.9 0.9 ? ?} {2 10000 0.9 ? 0.9 0.9} {3 876 1 1 1 1}}

test select {smedian cond allways} {
	csv_parse [exec cg select -f {chromosome begin 
		{test=smedian(1,$freq)} freq-*
	} data/vars-saggr.tsv] \t
} {{chromosome begin test freq-sample1 freq-sample2 freq-sample3} {1 259 0.11 0.1 0.11 0.2} {1 4001 0.2 0.2 0.1 0.2} {1 4050 0.3 0.3 ? ?} {1 5000 0.6 0.4 0.6 0.6} {1 5020 0.5 0.5 ? ?} {1 5020 0.45 ? 0.4 0.5} {2 4000 0.6 0.6 0.6 0.6} {2 4001 0.8 0.8 ? ?} {2 4001 0.355 0.7 0.01 ?} {2 4010 0.8 ? 0.8 0.8} {2 4010 0.7 ? 0.7 0.7} {2 5010 0.9 0.9 0.9 0.9} {2 10000 0.9 0.9 ? ?} {2 10000 0.9 ? 0.9 0.9} {3 876 1 1 1 1}}

test select {smode} {
	csv_parse [exec cg select -f {chromosome begin 
		{test=smode($freq)} freq-*
	} data/vars-saggr.tsv] \t
} {{chromosome begin test freq-sample1 freq-sample2 freq-sample3} {1 259 0.2,0.11,0.1 0.1 0.11 0.2} {1 4001 0.2 0.2 0.1 0.2} {1 4050 ? 0.3 ? ?} {1 5000 0.6 0.4 0.6 0.6} {1 5020 ? 0.5 ? ?} {1 5020 0.4,0.5,? ? 0.4 0.5} {2 4000 0.6 0.6 0.6 0.6} {2 4001 ? 0.8 ? ?} {2 4001 0.7,?,0.01 0.7 0.01 ?} {2 4010 0.8 ? 0.8 0.8} {2 4010 0.7 ? 0.7 0.7} {2 5010 0.9 0.9 0.9 0.9} {2 10000 ? 0.9 ? ?} {2 10000 0.9 ? 0.9 0.9} {3 876 1 1 1 1}}

test select {smode cond allways} {
	csv_parse [exec cg select -f {chromosome begin 
		{test=smode(1,$freq)} freq-*
	} data/vars-saggr.tsv] \t
} {{chromosome begin test freq-sample1 freq-sample2 freq-sample3} {1 259 0.2,0.11,0.1 0.1 0.11 0.2} {1 4001 0.2 0.2 0.1 0.2} {1 4050 ? 0.3 ? ?} {1 5000 0.6 0.4 0.6 0.6} {1 5020 ? 0.5 ? ?} {1 5020 0.4,0.5,? ? 0.4 0.5} {2 4000 0.6 0.6 0.6 0.6} {2 4001 ? 0.8 ? ?} {2 4001 0.7,?,0.01 0.7 0.01 ?} {2 4010 0.8 ? 0.8 0.8} {2 4010 0.7 ? 0.7 0.7} {2 5010 0.9 0.9 0.9 0.9} {2 10000 ? 0.9 ? ?} {2 10000 0.9 ? 0.9 0.9} {3 876 1 1 1 1}}

test select {spercent} {
	csv_parse [exec cg select -f {chromosome begin 
		{test=spercent($sequenced == "v", lmax($freq) > 0.1)} sequenced-* freq-*
	} data/vars-saggr.tsv] \t
} {{chromosome begin test sequenced-sample1 sequenced-sample2 sequenced-sample3 freq-sample1 freq-sample2 freq-sample3} {1 259 66.66666666666667 v v v 0.1 0.11 0.2} {1 4001 50.0 v v r 0.2 0.1 0.2} {1 4050 100.0 v u u 0.3 ? ?} {1 5000 100.0 v v r 0.4 0.6 0.6} {1 5020 100.0 v r r 0.5 ? ?} {1 5020 100.0 r v v ? 0.4 0.5} {2 4000 100.0 v v v 0.6 0.6 0.6} {2 4001 100.0 v r r 0.8 ? ?} {2 4001 100.0 v r r 0.7 0.01 ?} {2 4010 100.0 u v v ? 0.8 0.8} {2 4010 100.0 u v v ? 0.7 0.7} {2 5010 100.0 v v v 0.9 0.9 0.9} {2 10000 100.0 v r r 0.9 ? ?} {2 10000 100.0 r v v ? 0.9 0.9} {3 876 100.0 v v v 1 1 1}}

test select {scount sample} {
	csv_parse [exec cg select -f {chromosome begin 
		{test=scount(($sample eq "sample1" or $sample eq "sample2") and $sequenced == "v")} sequenced-*
	} data/vars-saggr.tsv] \t
} {{chromosome begin test sequenced-sample1 sequenced-sample2 sequenced-sample3} {1 259 2 v v v} {1 4001 2 v v r} {1 4050 1 v u u} {1 5000 2 v v r} {1 5020 1 v r r} {1 5020 1 r v v} {2 4000 2 v v v} {2 4001 1 v r r} {2 4001 1 v r r} {2 4010 1 u v v} {2 4010 1 u v v} {2 5010 2 v v v} {2 10000 1 v r r} {2 10000 1 r v v} {3 876 2 v v v}}

test select {scount sample in query} {
	exec cg select \
		-q {scount(($sample eq "sample1" or $sample eq "sample2") and $sequenced == "v") > 1} \
		-f {chromosome begin sequenced-*} data/vars-saggr.tsv
} {chromosome	begin	sequenced-sample1	sequenced-sample2	sequenced-sample3
1	259	v	v	v
1	4001	v	v	r
1	5000	v	v	r
2	4000	v	v	v
2	5010	v	v	v
3	876	v	v	v}

test select {scount sample in query and fields} {
	exec cg select \
		-q {scount(($sample eq "sample1" or $sample eq "sample2") and $sequenced == "v") > 1} \
		-f {chromosome begin {test=scount(($sample eq "sample1" or $sample eq "sample2") and $sequenced == "v")} sequenced-*
	} data/vars-saggr.tsv
} {chromosome	begin	test	sequenced-sample1	sequenced-sample2	sequenced-sample3
1	259	2	v	v	v
1	4001	2	v	v	r
1	5000	2	v	v	r
2	4000	2	v	v	v
2	5010	2	v	v	v
3	876	2	v	v	v}

test select {incorrect use of sample outside sample aggregates} {
	exec cg select \
		-q {$sample eq "sample1"} \
		-f {chromosome begin sequenced-*
	} data/vars-saggr.tsv
} {field "sample" not present in file and no sampleinfo file found} error

testsummarize
