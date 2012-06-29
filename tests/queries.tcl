#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl
cd /complgen/testqueries

test query {query 1} {
	exec cg select -q {
		$refGene_impact ~/CDS|UTR/ || $knownGene_impact ~/CDS|UTR/ || $ensGene_impact ~/CDS|UTR/
		|| $acembly_impact ~/CDS|UTR/ || $genscan_impact ~/CDS|UTR/
	} annottestcompar.tsv.rz temp.tsv
	exec diff temp.tsv answers/query1.tsv
} {}

test query {query 2} {
	exec cg select -q {
		count($*_impact, ~/CDS|UTR/) > 0
	} annottestcompar.tsv.rz temp.tsv
	exec diff temp.tsv answers/query1.tsv
} {}

test query {query 3} {
	exec cg select -q {$type=="ins"&&$coverage-testNA19240chr2122cg>=20} annottestcompar.tsv.rz temp.tsv
	exec diff temp.tsv answers/query3.tsv
} {}

test query {query 4} {
	exec cg select -q { $type=="snp" && $impact == "MISSENSE" } \
		-f {chromosome begin end} \
		annottestcompar.tsv.rz temp.tsv
	exec diff temp.tsv answers/query4.tsv
} {}

test query {query 5} {
	exec cg select -h	annottestcompar.tsv.rz > temp.tsv
	exec diff temp.tsv answers/query5.tsv
} {}

test query {query 6 region} {
	exec cg select -q {
		$chromosome == "22" && $begin < 16434940 && $end > 16054739
	} annottestcompar.tsv.rz temp.tsv
	exec diff temp.tsv answers/query6.tsv
} {}

test query {query 6 region 2} {
	exec cg select -q {
		region("22:16054739-16434940")
	} annottestcompar.tsv.rz temp.tsv
	exec diff temp.tsv answers/query6.tsv
} {}

test query {query 6 region 3} {
	exec cg select -q {
		region(22:16054739:16434940)
	} annottestcompar.tsv.rz temp.tsv
	exec diff temp.tsv answers/query6.tsv
} {}

test query {query 6 region 4} {
	exec cg select -q {
		region("chr22",16054739,16434940)
	} annottestcompar.tsv.rz temp.tsv
	exec diff temp.tsv answers/query6.tsv
} {}

test query {query regions} {
	exec cg select -q {
		region("chr22:16054739-16434940",21,9415561,9417128)
	} annottestcompar.tsv.rz temp.tsv
	exec diff temp.tsv answers/query6.tsv
} {2,4d1
< 21	9415561	9415562	del	G		?	-	-	?	?	-496	143			u	20901814		G	55	55	-663	167			v																																1	p11.2				1													MER41D											-	-	-	-	-	-	-	-	-	-	-	-	-	-
< 21	9417126	9417127	snp	C	A	18820400	A	C	96	96	-319	101			v	?	C	C	?	?	-192	172			r																																1	p11.2				1													HAL1											-	-	-	-	-	-	-	-	-	-	-	-	-	-
< 21	9417127	9417128	snp	G	A	18820400	A	G	96	96	-410	101			v	20901846	A	G	137	137	-298	170			v	dbsnp.92:rs1850103																															1	p11.2				1													HAL1											-	-	-	-	-	-	-	-	-	-	-	-	rs1850103	0.5
child process exited abnormally} error

test query {query 7} {
	exec cg select -q {
		sm(testNA19239chr2122cg, testNA19240chr2122cg)
	} annottestcompar.tsv.rz temp.tsv
	exec diff temp.tsv answers/query7.tsv
} {}

test query {query 8} {
	exec cg select -q {
		same(testNA19239chr2122cg, testNA19240chr2122cg)
	} annottestcompar.tsv.rz temp.tsv
	exec diff temp.tsv answers/query8.tsv
} {}

test query {query 9} {
	exec cg select -q {
		mm(testNA19239chr2122cg, testNA19240chr2122cg)
	} annottestcompar.tsv.rz temp.tsv
	exec diff temp.tsv answers/query9.tsv
} {}

test query {query 10} {
	exec cg select -q {
		$type == "snp" 
		&& $simpleRepeat == "" && $microsat == "" && $genomicSuperDups == "" 
	} annottestcompar.tsv.rz temp.tsv
	exec diff temp.tsv answers/query10.tsv
} {}

test query {query 11} {
	exec cg select -q {
		$sequenced-testNA19239chr2122cg == "v"
		&& $sequenced-testNA19240chr2122cg == "v"
	} annottestcompar.tsv.rz temp.tsv
	exec diff temp.tsv answers/query11.tsv
} {}

test query {query 12} {
	exec cg select -q {
		count($sequenced-*, == "v")  == 2
	} annottestcompar.tsv.rz temp.tsv
	exec diff temp.tsv answers/query12.tsv
} {}

test query {query 13} {
	exec cg select -q {
		$sequenced-testNA19239chr2122cg == "v"
		&& $type == "snp"
		&& $coverage-testNA19239chr2122cg >= 20 
		&& $coverage-testNA19239chr2122cg <= 100
		&& $cluster-testNA19239chr2122cg == ""
	} annottestcompar.tsv.rz temp.tsv
	exec diff temp.tsv answers/query13.tsv
} {}

test query {query 14} {
	exec cg select -q {
		count($sequenced-*, == "v")  == 4
		&& $type == "snp"
		&& count($coverage-*, <20) == 0
		&& count($coverage-*, > 100) == 0
		&& count($cluster-*, != "") == 0
	} annottestcompar.tsv.rz temp.tsv
	exec diff temp.tsv answers/query14.tsv
} {}

test query {query 15} {
	exec cg select -q {
		count($*_impact, ~/CDS|UTR/) > 0
	} annottestcompar.tsv.rz temp.tsv
	exec diff temp.tsv answers/query15.tsv
} {}

test query {query 16} {
	exec cg select -q {
		$ensGene_impact ~ /CDS/
	} annottestcompar.tsv.rz temp.tsv
	exec diff temp.tsv answers/query16.tsv
} {}

test query {query 17} {
	exec cg select -q {
		counthasone($ensGene_impact, == "CDSMIS") > 0
	} annottestcompar.tsv.rz temp.tsv
	exec diff temp.tsv answers/query17.tsv
} {}

test query {query 18} {
	exec cg select -q {
		counthasone($*_impact, ~ /CDS/) > 0
	} annottestcompar.tsv.rz temp.tsv
	exec diff temp.tsv answers/query18.tsv
} {}

test query {query 19} {
	exec cg select -q {
		counthasall($*_impact, ~ /CDS/) > 0
	} annottestcompar.tsv.rz temp.tsv
	exec diff temp.tsv answers/query19.tsv
} {}

test query {query 20} {
	exec cg select -f "chromosome begin end refGene_impact" -q {
		oneof($refGene_impact, "CDSMIS", "CDSINS") > 0
	} annottestcompar.tsv.rz temp.tsv
	exec diff temp.tsv answers/query20.tsv
} {}

test query {query 20 2} {
	exec cg select -f "chromosome begin end refGene_impact" -q {
		$refGene_impact in {"CDSMIS" "CDSINS"}
	} annottestcompar.tsv.rz temp.tsv
	exec diff temp.tsv answers/query20.tsv
} {}

test query {query 21 combine} {
	exec cg select -f "chromosome begin end 1000gCEU 1000glow" -q {
		count(lmax($1000gCEU),lmax($1000glow), > 0.98) > 0
	} annottestcompar.tsv.rz temp.tsv
	exec diff temp.tsv answers/query21.tsv
} {}

testsummarize
