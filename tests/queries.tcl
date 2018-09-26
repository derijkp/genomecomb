#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl
set keepdir [pwd]

test query {query 1} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -q {
		$refGene_impact ~/CDS|UTR/ || $knownGene_impact ~/CDS|UTR/ || $ensGene_impact ~/CDS|UTR/
		|| $acembly_impact ~/CDS|UTR/ || $genscan_impact ~/CDS|UTR/
	} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query1.tsv.lz4
} {}

test query {query 2} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -q {
		count($*_impact, ~/CDS|UTR/) > 0
	} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query1.tsv.lz4
} {}

test query {query 3} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -q {$type=="ins" && def($coverage-testNA19240chr2122cg,0)>=20} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query3.tsv.lz4
} {}

test query {query 4} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -q { $type=="snp" && $impact == "MISSENSE" } \
		-f {chromosome begin end} \
		annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query4.tsv.lz4
} {}

test query {query 5} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -h	annottestcompar.tsv.rz > temp.tsv
	cg tsvdiff temp.tsv answers/query5.tsv.lz4
} {}

test query {query 6 region} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -q {
		$chromosome == "22" && $begin < 16434940 && $end > 16054739
	} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query6.tsv.lz4
} {}

test query {query 6 region 2} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -q {
		region("22:16054739-16434940")
	} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query6.tsv.lz4
} {}

test query {query 6 region 3} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -q {
		region(22:16054739:16434940)
	} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query6.tsv.lz4
} {}

test query {query 6 region 4} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -q {
		region("chr22",16054739,16434940)
	} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query6.tsv.lz4
} {}

test query {query regions} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -q {
		region("chr22:16054739-16434940",21,9415561,9417128)
	} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query6.tsv.lz4
} {2,4d1
< 21	9415561	9415562	del	G		?	-	-	?	?	-496	143			u	20901814		G	55	55	-663	167			v																																1	p11.2				1													MER41D											-	-	-	-	-	-	-	-	-	-	-	-	-	-
< 21	9417126	9417127	snp	C	A	18820400	A	C	96	96	-319	101			v	?	C	C	?	?	-192	172			r																																1	p11.2				1													HAL1											-	-	-	-	-	-	-	-	-	-	-	-	-	-
< 21	9417127	9417128	snp	G	A	18820400	A	G	96	96	-410	101			v	20901846	A	G	137	137	-298	170			v	dbsnp.92:rs1850103																															1	p11.2				1													HAL1											-	-	-	-	-	-	-	-	-	-	-	-	rs1850103	0.5
child process exited abnormally} error

test query {query 7} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -q {
		sm(testNA19239chr2122cg, testNA19240chr2122cg)
	} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query7.tsv.lz4
} {}

test query {query 8} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -q {
		same(testNA19239chr2122cg, testNA19240chr2122cg)
	} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query8.tsv.lz4
} {}

test query {query 9} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -q {
		mm(testNA19239chr2122cg, testNA19240chr2122cg)
	} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query9.tsv.lz4
} {}

test query {query 10} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -q {
		$type == "snp" 
		&& $simpleRepeat == "" && $microsat == "" && $genomicSuperDups == "" 
	} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query10.tsv.lz4
} {}

test query {query 11} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -q {
		$sequenced-testNA19239chr2122cg == "v"
		&& $sequenced-testNA19240chr2122cg == "v"
	} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query11.tsv.lz4
} {}

test query {query 12} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -q {
		count($sequenced-*, == "v")  == 2
	} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query12.tsv.lz4
} {}

test query {query 13} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -q {
		$sequenced-testNA19239chr2122cg == "v"
		&& $type == "snp"
		&& $coverage-testNA19239chr2122cg >= 20 
		&& $coverage-testNA19239chr2122cg <= 100
		&& $cluster-testNA19239chr2122cg == ""
	} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query13.tsv.lz4
} {}

test query {query 14} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -q {
		count($sequenced-*, == "v")  == 4
		&& $type == "snp"
		&& count($coverage-*, <20) == 0
		&& count($coverage-*, > 100) == 0
		&& count($cluster-*, != "") == 0
	} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query14.tsv.lz4
} {}

test query {query 15} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -q {
		count($*_impact, ~/CDS|UTR/) > 0
	} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query15.tsv.lz4
} {}

test query {query 16} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -q {
		$ensGene_impact ~ /CDS/
	} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query16.tsv.lz4
} {}

test query {query 17} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -q {
		counthasone($ensGene_impact, == "CDSMIS") > 0
	} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query17.tsv.lz4
} {}

test query {query 18} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -q {
		counthasone($*_impact, ~ /CDS/) > 0
	} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query18.tsv.lz4
} {}

test query {query 19} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -q {
		counthasall($*_impact, ~ /CDS/) > 0
	} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query19.tsv.lz4
} {}

test query {query 20} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -f "chromosome begin end refGene_impact" -q {
		oneof($refGene_impact, "CDSMIS", "CDSINS") > 0
	} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query20.tsv.lz4
} {}

test query {query 20 2} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -f "chromosome begin end refGene_impact" -q {
		$refGene_impact in {"CDSMIS" "CDSINS"}
	} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query20.tsv.lz4
} {}

test query {query 21 combine} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -f "chromosome begin end 1000gCEU 1000glow" -q {
		count(lmax($1000gCEU),lmax($1000glow), > 0.98) > 0
	} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query21.tsv.lz4
} {}

test query {query 22} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -f "chromosome begin end refGene_impact knownGene_impact" -q {
		hasone($refGene_impact, "==", "CDSMIS") || hasone($knownGene_impact, "==", "CDSMIS")
	} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query22.tsv.lz4
} {}

test query {query 22 b} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -f "chromosome begin end refGene_impact knownGene_impact" -q {
		contains($refGene_impact, "CDSMIS") || ($knownGene_impact contains "CDSMIS")
	} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query22.tsv.lz4
} {}

test query {query 22 c} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -f "chromosome begin end refGene_impact knownGene_impact" -q {
		hasone($refGene_impact, == "CDSMIS") || hasone($knownGene_impact, == "CDSMIS")
	} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query22.tsv.lz4
} {}

test query {query 22 d} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -f "chromosome begin end refGene_impact knownGene_impact" -q {
		lone($refGene_impact @== "CDSMIS") || lone($knownGene_impact @== "CDSMIS")
	} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query22.tsv.lz4
} {}

test query {query 23 b} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv
	exec cg select -f "chromosome begin end refGene_impact knownGene_impact" -q {
		shares($refGene_impact, "CDSMIS CDSINS") && ($knownGene_impact shares "CDSMIS CDSINS")
	} annottestcompar.tsv.rz temp.tsv
	cg tsvdiff temp.tsv answers/query23.tsv.lz4
} {}

test query {liftover} {
	cd $::smalltestdir/testqueries ; file delete temp.tsv temp2.tsv
	exec cg liftover annottestcompar.tsv.rz temp.tsv /complgen/refseq/liftover/hg18ToHg19.over.chain
	cg select -q {$ROW < 100} temp.tsv temp2.tsv
	cg tsvdiff temp2.tsv answers/liftover.tsv.lz4
} {}

cd $keepdir

testsummarize
