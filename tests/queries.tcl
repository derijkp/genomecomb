#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl
set keepdir [pwd]

test query {query 1} {
	file delete tmp/temp.tsv.zst
	exec cg select -f ROW -q {
		$refGene_impact ~/CDS|UTR/ || $knownGene_impact ~/CDS|UTR/ || $ensGene_impact ~/CDS|UTR/
		|| $acembly_impact ~/CDS|UTR/ || $genscan_impact ~/CDS|UTR/
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query1.tsv.zst
} {}

test query {query 2} {
	file delete tmp/temp.tsv.zst
	exec cg select -f ROW -q {
		count($*_impact, ~/CDS|UTR/) > 0
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query1.tsv.zst
} {}

test query {query 3} {
	file delete tmp/temp.tsv.zst
	exec cg select -f ROW -q {
		$type=="ins" && def($coverage-testNA19240chr2122cg,0)>=20
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query3.tsv.zst
} {}

test query {query 4} {
	file delete tmp/temp.tsv.zst
	exec cg select -f ROW -q { $type=="snp" && $impact == "MISSENSE" } \
		-f {chromosome begin end} \
		data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query4.tsv.zst
} {}

test query {query 5} {
	file delete tmp/temp.tsv.zst
	exec cg select -h	data/testqueries/annottestcompar.tsv.zst > tmp/temp.tsv
	cg tsvdiff tmp/temp.tsv data/testqueries/answers/query5.tsv
} {}

test query {query 6 region} {
	file delete tmp/temp.tsv.zst
	exec cg select -q {
		$chromosome == "22" && $begin < 16434940 && $end > 16054739
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query6.tsv.zst
} {}

test query {query 6 region 2} {
	file delete tmp/temp.tsv.zst
	exec cg select -q {
		region("22:16054739-16434940")
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query6.tsv.zst
} {}

test query {query 6 region 3} {
	file delete tmp/temp.tsv.zst
	exec cg select -q {
		region(22:16054739:16434940)
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query6.tsv.zst
} {}

test query {query 6 region 4} {
	file delete tmp/temp.tsv.zst
	exec cg select -q {
		region("chr22",16054739,16434940)
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query6.tsv.zst
} {}

test query {query regions} {
	file delete tmp/temp.tsv.zst
	exec cg select -q {
		region("chr22:16054739-16434940",21,9415561,9417128)
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query6.tsv.zst
} {diff tmp/temp.tsv.zst data/testqueries/answers/query6.tsv.zst
header
  chromosome	begin	end	type	ref	alt	locus-testNA19239chr2122cg	alleleSeq1-testNA19239chr2122cg	alleleSeq2-testNA19239chr2122cg	totalScore1-testNA19239chr2122cg	totalScore2-testNA19239chr2122cg	refscore-testNA19239chr2122cg	coverage-testNA19239chr2122cg	refcons-testNA19239chr2122cg	cluster-testNA19239chr2122cg	sequenced-testNA19239chr2122cg	locus-testNA19240chr2122cg	alleleSeq1-testNA19240chr2122cg	alleleSeq2-testNA19240chr2122cg	totalScore1-testNA19240chr2122cg	totalScore2-testNA19240chr2122cg	refscore-testNA19240chr2122cg	coverage-testNA19240chr2122cg	refcons-testNA19240chr2122cg	cluster-testNA19240chr2122cg	sequenced-testNA19240chr2122cg	xRef	geneId	mrnaAcc	proteinAcc	symbol	orientation	component	componentIndex	hasCodingRegion	impact	nucleotidePos	proteinPos	annotationRefSequence	sampleSequence	genomeRefSequence	pfam	acembly_impact	acembly_gene	acembly_descr	ensGene_impact	ensGene_gene	ensGene_descr	genscan_impact	genscan_gene	genscan_descr	knownGene_impact	knownGene_gene	knownGene_descr	refGene_impact	refGene_gene	refGene_descr	chainSelf	cytoBand	dgv	evofold	gad	genomicSuperDups	gwasCatalog_name	gwasCatalog_score	homopolymer_base	homopolymer_size	microsat	oreganno	phastConsElements46way_name	phastConsElements46way_score	phastConsElements46wayPlacental_name	phastConsElements46wayPlacental_score	phastConsElements46wayPrimates_name	phastConsElements46wayPrimates_score	rmsk	simpleRepeat	targetScanS_name	targetScanS_score	tfbsConsSites_name	tfbsConsSites_score	tRNAs	vistaEnhancers_name	vistaEnhancers_score	wgRna_name	wgRna_score	1000gCEU	1000gCHBxJPT	1000glow	1000gYRI	dmgcmt_freq	dmgcmt_total	dmgep_freq	dmgep_total	dmgnbd_freq	dmgnbd_total	dmgpg_freq	dmgpg_total	snp135_name	snp135_freq
2,4d1
< 21	9415561	9415562	del	G		?	-	-	?	?	-496	143			u	20901814		G	55	55	-663	167			v																																1	p11.2				1													MER41D											-	-	-	-	-	-	-	-	-	-	-	-	-	-
< 21	9417126	9417127	snp	C	A	18820400	A	C	96	96	-319	101			v	?	C	C	?	?	-192	172			r																																1	p11.2				1													HAL1											-	-	-	-	-	-	-	-	-	-	-	-	-	-
< 21	9417127	9417128	snp	G	A	18820400	A	G	96	96	-410	101			v	20901846	A	G	137	137	-298	170			v	dbsnp.92:rs1850103																															1	p11.2				1													HAL1											-	-	-	-	-	-	-	-	-	-	-	-	rs1850103	0.5
child process exited abnormally} error

test query {query 7} {
	file delete tmp/temp.tsv.zst
	exec cg select -f ROW -q {
		sm(testNA19239chr2122cg, testNA19240chr2122cg)
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query7.tsv.zst
} {}

test query {query 8} {
	file delete tmp/temp.tsv.zst
	exec cg select -f ROW -q {
		same(testNA19239chr2122cg, testNA19240chr2122cg)
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query8.tsv.zst
} {}

test query {query 9} {
	file delete tmp/temp.tsv.zst
	exec cg select -f ROW -q {
		mm(testNA19239chr2122cg, testNA19240chr2122cg)
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query9.tsv.zst
} {}

test query {query 10} {
	file delete tmp/temp.tsv.zst
	exec cg select -f ROW -q {
		$type == "snp" 
		&& $simpleRepeat == "" && $microsat == "" && $genomicSuperDups == "" 
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query10.tsv.zst
} {}

test query {query 11} {
	file delete tmp/temp.tsv.zst
	exec cg select -f ROW -q {
		$sequenced-testNA19239chr2122cg == "v"
		&& $sequenced-testNA19240chr2122cg == "v"
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query11.tsv.zst
} {}

test query {query 12} {
	file delete tmp/temp.tsv.zst
	exec cg select -f ROW -q {
		count($sequenced-*, == "v")  == 2
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query12.tsv.zst
} {}

test query {query 13} {
	file delete tmp/temp.tsv.zst
	exec cg select -f ROW -q {
		$sequenced-testNA19239chr2122cg == "v"
		&& $type == "snp"
		&& $coverage-testNA19239chr2122cg >= 20 
		&& $coverage-testNA19239chr2122cg <= 100
		&& $cluster-testNA19239chr2122cg == ""
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query13.tsv.zst
} {}

test query {query 14} {
	file delete tmp/temp.tsv.zst
	exec cg select -f ROW -q {
		count($sequenced-*, == "v")  == 4
		&& $type == "snp"
		&& count($coverage-*, <20) == 0
		&& count($coverage-*, > 100) == 0
		&& count($cluster-*, != "") == 0
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query14.tsv.zst
} {}

test query {query 15} {
	file delete tmp/temp.tsv.zst
	exec cg select -f ROW -q {
		count($*_impact, ~/CDS|UTR/) > 0
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query15.tsv.zst
} {}

test query {query 16} {
	file delete tmp/temp.tsv.zst
	exec cg select -f ROW -q {
		$ensGene_impact ~ /CDS/
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query16.tsv.zst
} {}

test query {query 17} {
	file delete tmp/temp.tsv.zst
	exec cg select -f ROW -q {
		counthasone($ensGene_impact, == "CDSMIS") > 0
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query17.tsv.zst
} {}

test query {query 18} {
	file delete tmp/temp.tsv.zst
	exec cg select -f ROW -q {
		counthasone($*_impact, ~ /CDS/) > 0
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query18.tsv.zst
} {}

test query {query 19} {
	file delete tmp/temp.tsv.zst
	exec cg select -f ROW -q {
		counthasall($*_impact, ~ /CDS/) > 0
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query19.tsv.zst
} {}

test query {query 20} {
	file delete tmp/temp.tsv.zst
	exec cg select -f "chromosome begin end refGene_impact" -q {
		oneof($refGene_impact, "CDSMIS", "CDSINS") > 0
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query20.tsv.zst
} {}

test query {query 20 2} {
	file delete tmp/temp.tsv.zst
	exec cg select -f "chromosome begin end refGene_impact" -q {
		$refGene_impact in {"CDSMIS" "CDSINS"}
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query20.tsv.zst
} {}

test query {query 21 combine} {
	file delete tmp/temp.tsv.zst
	exec cg select -f "chromosome begin end 1000gCEU 1000glow" -q {
		count(lmax($1000gCEU),lmax($1000glow), > 0.98) > 0
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query21.tsv.zst
} {}

test query {query 22} {
	file delete tmp/temp.tsv.zst
	exec cg select -f "chromosome begin end refGene_impact knownGene_impact" -q {
		hasone($refGene_impact, "==", "CDSMIS") || hasone($knownGene_impact, "==", "CDSMIS")
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query22.tsv.zst
} {}

test query {query 22 b} {
	file delete tmp/temp.tsv.zst
	exec cg select -f "chromosome begin end refGene_impact knownGene_impact" -q {
		contains($refGene_impact, "CDSMIS") || ($knownGene_impact contains "CDSMIS")
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query22.tsv.zst
} {}

test query {query 22 c} {
	file delete tmp/temp.tsv.zst
	exec cg select -f "chromosome begin end refGene_impact knownGene_impact" -q {
		hasone($refGene_impact, == "CDSMIS") || hasone($knownGene_impact, == "CDSMIS")
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query22.tsv.zst
} {}

test query {query 22 d} {
	file delete tmp/temp.tsv.zst
	exec cg select -f "chromosome begin end refGene_impact knownGene_impact" -q {
		lone($refGene_impact @== "CDSMIS") || lone($knownGene_impact @== "CDSMIS")
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query22.tsv.zst
} {}

test query {query 23 b} {
	file delete tmp/temp.tsv.zst
	exec cg select -f "chromosome begin end refGene_impact knownGene_impact" -q {
		shares($refGene_impact, "CDSMIS CDSINS") && ($knownGene_impact shares "CDSMIS CDSINS")
	} data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst
	cg tsvdiff tmp/temp.tsv.zst data/testqueries/answers/query23.tsv.zst
} {}

test query {liftover} {
	file delete tmp/temp.tsv.zst tmp/temp2.tsv.zst
	exec cg liftover data/testqueries/annottestcompar.tsv.zst tmp/temp.tsv.zst $::refseqdir/liftover/hg18ToHg19.over.tsv
	cg select -overwrite 1 -q {$ROW < 100} tmp/temp.tsv.zst tmp/temp2.tsv.zst
	cg tsvdiff tmp/temp2.tsv.zst data/testqueries/answers/liftover.tsv.zst
} {}

cd $keepdir

testsummarize
