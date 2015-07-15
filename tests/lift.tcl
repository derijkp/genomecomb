#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test liftover {basic no correctvariants} {
	test_cleantmp
	cg select -rf {test} data/expected-var_lift-hg18tohg19.tsv tmp/expected.tsv
	exec cg liftover -correctvariants 0 data/var_lift.tsv tmp/temp.tsv /complgen/refseq/liftover/hg18ToHg19.over.tsv
	cg select -rf {test} tmp/temp.tsv tmp/temp2.tsv
	exec diff tmp/temp.tsv.unmapped data/expected-var_lift-hg18tohg19.tsv.unmapped
	exec diff tmp/temp2.tsv tmp/expected.tsv
} {12,13c12,13
< 1	2492275	2492276	snp	C	A	B2	r	o	19	G	G	v	t	18	A	C	1	2482141	2482142	C
< 1	2492275	2492276	snp	C	G	B1	v	m	18	G	G	r	o	18	A	C	1	2482141	2482142	C
---
> 1	2492275	2492276	snp	G	C	B1	v	m	18	C	C	r	o	18	T	G	1	2482141	2482142	C
> 1	2492275	2492276	snp	G	T	B2	r	o	19	C	C	v	t	18	T	G	1	2482141	2482142	C
} error regexp

test liftover {basic liftover with correctvariants} {
	test_cleantmp
	exec cg liftover data/var_lift.tsv tmp/temp.tsv /complgen/refseq/liftover/hg18ToHg19.over.tsv
	exec diff tmp/temp.tsv data/expected-var_lift-hg18tohg19.tsv
	exec diff tmp/temp.tsv.unmapped data/expected-var_lift-hg18tohg19.tsv.unmapped
} {}

test liftregion {basic} {
	test_cleantmp
	exec cg liftregion data/reg_lift.tsv tmp/temp.tsv /complgen/refseq/liftover/hg18ToHg19.over.tsv
	exec diff tmp/temp.tsv data/expected-reg_lift-hg18tohg19.tsv
	exec diff tmp/temp.tsv.unmapped data/expected-reg_lift-hg18tohg19.tsv.unmapped
} {}

test liftregion {changed refseq} {
	test_cleantmp
	# 1: alt becomes ref
	# 2,3: two alleles (one becomes ref)
	# 4: complement with refchange to complement!
	# 5: complement with refchange to other base
	write_tab tmp/temp.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	score-sample1	alleleSeq1-sample1	alleleSeq2-sample1	sequenced-sample2	zyg-sample2	score-sample2	alleleSeq1-sample2	alleleSeq2-sample2
		1	1609629	1609630	snp	T	C	v	m	1	C	C	r	r	2	T	T
		1	2474628	2474629	snp	G	C	v	t	2	G	C	r	o	4	T	T
		1	2474628	2474629	snp	G	T	r	o	3	G	C	v	m	6	T	T
		1	2489102	2489103	snp	A	T	v	m	4	T	T	r	r	5	A	A
		2	89594172	89594173	snp	T	A	v	t	5	T	A	r	r	8	T	T
		2	110743803	110743804	snp	A	T	v	t	6	A	T	r	r	10	A	A
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	score-sample1	alleleSeq1-sample1	alleleSeq2-sample1	sequenced-sample2	zyg-sample2	score-sample2	alleleSeq1-sample2	alleleSeq2-sample2	hg18_chromosome	hg18_begin	hg18_end	hg18_ref
		1	1619766	1619767	snp	C	T	r	r	1	C	C	v	m	2	T	T	1	1609629	1609630	T
		1	2484768	2484769	snp	T	C	v	c	2	G	C	r	r	4	T	T	1	2474628	2474629	G
		1	2484768	2484769	snp	T	G	v	c	2	G	C	r	r	4	T	T	1	2474628	2474629	G
		1	2499242	2499243	snp	C	A	r	o	4	T	T	v	m	5	A	A	1	2489102	2489103	A
		1	2499242	2499243	snp	C	T	v	m	4	T	T	r	o	5	A	A	1	2489102	2489103	A
		2	89563941	89563942	snp	T	A	v	t	5	A	T	v	m	8	A	A	2	89594172	89594173	T
		2	110577444	110577445	snp	C	A	v	c	6	T	A	r	o	10	T	T	2	110743803	110743804	A
		2	110577444	110577445	snp	C	T	v	c	6	T	A	v	m	10	T	T	2	110743803	110743804	A
	}
	file delete tmp/lifted.tsv
	exec cg liftover tmp/temp.tsv tmp/lifted.tsv /complgen/refseq/liftover/hg18ToHg19.over.tsv
	cg select -rc 1 tmp/lifted.tsv tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test liftregion {changed refseq, unsplit} {
	test_cleantmp
	# 1: alt becomes ref
	# 2,3: two alleles (one becomes ref)
	# 4: complement with refchange to complement!
	# 5: complement with refchange to other base
	write_tab tmp/temp.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	score-sample1	alleleSeq1-sample1	alleleSeq2-sample1	sequenced-sample2	zyg-sample2	score-sample2	alleleSeq1-sample2	alleleSeq2-sample2
		1	1609629	1609630	snp	T	C	v	m	1	C	C	r	r	2	T	T
		1	2474628	2474629	snp	G	C,T	v	t	2	G	C	v	m	4	T	T
		1	2489102	2489103	snp	A	T	v	m	4	T	T	r	r	5	A	A
		2	89594172	89594173	snp	T	A	v	t	5	T	A	r	r	8	T	T
		2	110743803	110743804	snp	A	T	v	t	6	A	T	r	r	10	A	A
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	score-sample1	alleleSeq1-sample1	alleleSeq2-sample1	sequenced-sample2	zyg-sample2	score-sample2	alleleSeq1-sample2	alleleSeq2-sample2	hg18_chromosome	hg18_begin	hg18_end hg18_ref
		1	1619766	1619767	snp	C	T	r	r	1	C	C	v	m	2	T	T	1	1609629	1609630	T
		1	2484768	2484769	snp	T	C,G	v	c	2	G	C	r	r	4	T	T	1	2474628	2474629	G
		1	2499242	2499243	snp	C	A,T	v	m	4	T	T	v	m	5	A	A	1	2489102	2489103	A
		2	89563941	89563942	snp	T	A	v	t	5	A	T	v	m	8	A	A	2	89594172	89594173	T
		2	110577444	110577445	snp	C	A,T	v	c	6	T	A	v	m	10	T	T	2	110743803	110743804	A
	}
	file delete tmp/lifted.tsv
	exec cg liftover -split 0 tmp/temp.tsv tmp/lifted.tsv /complgen/refseq/liftover/hg18ToHg19.over.tsv
	cg select -rc 1 tmp/lifted.tsv tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test liftover {variants with different ref ending up in same spot -s 1} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		chromosome	begin	end	type	ref	alt
		2       109820576       109820577       snp     G       T
		2	110858379	110858380	snp	A	A
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	hg18_chromosome	hg18_begin	hg18_end	hg18_ref
		2       110463287       110463288       snp     G       T	2	109820576       109820577	G
	}
	exec cg liftover tmp/temp.tsv tmp/lifted.tsv /complgen/refseq/liftover/hg18ToHg19.over.tsv
	exec cg correctvariants -f 1 -s 1 tmp/lifted.tsv tmp/result.tsv.temp /complgen/refseq/hg19
	cg select -rc 1 tmp/result.tsv.temp tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test liftover {variants with different alt ending up in same spot -s 0} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		chromosome	begin	end	type	ref	alt
		2       109820576       109820577       snp     G       C
		2	110858379	110858380	snp	A	T
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	hg18_chromosome	hg18_begin	hg18_end	hg18_ref
		2       110463287       110463288       snp     G       A,C,T	2	110858379	110858380	A
	}
	exec cg liftover tmp/temp.tsv tmp/lifted.tsv /complgen/refseq/liftover/hg18ToHg19.over.tsv
	file delete tmp/result.tsv.temp
	exec cg correctvariants -f 1 -s 0 tmp/lifted.tsv tmp/result.tsv.temp /complgen/refseq/hg19
	cg select -rc 1 tmp/result.tsv.temp tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test liftover {changed refseq, regionsfile} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		chromosome	begin	end	type	ref	alt
		1	1610084	1610085	snp	A	G
	}
	write_tab tmp/reg.tsv {
		chromosome	begin	end
		1	609620	2474628
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	hg18_chromosome	hg18_begin	hg18_end	hg18_ref
		1	1619766	1619767	snp	C	T	1	1609629	1609630	T
		1	1620224	1620225	snp	G	A	1	1610084	1610085	A
		1	1620884	1620885	snp	T	C	1	1610744	1610745	C
	}
	file delete tmp/lifted.tsv
	exec cg liftover -regionfile tmp/reg.tsv tmp/temp.tsv tmp/lifted.tsv /complgen/refseq/liftover/hg18ToHg19.over.tsv
	cg select -rc 1 tmp/lifted.tsv tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test liftover {changed refseq, regionsfile with multiple samples} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	alleleSeq1-sample1	alleleSeq2-sample1	sequenced-sample2	zyg-sample2	alleleSeq1-sample2	alleleSeq2-sample2
		1	1610084	1610085	snp	A	G	v	m	G	G	r	r	A	A
	}
	write_tab tmp/reg.tsv {
		chromosome	begin	end	sreg-sample1 sreg-sample2
		1	609620	2474628	1	0
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	sequenced-sample1	zyg-sample1	alleleSeq1-sample1	alleleSeq2-sample1	sequenced-sample2	zyg-sample2	alleleSeq1-sample2	alleleSeq2-sample2	hg18_chromosome	hg18_begin	hg18_end	hg18_ref
		1	1619766	1619767	snp	C	T	v	m	T	T	u	u	?	?	1	1609629	1609630	T
		1	1620224	1620225	snp	G	A	r	r	G	G	v	m	A	A	1	1610084	1610085	A
		1	1620884	1620885	snp	T	C	v	m	C	C	u	u	?	?	1	1610744	1610745	C
	}
	file delete tmp/lifted.tsv
	exec cg liftover -regionfile tmp/reg.tsv tmp/temp.tsv tmp/lifted.tsv /complgen/refseq/liftover/hg18ToHg19.over.tsv
	# remove comments to compare
	cg select -rc 1 tmp/lifted.tsv tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test liftover {no alt} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		chromosome	begin	end	type	reference	sequenced	zyg	alleleSeq1	alleleSeq2
		1	1610084	1610085	snp	A	v	m	G	G
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	reference	alt	sequenced	zyg	alleleSeq1	alleleSeq2	hg18_chromosome	hg18_begin	hg18_end	hg18_ref
		1	1620224	1620225	snp	G	A	r	r	G	G	1	1610084	1610085	A
	}
	file delete tmp/lifted.tsv
	exec cg liftover tmp/temp.tsv tmp/lifted.tsv /complgen/refseq/liftover/hg18ToHg19.over.tsv
	# remove comments to compare
	cg select -rc 1 tmp/lifted.tsv tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test liftover {no alt split} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		chromosome	begin	end	type	reference	sequenced	zyg	alleleSeq1	alleleSeq2
		1	1610084	1610085	snp	A	v	m	C	A
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	reference	alt	sequenced	zyg	alleleSeq1	alleleSeq2	hg18_chromosome	hg18_begin	hg18_end	hg18_ref
		1	1620224	1620225	snp	G	A	v	c	C	A	1	1610084	1610085	A
		1	1620224	1620225	snp	G	C	v	c	C	A	1	1610084	1610085	A
	}
	file delete tmp/lifted.tsv
	exec cg liftover -split 1 tmp/temp.tsv tmp/lifted.tsv /complgen/refseq/liftover/hg18ToHg19.over.tsv
	# remove comments to compare
	cg select -rc 1 tmp/lifted.tsv tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test liftover {no alt split, t} {
	test_cleantmp
	write_tab tmp/temp.tsv {
		chromosome	begin	end	type	reference	sequenced	zyg	alleleSeq1	alleleSeq2
		1	1610084	1610085	snp	A	v	m	G	A
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	reference	alt	sequenced	zyg	alleleSeq1	alleleSeq2	hg18_chromosome	hg18_begin	hg18_end	hg18_ref
		1	1620224	1620225	snp	G	A	v	t	G	A	1	1610084	1610085	A
	}
	file delete tmp/lifted.tsv
	exec cg liftover -split 1 tmp/temp.tsv tmp/lifted.tsv /complgen/refseq/liftover/hg18ToHg19.over.tsv
	# remove comments to compare
	cg select -rc 1 tmp/lifted.tsv tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test liftsample {basic} {
	test_cleantmp
	file mkdir tmp/sample
	file mkdir tmp/esample
	write_tab tmp/sample/fannotvar-sample.tsv {
		chromosome	begin	end	type	ref	alt
		1	1610084	1610085	snp	A	G
	}
	write_tab tmp/esample/fannotvar-liftedsample.tsv {
		#liftover_source	tmp/sample/fannotvar-sample.tsv
		#liftover	/complgen/refseq/liftover/hg18ToHg19.over.tsv
		#split	1
		#oldref	hg18
		#ref	hg19
		chromosome	begin	end	type	ref	alt	hg18_chromosome	hg18_begin	hg18_end	hg18_ref
		1	1619766	1619767	snp	C	T	1	1609629	1609630	T
		1	1620224	1620225	snp	G	A	1	1610084	1610085	A
		1	1620884	1620885	snp	T	C	1	1610744	1610745	C
	}
	mklink tmp/sample/fannotvar-sample.tsv tmp/sample/var-cg-cg-sample.tsv
	mklink tmp/esample/fannotvar-liftedsample.tsv tmp/esample/var-cg-cg-liftedsample.tsv
	mklink tmp/sample/doesnotexist.tsv tmp/sample/var-wronglink-sample.tsv
	write_tab tmp/sample/sreg-cg-cg-sample.tsv {
		chromosome	begin	end
		1	609620	2474628
	}
	write_tab tmp/esample/sreg-cg-cg-liftedsample.tsv {
		#liftover_source	tmp/sample/sreg-cg-cg-sample.tsv
		#liftover	/complgen/refseq/liftover/hg18ToHg19.over.tsv
		#oldref	hg18
		#ref	hg19
		chromosome	begin	end	hg18_chromosome	hg18_begin	hg18_end
		1	619757	1619814	1	609620	2474628
		1	1620086	1620859	1	609620	2474628
		1	1620860	1620903	1	609620	2474628
		1	1620904	2484768	1	609620	2474628
	}
	write_tab tmp/esample/sreg-cg-cg-liftedsample.tsv.unmapped {
		#liftover_source	tmp/sample/sreg-cg-cg-sample.tsv
		#liftover_unmapped	/complgen/refseq/liftover/hg18ToHg19.over.tsv
		chromosome	begin	end	old_begin	old_end
		1	1609677	1609759	609620	2474628
		1	1609849	1609946	609620	2474628
		1	1610719	1610720	609620	2474628
		1	1610763	1610764	609620	2474628
	}
	write_tab tmp/sample/reg_cluster-cg-cg-sample.tsv {
		chromosome	begin	end
		1	609620	609720
		1	609820	609930
	}
	write_tab tmp/esample/reg_cluster-cg-cg-liftedsample.tsv {
		#liftover_source	tmp/sample/reg_cluster-cg-cg-sample.tsv
		#liftover	/complgen/refseq/liftover/hg18ToHg19.over.tsv
		#oldref	hg18
		#ref	hg19
		chromosome	begin	end	hg18_chromosome	hg18_begin	hg18_end
		1	619757	619857	1	609620	609720	
		1	619957	620067	1	609820	609930
	}
	# test
	exec cg liftsample -silent 1 tmp/sample tmp/liftedsample /complgen/refseq/liftover/hg18ToHg19.over.tsv
	file delete tmp/liftedsample/fannotvar-liftedsample.tsv.unmapped
	file delete tmp/liftedsample/reg_cluster-cg-cg-liftedsample.tsv.unmapped
	file delete tmp/liftedsample/sampleinfo.tsv
	file delete -force tmp/liftedsample/log_jobs
	if {[file link tmp/liftedsample/var-wronglink-liftedsample.tsv] ne "doesnotexist.tsv"} {
		error "(dangling) links not correctly done"
	}
	file delete -force tmp/liftedsample/var-wronglink-liftedsample.tsv
	exec diff -r tmp/liftedsample tmp/esample
} {}

testsummarize
