#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test vcf2tsv {vcf2tsv} {
	exec cg vcf2tsv data/test1000glow.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test1000glow.vcf2tsv
} {}

test vcf2tsv {vcf2tsv} {
	file copy -force data/vars1.vcf tmp/vars1.vcf
	exec cg vcf2tsv -split 1 tmp/vars1.vcf tmp/temp.tsv
	exec cg checktsv tmp/temp.tsv
	cg splitalleles data/expected-vars1-var_annot.sft tmp/expected.tsv
	cg tsvdiff tmp/temp.tsv tmp/expected.tsv
} {diff tmp/temp.tsv tmp/expected.tsv
header diff
<extrafields: name quality filter zyg-sample1 phased-sample1 genotypes-sample1 zyg-sample2 phased-sample2 genotypes-sample2 totalcoverage allelecount totalallelecount
---
>extrafields: sequenced-sample1 sequenced-sample2 annot_name annot_freq*} match error

test vcf2tsv {vcf2tsv -split 0} {
	exec cg vcf2tsv -split 0 data/test.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test.vcf2tsv
} {}

test vcf2tsv {vcf2tsv -split 0 ins and del} {
	exec cg vcf2tsv -split 0 data/test2.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test2.vcf2tsv
} {}

test vcf2tsv {vcf2tsv -split 0 1000glow} {
	exec cg vcf2tsv -split 0 data/test1000glow.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test1000glow.vcf2tsv
} {3c3
< #split	0
---
> #split	1
child process exited abnormally} error

test vcf2tsv {vcf2tsv -split ori} {
	exec cg vcf2tsv -split ori data/test.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-testori.vcf2tsv
} {}

test vcf2tsv {vcf2tsv ins and del -split ori} {
	exec cg vcf2tsv -split ori data/test2.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test2ori.vcf2tsv
} {}

test vcf2tsv {vcf2tsv} {
	exec cg vcf2tsv -s 1 data/test.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-tests.vcf2tsv
} {}

test vcf2tsv {vcf2tsv combination alleles -> position insdel} {
	write_vcf tmp/test.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      SAMPLE
		chrtest	10	.	CTG	C	10	.	.	GT	0/1
		chrtest	20	.	CTGA	CT	20	.	.	GT	0/1
		chrtest	30	.	A	ATG,G	30	.	.	GT	0/1
		chrtest	40	.	GC	GCC,G	40	.	.	GT	1/2
		chrtest	50	.	AT	ATG,AG,CT	50	.	.	GT	1/2
		chrtest	60	.	ATTT	AT,ATT,ATTTT,ATTTTT,A	60	.	.	GT	1/2
		chrtest	70	.	GGT	G,GCCGT,GCCGTGGGGCTGGGGGCGTTGT	70	.	.	GT	1/2
		chrtest	80	.	TTTCTATTCTA	T,TTTCTA,TTTCTATTCTATTCTA	80	.	.	GT	1/2
		chrtest	90	.	ATG	ATGA,A,AT,CTG,AGG,ATC	90	.	.	GT	1/2
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	name	quality	filter	alleleSeq1	alleleSeq2	zyg	phased	genotypes	alleledepth_ref	alleledepth	TE	genoqual	coverage	haploqual	NS	totalcoverage	frequency	Ancestralallele	dbsnp	Hapmap2
		chrtest	10	12	del	TG	{}	.	10	.	TG	{}	t	0	0;1	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}
		chrtest	21	23	del	GA	{}	.	20	.	GA	{}	t	0	0;1	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}
		chrtest	29	30	snp	A	G	.	30	.	A	A	r	0	0;0	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}
		chrtest	30	30	ins	{}	TG	.	30	.	{}	TG	t	0	0;1	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}
		chrtest	40	40	ins	{}	C	.	40	.	C	{}	t	0	1;0	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}
		chrtest	40	41	del	C	{}	.	40	.	C	{}	t	0	0;1	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}
		chrtest	49	50	snp	A	C	.	50	.	A	A	r	0	0;0	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}
		chrtest	50	51	snp	T	G	.	50	.	T	G	t	0	0;1	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}
		chrtest	51	51	ins	{}	G	.	50	.	G	{}	t	0	1;0	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}
		chrtest	60	60	ins	{}	T	.	60	.	{}	{}	r	0	0;0	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}
		chrtest	60	60	ins	{}	TT	.	60	.	{}	{}	r	0	0;0	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}
		chrtest	60	61	del	T	{}	.	60	.	@	{}	c	0	2;1	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}
		chrtest	60	62	del	TT	{}	.	60	.	{}	@	c	0	1;2	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}
		chrtest	60	63	del	TTT	{}	.	60	.	@	@	o	0	2;2	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}
		chrtest	70	70	ins	{}	CC	.	70	.	{}	CC	t	0	0;1	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}
		chrtest	70	70	ins	{}	CCGTGGGGCTGGGGGCGTT	.	70	.	{}	CC	o	0	0;2	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}
		chrtest	70	72	del	GT	{}	.	70	.	{}	GT	t	0	1;0	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}
		chrtest	80	80	ins	{}	TTCTA	.	80	.	{}	{}	r	0	0;0	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}
		chrtest	80	85	del	TTCTA	{}	.	80	.	@	{}	c	0	2;1	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}
		chrtest	80	90	del	TTCTATTCTA	{}	.	80	.	{}	@	c	0	1;2	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}
		chrtest	89	90	snp	A	C	.	90	.	A	A	r	0	0;0	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}
		chrtest	90	91	snp	T	G	.	90	.	T	@	o	0	0;2	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}
		chrtest	90	92	del	TG	{}	.	90	.	TG	{}	t	0	0;1	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}
		chrtest	91	92	del	G	{}	.	90	.	G	@	o	0	0;2	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}
		chrtest	91	92	snp	G	C	.	90	.	G	@	o	0	0;2	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}
		chrtest	92	92	ins	{}	A	.	90	.	A	{}	t	0	1;0	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}
	}
	cg vcf2tsv tmp/test.vcf | cg select -rc 1 > tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test vcf2tsv {vcf2tsv ins and del} {
	exec cg vcf2tsv -s 1 data/test2.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test2s.vcf2tsv
} {}

test vcf2tsv {vcf2tsv ins and del -typelist .} {
	exec cg vcf2tsv -typelist . -s 1 data/test2.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test2s.vcf2tsv
} {28c28
< #fields	frequency	.	Float	Allele Frequency	info
---
> #fields	frequency	A	Float	Allele Frequency	info
46,47c46,47
< 20	1110695	1110696	snp	A	G	rs6040355	67	PASS	G	T	c	1	1,2	21	6	23,27	T	G	c	1	2,1	2	0	18,2	T	T	o	0	2;2	35	4		2	10	0.333,0.667	T	1	
< 20	1110695	1110696	snp	A	T	rs6040355	67	PASS	G	T	c	1	2,1	21	6	23,27	T	G	c	1	1,2	2	0	18,2	T	T	m	0	1;1	35	4		2	10	0.333,0.667	T	1	
---
> 20	1110695	1110696	snp	A	G	rs6040355	67	PASS	G	T	c	1	1,2	21	6	23,27	T	G	c	1	2,1	2	0	18,2	T	T	o	0	2;2	35	4		2	10	0.333	T	1	
> 20	1110695	1110696	snp	A	T	rs6040355	67	PASS	G	T	c	1	2,1	21	6	23,27	T	G	c	1	1,2	2	0	18,2	T	T	m	0	1;1	35	4		2	10	0.667	T	1	
child process exited abnormally} error

test vcf2tsv {vcf2tsv vars_mirna.vcf} {
	exec cg vcf2tsv -s 1 data/vars_mirna.vcf tmp/temp.tsv
	cg select -rc 1 -rf {name quality filter totalcoverage	allelecount	totalallelecount} tmp/temp.tsv tmp/temp2.tsv
	cg select -rc 1 -rf {name} data/vars_mirna.tsv tmp/expected.tsv
	exec diff tmp/temp2.tsv tmp/expected.tsv
} {22,23c22,23
< chr1	1102505	1102508	del	NNN	
< chr1	1102520	1102537	del	NNNNNNNNNNNNNNNNN	
---
> chr1	1102505	1102508	del	3	
> chr1	1102520	1102537	del	17	
child process exited abnormally} error

test vcf2tsv {vcf2tsv ins and del -split 1} {
	exec cg vcf2tsv -s 1 data/test2.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test2s.vcf2tsv
} {}

test vcf2tsv {cg_vcf2tsv make uppercase -split ori} {
	write_vcf tmp/temp.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
		20	14370	rs6054257	g	a	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	name	quality	filter	alleleSeq1-NA00001	alleleSeq2-NA00001	zyg-NA00001	phased-NA00001	genotypes-NA00001	alleledepth_ref-NA00001	alleledepth-NA00001	TE-NA00001	genoqual-NA00001	coverage-NA00001	haploqual-NA00001	alleleSeq1-NA00002	alleleSeq2-NA00002	zyg-NA00002	phased-NA00002	genotypes-NA00002	alleledepth_ref-NA00002	alleledepth-NA00002	TE-NA00002	genoqual-NA00002	coverage-NA00002	haploqual-NA00002	alleleSeq1-NA00003	alleleSeq2-NA00003	zyg-NA00003	phased-NA00003	genotypes-NA00003	alleledepth_ref-NA00003	alleledepth-NA00003	TE-NA00003	genoqual-NA00003	coverage-NA00003	haploqual-NA00003	NS	totalcoverage	frequency	Ancestralallele	dbsnp	Hapmap2
		20	14369	14370	snp	G	A	rs6054257	29	PASS	G	G	r	1	0,0	{}	{}	{}	48	1	51,51	A	G	t	1	1,0	{}	{}	{}	48	8	51,51	A	A	m	0	1;1	{}	{}	{}	43	5	.,.	3	14	0.5	{}	1	1
	}
	cg vcf2tsv -split ori tmp/temp.vcf tmp/temp.tsv
	cg select -rc 1 tmp/temp.tsv tmp/temp2.tsv
	exec diff tmp/temp2.tsv tmp/expected.tsv
#	cg tsvdiff tmp/temp.tsv tmp/expected.tsv
} {}

test vcf2tsv {cg_vcf2tsv make uppercase -split 1} {
	write_vcf tmp/temp.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
		20	14370	rs6054257	g	a	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	name	quality	filter	alleleSeq1-NA00001	alleleSeq2-NA00001	zyg-NA00001	phased-NA00001	genotypes-NA00001	alleledepth_ref-NA00001	alleledepth-NA00001	TE-NA00001	genoqual-NA00001	coverage-NA00001	haploqual-NA00001	alleleSeq1-NA00002	alleleSeq2-NA00002	zyg-NA00002	phased-NA00002	genotypes-NA00002	alleledepth_ref-NA00002	alleledepth-NA00002	TE-NA00002	genoqual-NA00002	coverage-NA00002	haploqual-NA00002	alleleSeq1-NA00003	alleleSeq2-NA00003	zyg-NA00003	phased-NA00003	genotypes-NA00003	alleledepth_ref-NA00003	alleledepth-NA00003	TE-NA00003	genoqual-NA00003	coverage-NA00003	haploqual-NA00003	NS	totalcoverage	frequency	Ancestralallele	dbsnp	Hapmap2
		20	14369	14370	snp	G	A	rs6054257	29	PASS	G	G	r	1	0,0	{}	{}	{}	48	1	51,51	A	G	t	1	1,0	{}	{}	{}	48	8	51,51	A	A	m	0	1;1	{}	{}	{}	43	5	.,.	3	14	0.5	{}	1	1
	}
	cg vcf2tsv -split 1 tmp/temp.vcf tmp/temp.tsv
	cg select -rc 1 tmp/temp.tsv tmp/temp2.tsv
	exec diff tmp/temp2.tsv tmp/expected.tsv
} {}

test vcf2tsv {cg_vcf2tsv make uppercase -split 1} {
	write_vcf tmp/temp.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
		20	14370	rs6054257	ttaa	t	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	name	quality	filter	alleleSeq1-NA00001	alleleSeq2-NA00001	zyg-NA00001	phased-NA00001	genotypes-NA00001	alleledepth_ref-NA00001	alleledepth-NA00001	TE-NA00001	genoqual-NA00001	coverage-NA00001	haploqual-NA00001	alleleSeq1-NA00002	alleleSeq2-NA00002	zyg-NA00002	phased-NA00002	genotypes-NA00002	alleledepth_ref-NA00002	alleledepth-NA00002	TE-NA00002	genoqual-NA00002	coverage-NA00002	haploqual-NA00002	alleleSeq1-NA00003	alleleSeq2-NA00003	zyg-NA00003	phased-NA00003	genotypes-NA00003	alleledepth_ref-NA00003	alleledepth-NA00003	TE-NA00003	genoqual-NA00003	coverage-NA00003	haploqual-NA00003	NS	totalcoverage	frequency	Ancestralallele	dbsnp	Hapmap2
		20	14370	14373	del	TAA	{}	rs6054257	29	PASS	TAA	TAA	r	1	0,0	{}	{}	{}	48	1	51,51	{}	TAA	t	1	1,0	{}	{}	{}	48	8	51,51	{}	{}	m	0	1;1	{}	{}	{}	43	5	.,.	3	14	0.5	{}	1	1
	}
	cg vcf2tsv -split 1 tmp/temp.vcf tmp/temp.tsv
	cg select -rc 1 tmp/temp.tsv tmp/temp2.tsv
	exec diff tmp/temp2.tsv tmp/expected.tsv
} {}

test vcf2tsv {cg_vcf2tsv -removefields} {
	write_vcf tmp/temp.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
		20	14370	rs6054257	ttaa	t	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	quality	alleleSeq1-NA00001	alleleSeq2-NA00001	zyg-NA00001	phased-NA00001	genotypes-NA00001	alleledepth_ref-NA00001	alleledepth-NA00001	TE-NA00001	genoqual-NA00001	coverage-NA00001	alleleSeq1-NA00002	alleleSeq2-NA00002	zyg-NA00002	phased-NA00002	genotypes-NA00002	alleledepth_ref-NA00002	alleledepth-NA00002	TE-NA00002	genoqual-NA00002	coverage-NA00002	alleleSeq1-NA00003	alleleSeq2-NA00003	zyg-NA00003	phased-NA00003	genotypes-NA00003	alleledepth_ref-NA00003	alleledepth-NA00003	TE-NA00003	genoqual-NA00003	coverage-NA00003	NS	totalcoverage	frequency	Ancestralallele
		20	14370	14373	del	TAA	{}	29	TAA	TAA	r	1	0,0	{}	{}	{}	48	1	{}	TAA	t	1	1,0	{}	{}	{}	48	8	{}	{}	m	0	1;1	{}	{}	{}	43	5	3	14	0.5	{}
	}
	cg vcf2tsv -split 1 -removefields {name filter DB H2 HQ} tmp/temp.vcf tmp/temp.tsv
	cg select -rc 1 tmp/temp.tsv tmp/temp2.tsv
	exec diff tmp/temp2.tsv tmp/expected.tsv
} {}

test vcf2tsv {cg_vcf2tsv -removefields only info} {
	write_vcf tmp/temp.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
		20	14370	rs6054257	ttaa	t	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	name	quality	filter	alleleSeq1-NA00001	alleleSeq2-NA00001	zyg-NA00001	phased-NA00001	genotypes-NA00001	alleledepth_ref-NA00001	alleledepth-NA00001	TE-NA00001	genoqual-NA00001	coverage-NA00001	haploqual-NA00001	alleleSeq1-NA00002	alleleSeq2-NA00002	zyg-NA00002	phased-NA00002	genotypes-NA00002	alleledepth_ref-NA00002	alleledepth-NA00002	TE-NA00002	genoqual-NA00002	coverage-NA00002	haploqual-NA00002	alleleSeq1-NA00003	alleleSeq2-NA00003	zyg-NA00003	phased-NA00003	genotypes-NA00003	alleledepth_ref-NA00003	alleledepth-NA00003	TE-NA00003	genoqual-NA00003	coverage-NA00003	haploqual-NA00003	totalcoverage
		20	14370	14373	del	TAA	{}	rs6054257	29	PASS	TAA	TAA	r	1	0,0	{}	{}	{}	48	1	51,51	{}	TAA	t	1	1,0	{}	{}	{}	48	8	51,51	{}	{}	m	0	1;1	{}	{}	{}	43	5	.,.	14
	}
	cg vcf2tsv -split 1 -removefields {AA NS AF DB H2} tmp/temp.vcf tmp/temp.tsv
	cg select -rc 1 tmp/temp.tsv tmp/temp2.tsv
	exec diff tmp/temp2.tsv tmp/expected.tsv
} {}

test vcf2tsv {cg_vcf2tsv -removefields only info, ref} {
	write_vcf tmp/temp.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001
		20	14370	.	A	.	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:DP	0/0:48
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	name	quality	filter	alleleSeq1	alleleSeq2	zyg	phased	genotypes	alleledepth_ref	alleledepth	TE	genoqual	coverage	haploqual totalcoverage
		20	14369	14370	snp	A	.	.	29	PASS	A	A	r	0	0;0	{}	{}	{}	{}	48	{}	14
	}
	cg vcf2tsv -split 1 -removefields {AA NS AF DB H2} tmp/temp.vcf tmp/temp.tsv
	cg select -rc 1 tmp/temp.tsv tmp/temp2.tsv
	exec diff tmp/temp2.tsv tmp/expected.tsv
} {}

test vcf2tsv {cg_vcf2tsv -keepfields} {
	write_vcf tmp/temp.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
		20	14370	rs6054257	ttaa	t	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	quality	alleleSeq1-NA00001	alleleSeq2-NA00001	zyg-NA00001	phased-NA00001	genotypes-NA00001	alleledepth_ref-NA00001	alleledepth-NA00001	genoqual-NA00001	alleleSeq1-NA00002	alleleSeq2-NA00002	zyg-NA00002	phased-NA00002	genotypes-NA00002	alleledepth_ref-NA00002	alleledepth-NA00002	genoqual-NA00002	alleleSeq1-NA00003	alleleSeq2-NA00003	zyg-NA00003	phased-NA00003	genotypes-NA00003	alleledepth_ref-NA00003	alleledepth-NA00003	genoqual-NA00003	totalcoverage
		20	14370	14373	del	TAA	{}	29	TAA	TAA	r	1	0,0	{}	{}	48	{}	TAA	t	1	1,0	{}	{}	48	{}	{}	m	0	1;1	{}	{}	43	14
	}
	cg vcf2tsv -split 1 -keepfields {genoqual alleledepth} tmp/temp.vcf tmp/temp.tsv
	cg select -rc 1 tmp/temp.tsv tmp/temp2.tsv
	exec diff tmp/temp2.tsv tmp/expected.tsv
} {}

test vcf2tsv {cg_vcf2tsv -keepfields including name filter and info field} {
	write_vcf tmp/temp.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
		20	14370	rs6054257	ttaa	t	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	name quality	filter	alleleSeq1-NA00001	alleleSeq2-NA00001	zyg-NA00001	phased-NA00001	genotypes-NA00001	alleledepth_ref-NA00001	alleledepth-NA00001	genoqual-NA00001	alleleSeq1-NA00002	alleleSeq2-NA00002	zyg-NA00002	phased-NA00002	genotypes-NA00002	alleledepth_ref-NA00002	alleledepth-NA00002	genoqual-NA00002	alleleSeq1-NA00003	alleleSeq2-NA00003	zyg-NA00003	phased-NA00003	genotypes-NA00003	alleledepth_ref-NA00003	alleledepth-NA00003	genoqual-NA00003	totalcoverage	Hapmap2
		20	14370	14373	del	TAA	{}	rs6054257	29	PASS	TAA	TAA	r	1	0,0	{}	{}	48	{}	TAA	t	1	1,0	{}	{}	48	{}	{}	m	0	1;1	{}	{}	43	14	1
	}
	cg vcf2tsv -split 1 -keepfields {name filter genoqual alleledepth totalcoverage Hapmap2} tmp/temp.vcf tmp/temp.tsv
	cg select -rc 1 tmp/temp.tsv tmp/temp2.tsv
	exec diff tmp/temp2.tsv tmp/expected.tsv
} {}

test vcf2tsv {vcf2tsv extraalleles} {
	exec cg vcf2tsv -split 1 data/extraalleles.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test_extraalleles.vcf2tsv
} {}

test vcf2tsv {vcf2tsv -split 0 space in name} {
	file copy -force data/test.vcf "tmp/test it.vcf"
	exec cg vcf2tsv -split 0 "tmp/test it.vcf" "tmp/test it.tsv"
	exec diff "tmp/test it.tsv" data/expected-test.vcf2tsv
} {}

test vcf2tsv {missing geno fields} {
	write_vcf tmp/temp.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
		20	14370	rs6054257	g	a	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0	.
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	name	quality	filter	alleleSeq1-NA00001	alleleSeq2-NA00001	zyg-NA00001	phased-NA00001	genotypes-NA00001	alleledepth_ref-NA00001	alleledepth-NA00001	TE-NA00001	genoqual-NA00001	coverage-NA00001	haploqual-NA00001	alleleSeq1-NA00002	alleleSeq2-NA00002	zyg-NA00002	phased-NA00002	genotypes-NA00002	alleledepth_ref-NA00002	alleledepth-NA00002	TE-NA00002	genoqual-NA00002	coverage-NA00002	haploqual-NA00002	alleleSeq1-NA00003	alleleSeq2-NA00003	zyg-NA00003	phased-NA00003	genotypes-NA00003	alleledepth_ref-NA00003	alleledepth-NA00003	TE-NA00003	genoqual-NA00003	coverage-NA00003	haploqual-NA00003	NS	totalcoverage	frequency	Ancestralallele	dbsnp	Hapmap2
		20	14369	14370	snp	G	A	rs6054257	29	PASS	G	G	r	1	0,0	{}	{}	{}	48	1	51,51	A	G	t	1	1,0	{}	{}	{}	{}	{}	{}	?	?	?	?	?	{}	{}	{}	{}	{}	{}	3	14	0.5	{}	1	1
	}
	cg vcf2tsv tmp/temp.vcf tmp/temp.tsv
	cg select -rc 1 tmp/temp.tsv tmp/temp2.tsv
	exec diff tmp/temp2.tsv tmp/expected.tsv
#	cg tsvdiff tmp/temp.tsv tmp/expected.tsv
} {}

test vcf2tsv {duplicate fields (in info and format)} {
	write_vcf tmp/temp.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
		20	14370	rs6054257	g	a	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ:DB	0|0:48:1:51,51:t	1|0	.
	} "##FORMAT=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">\n"
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	name	quality	filter	alleleSeq1-NA00001	alleleSeq2-NA00001	zyg-NA00001	phased-NA00001	genotypes-NA00001	alleledepth_ref-NA00001	alleledepth-NA00001	TE-NA00001	genoqual-NA00001	coverage-NA00001	haploqual-NA00001	dbsnp-NA00001	alleleSeq1-NA00002	alleleSeq2-NA00002	zyg-NA00002	phased-NA00002	genotypes-NA00002	alleledepth_ref-NA00002	alleledepth-NA00002	TE-NA00002	genoqual-NA00002	coverage-NA00002	haploqual-NA00002	dbsnp-NA00002	alleleSeq1-NA00003	alleleSeq2-NA00003	zyg-NA00003	phased-NA00003	genotypes-NA00003	alleledepth_ref-NA00003	alleledepth-NA00003	TE-NA00003	genoqual-NA00003	coverage-NA00003	haploqual-NA00003	dbsnp-NA00003	NS	totalcoverage	frequency	Ancestralallele	info_dbsnp	Hapmap2
		20	14369	14370	snp	G	A	rs6054257	29	PASS	G	G	r	1	0,0	{}	{}	{}	48	1	51,51	t	A	G	t	1	1,0	{}	{}	{}	{}	{}	{}	{}	?	?	?	?	?	{}	{}	{}	{}	{}	{}	{}	3	14	0.5	{}	1	1
	}
	cg vcf2tsv tmp/temp.vcf tmp/temp.tsv
	cg select -rc 1 tmp/temp.tsv tmp/temp2.tsv
	exec diff tmp/temp2.tsv tmp/expected.tsv
} {}

test vcf2tsv {wrong nr lines with -removefields} {
	cg vcf2tsv -removefields {
		name filter AN AC AF AA ExcessHet InbreedingCoeff MLEAC MLEAF NDA RPA RU STR
	} data/varall-freebayes-rdsbwa-NA19238mx2-part.vcf tmp/temp.tsv
	cg checktsv tmp/temp.tsv
} {}

# structural variants sv

test vcf2tsv {svtest.vcf} {
	cg vcf2tsv -split 1 data/svtest.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/svtest.tsv
} {}

test vcf2tsv {vcf2tsv -split 0 svtest.vcf} {
	cg vcf2tsv -split 0 data/svtest.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/svtest-unsplit.tsv
} {}

test vcf2tsv {vcf2tsv -split ori svtest.vcf} {
	cg vcf2tsv -split ori data/svtest.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/svtest-ori.tsv
} {}

test vcf2tsv {AD 0} {
	write_vcf tmp/temp.vcf {
		CHROM    POS    ID      REF     ALT     QUAL    FILTER  INFO    FORMAT   NA00001
		20	14370	.	G	.	35	PASS	{}	GT:GQ:AD	0|0:48:0
		20	14371	.	A	C	36	PASS	{}	GT:GQ:AD	1|1:48:0,18
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	name	quality	filter	alleleSeq1	alleleSeq2	zyg	phased	genotypes	alleledepth_ref	alleledepth	TE genoqual	coverage	haploqual NS	totalcoverage	frequency	Ancestralallele	dbsnp	Hapmap2
		20	14369	14370	snp	G	.	.	35	PASS	G          G          r   1      0,0       {}              {}          {} 48       0        {}        {} {}            {}        {}              {}         {}
		20	14370	14371	snp	A	C	.	36	PASS	C          C          m   1      1,1       0               18          {} 48       18        {}        {} {}            {}        {}              {}         {}
	}
	cg vcf2tsv tmp/temp.vcf tmp/result.tsv
	cg select -overwrite 1 -rc 1 tmp/result.tsv tmp/cresult.tsv
	exec diff tmp/cresult.tsv tmp/expected.tsv
} {}

test vcf2tsv {blocked gvcf} {
	write_vcf tmp/test.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      SAMPLE
		chr21	42775286	.	T	<NON_REF>	.	.	END=42775359	GT:AD:DP:GQ:PL	0/0:7,0:7:18:0,18,270
		chr21	42775360	.	TGTTTAA	T,<NON_REF>	0	.	RAW_MQ=70729.00	GT:AD:GQ:PL:SB	0/0:29,0,0:87:0,87,1242,87,1242,1242:18,11,0,0
		chr21	42775361	.	G	<NON_REF>	.	.	END=42775362	GT:AD:DP:GQ:PL	0/0:5,0:5:12:0,12,180		
	} {} {
		##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
	}
	file_write tmp/expected.tsv [deindent {
		chromosome	begin	end	type	ref	alt	name	quality	filter	alleleSeq1	alleleSeq2	zyg	phased	genotypes	alleledepth_ref	alleledepth	TE	genoqual	coverage	haploqual	NS	totalcoverage	frequency	Ancestralallele	dbsnp	Hapmap2
		chr21	42775285	42775359	ref	74	.	.	.	.	74	74	r	0	0;0	7	0		18	7							
		chr21	42775359	42775360	snp	T	.	.	.	.	T	T	r	0	0;0	29	0		87	29							
		chr21	42775360	42775362	ref	2	.	.	.	.	2	2	r	0	0;0	5	0		12	5							
		chr21	42775360	42775366	del	GTTTAA		.	0	.	GTTTAA	GTTTAA	r	0	0;0	29	0		87	29							
	}]\n
	cg vcf2tsv -refout 1 tmp/test.vcf tmp/result.tsv
	cg select -overwrite 1 -rc 1 tmp/result.tsv tmp/cresult.tsv
	exec diff tmp/cresult.tsv tmp/expected.tsv
} {}

test vcf2tsv {blocked gvcf -split ori} {
	write_vcf tmp/test.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      SAMPLE
		chr21	42775286	.	T	<NON_REF>	.	.	END=42775359	GT:AD:DP:GQ:PL	0/0:7,0:7:18:0,18,270
		chr21	42775360	.	TGTTTAA	T,<NON_REF>	0	.	RAW_MQ=70729.00	GT:AD:GQ:PL:SB	0/0:29,0,0:87:0,87,1242,87,1242,1242:18,11,0,0
		chr21	42775361	.	G	<NON_REF>	.	.	END=42775362	GT:AD:DP:GQ:PL	0/0:5,0:5:12:0,12,180		
		chr21	42775363	.	TGT	T,<NON_REF>	0	.	RAW_MQ=70730.00	GT:AD:DP:GQ:PL:SB	0/1:29,5,0:34:87:0,87,1242,87,1242,1242:18,11,0,0
	} {} {
		##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
	}
	file_write tmp/expected.tsv [deindent {
		chromosome	begin	end	type	ref	alt	name	quality	filter	alleleSeq1	alleleSeq2	zyg	phased	genotypes	alleledepth_ref	alleledepth	TE	genoqual	coverage	haploqual	NS	totalcoverage	frequency	Ancestralallele	dbsnp	Hapmap2
		chr21	42775285	42775359	ref	74	.	.	.	.	74	74	r	0	0;0	7	0		18	7							
		chr21	42775359	42775366	sub	TGTTTAA	T,.	.	0	.	TGTTTAA	TGTTTAA	r	0	0;0	29	0,0		87								
		chr21	42775360	42775362	ref	2	.	.	.	.	2	2	r	0	0;0	5	0		12	5							
		chr21	42775362	42775365	sub	TGT	T,.	.	0	.	TGT	T	t	0	0;1	29	5,0		87	34							
	}]\n
	cg vcf2tsv -refout 1 -split ori tmp/test.vcf tmp/result.tsv
	cg select -overwrite 1 -rc 1 tmp/result.tsv tmp/cresult.tsv
	exec diff tmp/cresult.tsv tmp/expected.tsv
} {}

test vcf2tsv {blocked gvcf 2} {
	file_write tmp/temp.vcf [bgvcfheader]\n[deindent {
		chr21	9452821	.	T	<NON_REF>	.	.	END=9452821	GT:DP:GQ:MIN_DP:PL	0/0:7:12:7:0,12,180
		chr21	9452822	.	T	<NON_REF>	.	.	END=9452826	GT:DP:GQ:MIN_DP:PL	0/0:7:9:7:0,9,135
		chr21	9452827	.	A	<NON_REF>	.	.	END=9452831	GT:DP:GQ:MIN_DP:PL	0/0:8:12:8:0,12,180
		chr21	9452832	.	TATC	T,<NON_REF>	76.73	.	AS_RAW_BaseQRankSum=|1.5,1|NaN;AS_RAW_MQ=11797.00|5840.00|0.00;AS_RAW_MQRankSum=|-1.4,1|NaN;AS_RAW_ReadPosRankSum=|1.0,1|NaN;AS_SB_TABLE=0,4|0,3|0,0;BaseQRankSum=1.579;DP=8;ExcessHet=3.0103;MLEAC=1,0;MLEAF=0.500,0.00;MQRankSum=-1.368;NDA=2;RAW_MQ=20038.00;ReadPosRankSum=1.006	GT:AD:GQ:PL:SB	0/1:4,3,0:99:114,0,159,126,168,294:0,4,0,3
		chr21	9452836	.	A	<NON_REF>	.	.	END=9452836	GT:DP:GQ:MIN_DP:PL	0/0:7:0:7:0,0,28
		chr21	9452837	.	T	<NON_REF>	.	.	END=9452838	GT:DP:GQ:MIN_DP:PL	0/0:7:21:7:0,21,296
	}]\n
	file_write tmp/expected.tsv [deindent {
		chromosome	begin	end	type	ref	alt	name	quality	filter	alleleSeq1	alleleSeq2	zyg	phased	genotypes	alleledepth_ref	alleledepth	coverage	genoqual	GQX	MIN_DP	PGT	PID	PL	SB	AS_InbreedingCoeff	AS_QD	AS_RAW_BaseQRankSum	AS_RAW_MQ	AS_RAW_MQRankSum	AS_RAW_ReadPosRankSum	AS_SB_TABLE	BaseQRankSum	totalcoverage	DS	ExcessHet	InbreedingCoeff	MLEAC	MLEAF	MQ	MQRankSum	NDA	RAW_MQ	ReadPosRankSum
		chr21	9452820	9452821	snp	T	.	.	.	.	T	T	r	0	0;0			7	12		7			0,12,180																				
		chr21	9452821	9452826	ref	5	.	.	.	.	5	5	r	0	0;0			7	9		7			0,9,135																				
		chr21	9452826	9452831	ref	5	.	.	.	.	5	5	r	0	0;0			8	12		8			0,12,180																				
		chr21	9452832	9452835	del	ATC		.	76.73	.	ATC		t	0	0;1	4	3	7	99					114,0,159,126,168,294	0,4,0,3			|1.5,1|NaN	11797.00|5840.00|0.00	|-1.4,1|NaN	|1.0,1|NaN	0,4|0,3|0,0	1.579	8		3.0103		1	0.500		-1.368	2	20038.00	1.006
		chr21	9452835	9452836	snp	A	.	.	.	.	A	A	r	0	0;0			7	0		7			0,0,28																				
		chr21	9452836	9452838	ref	2	.	.	.	.	2	2	r	0	0;0			7	21		7			0,21,296																				
	}]\n
	cg vcf2tsv tmp/temp.vcf tmp/temp.tsv
	cg select -overwrite 1 -rc 1 tmp/temp.tsv tmp/temp2.tsv
	exec diff tmp/temp2.tsv tmp/expected.tsv
} {}

test vcf2tsv {blocked gvcf no GQ in ref (strelka)} {
	file_write tmp/temp.vcf [bgvcfheader]\n[deindent {
		chr21	42775175	.	T	.	.	PASS	END=42775179;BLOCKAVG_min30p3a	GT:GQX:DP:DPF:MIN_DP	0/0:15:6:0:6
		chr21	42775180	.	C	T	112	LowGQX	SNVHPOL=8;MQ=60	GT:GQ:GQX:DP:DPF:AD:ADF:ADR:SB:FT:PL	1/1:15:0:6:0:0,6:0,5:0,1:-9.2:LowGQX:149,18,0
		chr21	42775181	.	T	.	.	PASS	END=42775188;BLOCKAVG_min30p3a	GT:GQX:DP:DPF:MIN_DP	0/0:15:6:0:6
	}]\n
	file_write tmp/expected.tsv [deindent {
		chromosome	begin	end	type	ref	alt	name	quality	filter	alleleSeq1	alleleSeq2	zyg	phased	genotypes	alleledepth_ref	alleledepth	coverage	genoqual	GQX	MIN_DP	PGT	PID	PL	SB	AS_InbreedingCoeff	AS_QD	AS_RAW_BaseQRankSum	AS_RAW_MQ	AS_RAW_MQRankSum	AS_RAW_ReadPosRankSum	AS_SB_TABLE	BaseQRankSum	totalcoverage	DS	ExcessHet	InbreedingCoeff	MLEAC	MLEAF	MQ	MQRankSum	NDA	RAW_MQ	ReadPosRankSum
		chr21	42775174	42775179	ref	5	.	.	.	PASS	5	5	r	0	0;0			6	15	15	6																							
		chr21	42775179	42775180	snp	C	T	.	112	LowGQX	T	T	m	0	1;1	0	6	6	15	0				149,18,0	-9.2															60				
		chr21	42775180	42775188	ref	8	.	.	.	PASS	8	8	r	0	0;0			6	15	15	6																							
	}]\n
	cg vcf2tsv tmp/temp.vcf tmp/temp.tsv
	cg select -overwrite 1 -rc 1 tmp/temp.tsv tmp/temp2.tsv
	exec diff tmp/temp2.tsv tmp/expected.tsv
} {}

test vcf2tsv {gvcf GVCF} {
	cg vcf2tsv -split 0 data/gatkh.gvcf tmp/test.tsv
	cg tsvdiff tmp/test.tsv data/gatkh-expected.tsv
} {}

test vcf2tsv {gvcf BP_RESOLUTION} {
	cg vcf2tsv -split 0 data/gatkh_bp.gvcf tmp/test.tsv
	cg tsvdiff tmp/test.tsv data/gatkh_bp-expected.tsv
} {}

test vcf2tsv {converted from gvcf BP_RESOLUTION} {
	cg gatk_genotypevcfs -dbdir $::refseqdir/hg19 data/varall-gatkh-bwa-sample1.gvcf tmp/test.vcf
	cg vcf2tsv tmp/test.vcf tmp/test.tsv
	cg tsvdiff tmp/test.tsv data/varall-gatkh-bwa-sample1.tsv
} {}

test vcf2tsv {vcf2tsv -split 1 -sort 0 ins and del} {
	exec cg vcf2tsv -split 0 -sort 0 data/test2.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test2.vcf2tsv
} {}

test vcf2tsv {vcf2tsv -sort 0 gatkh} {
	file copy data/varall-gatkh-bwa-sample2.gvcf tmp
	cg vcf2tsv -sort 0 tmp/varall-gatkh-bwa-sample2.gvcf tmp/result.tsv
	cg checksort tmp/result.tsv
} {}

test vcf2tsv {vcf2tsv -sort 0 testcases} {
	write_vcf tmp/test.gvcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      SAMPLE
		chrtest	1	.	T	<NON_REF>	.	.	.	GT:AD:DP:GQ:PL	0/0:7,0:7:18:0,18,270
		chrtest	2	.	C	T,<NON_REF>	71.03	.	AS_RAW_BaseQRankSum=||;AS_RAW_MQ=0.00|8002.00|0.00;AS_RAW_MQRankSum=||;AS_RAW_ReadPosRankSum=||;AS_SB_TABLE=0,0|2,2|0,0;DP=4;ExcessHet=3.0103;MLEAC=2,0;MLEAF=1.00,0.00;NDA=2;RAW_MQ=8002.00	GT:AD:GQ:PL:SB	1/1:0,4,0:12:99,12,0,99,12,99:0,0,2,2
		chrtest	3	.	CTG	C,<NON_REF>	871.73	.	AS_RAW_BaseQRankSum=|0.1,1|NaN;AS_RAW_MQ=96745.00|73763.00|0.00;AS_RAW_MQRankSum=|-1.8,1|NaN;AS_RAW_ReadPosRankSum=|-1.1,1|NaN;AS_SB_TABLE=25,3|23,1|0,0;BaseQRankSum=0.159;DP=53;ExcessHet=3.0103;MLEAC=1,0;MLEAF=0.500,0.00;MQRankSum=-1.711;NDA=2;RAW_MQ=173757.00;ReadPosRankSum=-1.084	GT:AD:GQ:PL:SB	0/1:28,24,0:99:909,0,1086,993,1159,2152:25,3,23,1
		chrtest	4	.	A	G,<NON_REF>	71.03	.	AS_RAW_BaseQRankSum=||;AS_RAW_MQ=0.00|8002.00|0.00;AS_RAW_MQRankSum=||;AS_RAW_ReadPosRankSum=||;AS_SB_TABLE=0,0|2,2|0,0;DP=4;ExcessHet=3.0103;MLEAC=2,0;MLEAF=1.00,0.00;NDA=2;RAW_MQ=8002.00	GT:AD:GQ:PL:SB	1/1:0,4,0:12:99,12,0,99,12,99:0,0,2,2
		chrtest	5	.	TTTTTAA	T,<NON_REF>	0	.	RAW_MQ=70729.00	GT:AD:GQ:PL:SB	0/0:29,0,0:87:0,87,1242,87,1242,1242:18,11,0,0
		chrtest	6	.	C	T,CGCT,<NON_REF>	362.73	.	AS_RAW_BaseQRankSum=|-0.2,1|0.2,1|NaN;AS_RAW_MQ=97200.00|54000.00|39600.00|0.00;AS_RAW_MQRankSum=|0.0,1|0.0,1|NaN;AS_RAW_ReadPosRankSum=|-0.9,1|-0.3,1|NaN;AS_SB_TABLE=16,11|9,6|4,7|0,0;BaseQRankSum=0.036;DP=55;ExcessHet=3.0103;MLEAC=1,0,0;MLEAF=0.500,0.00,0.00;MQRankSum=0.000;NDA=3;RAW_MQ=198000.00;ReadPosRankSum=-0.677	GT:AD:GQ:PL:SB	0/1:27,15,11,0:46:400,0,918,46,481,1238,507,968,1154,1592:16,11,13,13
		chrtest	7	.	G	A,<NON_REF>	71.03	.	AS_RAW_BaseQRankSum=||;AS_RAW_MQ=0.00|8002.00|0.00;AS_RAW_MQRankSum=||;AS_RAW_ReadPosRankSum=||;AS_SB_TABLE=0,0|2,2|0,0;DP=4;ExcessHet=3.0103;MLEAC=2,0;MLEAF=1.00,0.00;NDA=2;RAW_MQ=8002.00	GT:AD:GQ:PL:SB	1/1:0,4,0:12:99,12,0,99,12,99:0,0,2,2
		chrtest	8	.	A	ATG,G,<NON_REF>	0.75	.	AS_RAW_BaseQRankSum=|-0.6,1|0.6,1|NaN;AS_RAW_MQ=49754.00|2329.00|1600.00|0.00;AS_RAW_MQRankSum=|-2.1,1|-1.2,1|NaN;AS_RAW_ReadPosRankSum=|-0.2,1|0.8,1|NaN;AS_SB_TABLE=17,1|0,2|0,1|0,0;BaseQRankSum=-0.102;DP=21;ExcessHet=3.0103;MLEAC=1,0,0;MLEAF=0.500,0.00,0.00;MQRankSum=-2.304;NDA=3;RAW_MQ=53683.00;ReadPosRankSum=0.252	GT:AD:GQ:PL:SB	0/1:18,2,1,0:30:30,0,750,42,714,795,84,756,798,840:17,1,0,3
		chrtest	9	.	CT	C,CTTTT,CTTTTT,<NON_REF>	244.73	.	AS_RAW_BaseQRankSum=|1.0,1|0.2,1|-0.1,1|NaN;AS_RAW_MQ=68049.00|21600.00|6736.00|8321.00|0.00;AS_RAW_MQRankSum=|0.6,1|-4.3,1|-3.4,1|NaN;AS_RAW_ReadPosRankSum=|-0.4,1|0.6,1|0.2,1|NaN;AS_SB_TABLE=4,15|3,3|0,4|0,4|0,0;BaseQRankSum=0.664;DP=55;ExcessHet=3.0103;MLEAC=0,1,0,0;MLEAF=0.00,0.500,0.00,0.00;MQRankSum=-3.036;NDA=4;RAW_MQ=181910.00;ReadPosRankSum=0.173	GT:AD:GQ:PL:SB	0/2:19,6,4,4,0:2:280,235,731,0,343,848,2,351,843,1047,355,719,889,943,1245:4,15,3,11
		chrtest	10	.	G	C,T,<NON_REF>	0	.	AS_RAW_BaseQRankSum=|-0.4,1|-1.6,1|NaN;AS_RAW_MQ=186624.00|8036.00|16068.00|0.00;AS_RAW_MQRankSum=|-1.3,1|-1.2,1|NaN;AS_RAW_ReadPosRankSum=|0.4,1|1.0,1|NaN;AS_SB_TABLE=59,95|3,5|3,12|0,0;BaseQRankSum=-1.454;DP=177;ExcessHet=3.0103;MLEAC=0,0,0;MLEAF=0.00,0.00,0.00;MQRankSum=-1.587;NDA=3;RAW_MQ=210728.00;ReadPosRankSum=1.021	GT:AD:GQ:PL:SB	0/0:154,8,15,0:26:0,225,5596,26,5095,5275,466,5583,5334,5810:59,95,6,17
		chrtest	11	.	GAAAAAAAA	G,GAA,GAAA,GAAAAAA,GAAAAAAA,GAAAAAAAAA,<NON_REF>	239.73	.	AS_RAW_BaseQRankSum=|2.2,1|2.2,1|-0.4,1|1.2,1|1.3,1|0.4,1|NaN;AS_RAW_MQ=16767.00|11114.00|4886.00|2316.00|13634.00|10378.00|7297.00|0.00;AS_RAW_MQRankSum=|-2.9,1|-2.5,1|-3.1,1|-1.1,1|-1.8,1|-0.7,1|NaN;AS_RAW_ReadPosRankSum=|0.9,1|0.2,1|1.8,1|0.5,1|0.4,1|0.6,1|NaN;AS_SB_TABLE=2,6|5,7|1,5|1,3|2,6|4,4|1,3|0,0;BaseQRankSum=1.884;DP=99;ExcessHet=3.0103;MLEAC=1,0,0,0,0,0,0;MLEAF=0.500,0.00,0.00,0.00,0.00,0.00,0.00;MQRankSum=-2.735;NDA=7;RAW_MQ=117810.00;ReadPosRankSum=1.009	GT:AD:GQ:PL:SB	0/1:8,12,6,4,8,8,4,0:99:277,0,1968,105,1521,1837,176,1445,1654,1786,132,460,565,636,742,136,145,250,321,356,420,280,233,338,408,258,214,528,433,1284,1388,1457,880,578,660,1605:2,6,14,28
		chrtest	12	.	CCTGCTGTGACAGTTCCCTGCATGCAGGGCAGGAGTGTGTGCTTCTTCCCAGCAAAGGCAGAGTCAGGGCCTACAGAAACTGTGCCCACAGCCTAT	C,TCTGCTGTGACAGTTCCCTGCATGCAGGGCAGGAGTGTGTGCTTCTTCCCAGCAAAGGCAGAGTCAGGGCCTACAGAAACTGTGCCCACAGCCTAT,CAGTGATGGGGTCTGCTGTGACAGTTCCCTGCATGCAGGGCAGGAGTGTGTGCTTCTTCCCAGCAAAGGCAGAGTCAGGGCCTACAGAAACTGTGCCCACAGCCTAT,<NON_REF>	4.59	.	AS_RAW_BaseQRankSum=|-1.2,1|-1.9,1|NaN|NaN;AS_RAW_MQ=39710.00|6329.00|5200.00|0.00|0.00;AS_RAW_MQRankSum=|-1.8,1|-0.5,1|NaN|NaN;AS_RAW_ReadPosRankSum=|-2.3,1|0.5,1|NaN|NaN;AS_SB_TABLE=4,9|3,1|1,1|0,0|0,0;BaseQRankSum=-1.884;DP=25;ExcessHet=3.0103;MLEAC=0,1,0,0;MLEAF=0.00,0.500,0.00,0.00;MQRankSum=-1.661;NDA=4;RAW_MQ=70064.00;ReadPosRankSum=-1.493	GT:AD:GQ:PL:SB	0/2:13,4,2,0,0:9:39,9,889,0,538,662,39,645,570,674,96,781,627,695,851:4,9,4,2
	}
	cg vcf2tsv tmp/test.gvcf tmp/expected.tsv
	cg vcf2tsv -sort 0 tmp/test.gvcf tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test vcf2tsv {vcf2tsv -sort 0 split vcf} {
	write_vcf tmp/test.gvcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      SAMPLE
		chrtest	1	.	T	<NON_REF>	.	.	.	GT:AD:DP:GQ:PL	0/0:7,0:7:18:0,18,270
		chrtest	2	.	C	T,<NON_REF>	71.03	.	AS_RAW_BaseQRankSum=||;AS_RAW_MQ=0.00|8002.00|0.00;AS_RAW_MQRankSum=||;AS_RAW_ReadPosRankSum=||;AS_SB_TABLE=0,0|2,2|0,0;DP=4;ExcessHet=3.0103;MLEAC=2,0;MLEAF=1.00,0.00;NDA=2;RAW_MQ=8002.00	GT:AD:GQ:PL:SB	1/1:0,4,0:12:99,12,0,99,12,99:0,0,2,2
		chrtest	2	.	C	CT,<NON_REF>	71.04	.	AS_RAW_BaseQRankSum=||;AS_RAW_MQ=0.00|8002.00|0.00;AS_RAW_MQRankSum=||;AS_RAW_ReadPosRankSum=||;AS_SB_TABLE=0,0|2,2|0,0;DP=4;ExcessHet=3.0103;MLEAC=2,0;MLEAF=1.00,0.00;NDA=2;RAW_MQ=8002.00	GT:AD:GQ:PL:SB	1/1:0,4,0:12:99,12,0,99,12,99:0,0,2,2
		chrtest	2	.	CTG	C,<NON_REF>	871.73	.	AS_RAW_BaseQRankSum=|0.1,1|NaN;AS_RAW_MQ=96745.00|73763.00|0.00;AS_RAW_MQRankSum=|-1.8,1|NaN;AS_RAW_ReadPosRankSum=|-1.1,1|NaN;AS_SB_TABLE=25,3|23,1|0,0;BaseQRankSum=0.159;DP=53;ExcessHet=3.0103;MLEAC=1,0;MLEAF=0.500,0.00;MQRankSum=-1.711;NDA=2;RAW_MQ=173757.00;ReadPosRankSum=-1.084	GT:AD:GQ:PL:SB	0/1:28,24,0:99:909,0,1086,993,1159,2152:25,3,23,1
		chrtest	2	.	A	G,<NON_REF>	71.03	.	AS_RAW_BaseQRankSum=||;AS_RAW_MQ=0.00|8002.00|0.00;AS_RAW_MQRankSum=||;AS_RAW_ReadPosRankSum=||;AS_SB_TABLE=0,0|2,2|0,0;DP=4;ExcessHet=3.0103;MLEAC=2,0;MLEAF=1.00,0.00;NDA=2;RAW_MQ=8002.00	GT:AD:GQ:PL:SB	1/1:0,4,0:12:99,12,0,99,12,99:0,0,2,2
		chrtest	3	.	T	A,<NON_REF>	80.01	.	AS_RAW_BaseQRankSum=||;AS_RAW_MQ=0.00|8002.00|0.00;AS_RAW_MQRankSum=||;AS_RAW_ReadPosRankSum=||;AS_SB_TABLE=0,0|2,2|0,0;DP=4;ExcessHet=3.0103;MLEAC=2,0;MLEAF=1.00,0.00;NDA=2;RAW_MQ=8002.00	GT:AD:GQ:PL:SB	1/1:0,4,0:12:99,12,0,99,12,99:0,0,2,2
		chrtest	5	.	TTTTTAA	T,<NON_REF>	0	.	RAW_MQ=70729.00	GT:AD:GQ:PL:SB	0/0:29,0,0:87:0,87,1242,87,1242,1242:18,11,0,0
	}
	cg vcf2tsv tmp/test.gvcf tmp/expected.tsv
	cg vcf2tsv -sort 0 tmp/test.gvcf tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test vcf2tsv {vcf2tsv -refout 1 gatkh} {
	write_vcf tmp/test.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      SAMPLE
		chr21	42775286	.	T	<NON_REF>	.	.	.	GT:AD:DP:GQ:PL	0/0:7,0:7:18:0,18,270
		chr21	42775287	.	C	CGAGCT,<NON_REF>	63.73	.	ReadPosRankSum=1.027	GT:AD:GQ:PL:SB	0/1:307,26,0:99:101,0,12932,1041,13014,14055:169,138,7,19
		chr21	42775288	.	A	<NON_REF>	.	.	.	GT:AD:DP:GQ:PL	0/0:7,0:7:18:0,18,270
		chr21	42775359	.	A	<NON_REF>	.	.	.	GT:AD:DP:GQ:PL	0/0:5,0:5:12:0,12,180
		chr21	42775360	.	TGTTTAA	T,<NON_REF>	0	.	RAW_MQ=70729.00	GT:AD:GQ:PL:SB	0/0:29,0,0:87:0,87,1242,87,1242,1242:18,11,0,0
		chr21	42775361	.	G	<NON_REF>	.	.	.	GT:AD:DP:GQ:PL	0/0:5,0:5:12:0,12,180		
		chr21	42775362	.	CT	C,CG,CTTTTT,<NON_REF>	244.73	.	AS_RAW_BaseQRankSum=|1.0,1|0.2,1|-0.1,1|NaN;AS_RAW_MQ=68049.00|21600.00|6736.00|8321.00|0.00;AS_RAW_MQRankSum=|0.6,1|-4.3,1|-3.4,1|NaN;AS_RAW_ReadPosRankSum=|-0.4,1|0.6,1|0.2,1|NaN;AS_SB_TABLE=4,15|3,3|0,4|0,4|0,0;BaseQRankSum=0.664;DP=55;ExcessHet=3.0103;MLEAC=0,1,0,0;MLEAF=0.00,0.500,0.00,0.00;MQRankSum=-3.036;NDA=4;RAW_MQ=181910.00;ReadPosRankSum=0.173	GT:AD:GQ:PL:SB	0/2:19,6,4,4,0:2:280,235,731,0,343,848,2,351,843,1047,355,719,889,943,1245:4,15,3,11
		chr21	42775364	.	C	<NON_REF>	.	.	.	GT:AD:DP:GQ:PL	0/0:5,0:5:12:0,12,180		
		chr21	42775365	.	TA	T,CA,TAC,<NON_REF>	500.00	.	AS_RAW_BaseQRankSum=|1.0,1|0.2,1|-0.1,1|NaN;AS_RAW_MQ=68049.00|21600.00|6736.00|8321.00|0.00;AS_RAW_MQRankSum=|0.6,1|-4.3,1|-3.4,1|NaN;AS_RAW_ReadPosRankSum=|-0.4,1|0.6,1|0.2,1|NaN;AS_SB_TABLE=4,15|3,3|0,4|0,4|0,0;BaseQRankSum=0.664;DP=55;ExcessHet=3.0103;MLEAC=0,1,0,0;MLEAF=0.00,0.500,0.00,0.00;MQRankSum=-3.036;NDA=4;RAW_MQ=181910.00;ReadPosRankSum=0.173	GT:AD:GQ:PL:SB	0/2:19,6,4,4,0:2:280,235,731,0,343,848,2,351,843,1047,355,719,889,943,1245:4,15,3,11
	}
	file_write tmp/expected.tsv [deindent {
		chromosome	begin	end	type	ref	alt	name	quality	filter	alleleSeq1	alleleSeq2	zyg	phased	genotypes	alleledepth_ref	alleledepth	TE	genoqual	coverage	haploqual	NS	totalcoverage	frequency	Ancestralallele	dbsnp	Hapmap2
		chr21	42775285	42775286	snp	T	.	.	.	.	T	T	r	0	0;0	7	0		18	7							
		chr21	42775286	42775287	snp	C	.	.	.	.	C	C	r	0	0;0	307	0		99	333							
		chr21	42775287	42775287	ins		GAGCT	.	63.73	.		GAGCT	t	0	0;1	307	26		99	333							
		chr21	42775287	42775288	snp	A	.	.	.	.	A	A	r	0	0;0	7	0		18	7							
		chr21	42775358	42775359	snp	A	.	.	.	.	A	A	r	0	0;0	5	0		12	5							
		chr21	42775359	42775360	snp	T	.	.	.	.	T	T	r	0	0;0	29	0		87	29							
		chr21	42775360	42775361	snp	G	.	.	.	.	G	G	r	0	0;0	5	0		12	5							
		chr21	42775360	42775366	del	GTTTAA		.	0	.	GTTTAA	GTTTAA	r	0	0;0	29	0		87	29							
		chr21	42775361	42775362	snp	C	.	.	.	.	C	C	r	0	0;0	19	0		2	33			55				
		chr21	42775362	42775362	ins		TTTT	.	244.73	.			r	0	0;0	19	4		2	33			55				
		chr21	42775362	42775363	del	T		.	244.73	.	T	@	o	0	0;2	19	6		2	33			55				
		chr21	42775362	42775363	snp	T	G	.	244.73	.	T	G	t	0	0;1	19	4		2	33			55				
		chr21	42775363	42775364	snp	C	.	.	.	.	C	C	r	0	0;0	5	0		12	5							
		chr21	42775364	42775365	snp	T	C	.	500.00	.	T	C	t	0	0;1	19	4		2	33			55				
		chr21	42775365	42775366	del	A		.	500.00	.	A	A	r	0	0;0	19	6		2	33			55				
		chr21	42775366	42775366	ins		C	.	500.00	.			r	0	0;0	19	4		2	33			55				
	}]\n
	cg vcf2tsv -refout 1 tmp/test.vcf tmp/result.tsv
	cg select -overwrite 1 -rc 1 tmp/result.tsv tmp/cresult.tsv
	exec diff tmp/cresult.tsv tmp/expected.tsv
} {}

test vcf2tsv {vcf2tsv -sort 0 all vcf files} {
	foreach file [glob data/*vcf] {
		cg vcf2tsv -split 1 -sort 0 $file tmp/test.tsv
		cg vcf2tsv -split 1 $file tmp/expected.tsv
		exec diff tmp/test.tsv tmp/expected.tsv
	}
} {}

test vcf2tsv {calc DP from AD if not present} {
	write_vcf tmp/test.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      SAMPLE
		chr21	42775287	.	C	CGAGCT,<NON_REF>	63.73	.	ReadPosRankSum=1.027	GT:AD:GQ:PL:SB	0/1:307,26,0:99:101,0,12932,1041,13014,14055:169,138,7,19
		chr21	42775360	.	TTTTTAA	T,<NON_REF>	0	.	AS_RAW_BaseQRankSum=|NaN|NaN;AS_RAW_MQ=62329.00|0.00|0.00;AS_RAW_MQRankSum=|NaN|NaN;AS_RAW_ReadPosRankSum=|NaN|NaN;AS_SB_TABLE=18,11|0,0|0,0;DP=33;ExcessHet=3.0103;MLEAC=0,0;MLEAF=0.00,0.00;NDA=2;RAW_MQ=70729.00	GT:AD:GQ:PL:SB	0/0:29,0,0:87:0,87,1242,87,1242,1242:18,11,0,0
	}
	file_write tmp/expected.tsv [deindent {
		chromosome	begin	end	type	ref	alt	name	quality	filter	alleleSeq1	alleleSeq2	zyg	phased	genotypes	alleledepth_ref	alleledepth	TE	genoqual	coverage	haploqual	NS	totalcoverage	frequency	Ancestralallele	dbsnp	Hapmap2
		chr21	42775286	42775287	snp	C	.	.	.	.	C	C	r	0	0;0	307	0		99	333							
		chr21	42775287	42775287	ins		GAGCT	.	63.73	.		GAGCT	t	0	0;1	307	26		99	333							
		chr21	42775359	42775360	snp	T	.	.	.	.	T	T	r	0	0;0	29	0		87	29			33				
		chr21	42775360	42775366	del	TTTTAA		.	0	.	TTTTAA	TTTTAA	r	0	0;0	29	0		87	29			33				
	}]\n
	cg vcf2tsv -refout 1 tmp/test.vcf tmp/result.tsv
	cg select -overwrite 1 -rc 1 tmp/result.tsv tmp/cresult.tsv
	exec diff tmp/cresult.tsv tmp/expected.tsv
} {}

test vcf2tsv {sv with end < begin} {
	file_write tmp/test.vcf [deindent {
		##fileformat=VCFv4.2
		##source=Sniffles
		##fileDate=20180915
		##contig=<ID=chr3,length=198022430>
		##ALT=<ID=DEL,Description="Deletion">
		##ALT=<ID=DUP,Description="Duplication">
		##ALT=<ID=INV,Description="Inversion">
		##ALT=<ID=INVDUP,Description="InvertedDUP with unknown boundaries">
		##ALT=<ID=TRA,Description="Translocation">
		##ALT=<ID=INS,Description="Insertion">
		##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">
		##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
		##INFO=<ID=MAPQ,Number=1,Type=Integer,Description="Median mapping quality of paired-ends">
		##INFO=<ID=RE,Number=1,Type=Integer,Description="read support">
		##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
		##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variation">
		##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">
		##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">
		##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
		##INFO=<ID=STD_quant_start,Number=A,Type=Integer,Description="STD of the start breakpoints across the reads.">
		##INFO=<ID=STD_quant_stop,Number=A,Type=Integer,Description="STD of the stop breakpoints across the reads.">
		##INFO=<ID=Kurtosis_quant_start,Number=A,Type=Integer,Description="Kurtosis value of the start breakpoints accross the reads.">
		##INFO=<ID=Kurtosis_quant_stop,Number=A,Type=Integer,Description="Kurtosis value of the stop breakpoints accross the reads.">
		##INFO=<ID=SUPTYPE,Number=1,Type=String,Description="Type by which the variant is supported.(SR,ALN)">
		##INFO=<ID=SUPTYPE,Number=1,Type=String,Description="Type by which the variant is supported.(SR,ALN)">
		##INFO=<ID=STRANDS,Number=A,Type=String,Description="Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)">
		##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency.">
		##INFO=<ID=ZMW,Number=A,Type=Integer,Description="Number of ZMWs (Pacbio) supporting SV.">
		##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
		##FORMAT=<ID=DR,Number=1,Type=Integer,Description="# high-quality reference reads">
		##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# high-quality variant reads">
		#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	test
		chr3	184156712	22375	N	<INS>	.	UNRESOLVED	IMPRECISE;SVMETHOD=Snifflesv1.0.8;CHR2=chr3;END=184156710;STD_quant_start=3.987480;STD_quant_stop=491.614280;Kurtosis_quant_start=-1.349856;Kurtosis_quant_stop=-1.173662;SVTYPE=INS;SUPTYPE=AL,SR;SVLEN=999999999;STRANDS=+-;RE=17;AF=1	GT:DR:DV	1/1:0:17
	}]\n
	set result {}
	lappend result [catch {
		cg vcf2tsv tmp/test.vcf tmp/result1.tsv
	} msg] $msg
	lappend result [catch {
		cg vcf2tsv -locerror test tmp/test.vcf tmp/result1.tsv
	} msg] $msg
	cg vcf2tsv -locerror keep tmp/test.vcf tmp/result2.tsv
	lappend result [cg select -rc 1 -f {begin end} tmp/result2.tsv]
	cg vcf2tsv -locerror correct tmp/test.vcf tmp/result3.tsv
	lappend result [cg select -rc 1 -f {begin end} tmp/result3.tsv]
	set result
} {1 {END position 184156710 < begin 184156712
error converting vcf file: child process exited abnormally} 0 {} {begin	end
184156712	184156710} {begin	end
184156712	184156712}}

test vcf2tsv {vcfheader2tsv vs vcf2tsv} {
	foreach file [glob data/*.vcf] {
		exec cg vcfheader2tsv $file tmp/tempheader.tsv
		exec cg vcf2tsv $file tmp/temp.tsv
		exec cg select -overwrite 1 -header tmp/tempheader.tsv tmp/header.header
		exec cg select -overwrite 1 -header tmp/temp.tsv tmp/temp.header
		exec diff tmp/header.header tmp/temp.header
	}
} {}

test vcf2tsv {output various varcallers} {
	foreach file [glob data/var-*-bwa-test.vcf] {
		puts $file
		cg vcf2tsv -split 1 $file tmp/test.tsv
		cg tsvdiff tmp/test.tsv [file root $file].tsv
		exec diff tmp/test.tsv [file root $file].tsv
	}
} {}

test vcf2tsv {vcf2tsv error not enough fields in line} {
	write_vcf tmp/test.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      SAMPLE
		chrtest	10	.	CTG	C	10	.	.	GT	0/1
		chrtest	20	.	CTGA	CT	20	.	.	GT
		chrtest	30	.	A	ATG,G	30	.	.	GT	0/1
	}
	cg vcf2tsv tmp/test.vcf | cg select -rc 1 > tmp/result.tsv
} {not enough fields in line:
chrtest	20	.	CTGA	CT	20	.	.	GT
*} error match

test vcf2tsv {genotype not specified} {
	write_vcf tmp/test.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      SAMPLE
		chrtest	10	.	CTG	C	10	.	.	GT	./.
	}
	file_write tmp/expected.tsv [deindent {
		chromosome	begin	end	type	ref	alt	name	quality	filter	alleleSeq1	alleleSeq2	zyg	phased	genotypes	alleledepth_ref	alleledepth	TE	genoqual	coverage	haploqual	NS	totalcoverage	frequency	Ancestralallele	dbsnp	Hapmap2
		chrtest	10	12	del	TG		.	10	.	?	?	v	0	?;?.												
	}]\n
	cg vcf2tsv tmp/test.vcf tmp/result.tsv
	cg select -overwrite 1 -rc 1 tmp/result.tsv tmp/cresult.tsv
	exec diff tmp/cresult.tsv tmp/expected.tsv
} {}

test vcf2tsv {exec vcf2tsv} {
	set tempfile [tempfile]
	# list vcf2tsv 1 {. AD R RPA R AC A AF A} - - {} 0 * error 0 $tempfile {} < data/test1000glow.vcf > tmp/temp.tsv
	exec vcf2tsv 1 {. AD R RPA R AC A AF A} - - {} 0 * error 0 $tempfile {} < data/test1000glow.vcf > tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test1000glow.vcf2tsv
} {}

testsummarize
