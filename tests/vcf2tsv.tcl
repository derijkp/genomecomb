#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

proc write_vcf {file data {extracomment {}}} {
	set data [split [string trim $data] \n]
	set f [open $file w]
	puts $f [deindent {
		##fileformat=VCFv4.0
		##fileDate=20090805
		##source=myImputationProgramV3.1
		##reference=1000GenomesPilot-NCBI36
		##phasing=partial
		##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
		##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
		##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
		##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
		##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
		##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
		##FILTER=<ID=q10,Description="Quality below 10">
		##FILTER=<ID=s50,Description="Less than 50% of samples have data">
		##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
		##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
		##FORMAT=<ID=TE,Number=A,Type=Integer,Description="test for alt alleles in the order listed">
		##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
		##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
		##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
	}]
	if {$extracomment ne ""} {puts -nonewline $f $extracomment}
	set header [lindex $data 0]
	set data [lrange $data 1 end]
	puts $f \#[join $header \t]
	foreach line $data {
		puts $f [join $line \t]
	}
	close $f
}

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
} {27,28c27,28
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

test vcf2tsv {svtest.vcf} {
	cg vcf2tsv -split 1 data/svtest.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/svtest.tsv
} {}

test vcf2tsv {vcf2tsv -split 0 svtest.vcf} {
	cg vcf2tsv -split 0 data/svtest.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/svtest-unsplit.tsv
} {}

test vcf2tsv {svtest.vcf} {
	cg vcf2tsv -split ori data/svtest.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/svtest.tsv
} {*structural variants not supported in -split ori, vcf has SVLEN*} match error

test vcf2tsv {AD 0} {
	write_vcf tmp/temp.vcf {
		CHROM    POS    ID      REF     ALT     QUAL    FILTER  INFO    FORMAT   NA00001
		20	14370	.	G	.	35	PASS	{}	GT:GQ:AD	0|0:48:0
		20	14371	.	A	C	36	PASS	{}	GT:GQ:AD	1|1:48:0,18
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	type	ref	alt	name	quality	filter	alleleSeq1	alleleSeq2	zyg	phased	genotypes	alleledepth_ref	alleledepth	TE genoqual	coverage	haploqual NS	totalcoverage	frequency	Ancestralallele	dbsnp	Hapmap2
		20	14369	14370	snp	G	.	.	35	PASS	G          G          r   1      0,0       {}              {}          {} 48       {}        {}        {} {}            {}        {}              {}         {}
		20	14370	14371	snp	A	C	.	36	PASS	C          C          m   1      1,1       0               18          {} 48       {}        {}        {} {}            {}        {}              {}         {}
	}
	cg vcf2tsv tmp/temp.vcf tmp/temp.tsv
	cg select -overwrite 1 -rc 1 tmp/temp.tsv tmp/temp2.tsv
	exec diff tmp/temp2.tsv tmp/expected.tsv
} {}

test vcf2tsv {blocked gvcf} {
	file_write tmp/temp.vcf [deindent {
		##fileformat=VCFv4.2
		##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
		##FILTER=<ID=LowQual,Description="Low quality">
		##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
		##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
		##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
		##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
		##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
		##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another">
		##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
		##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
		##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.">
		##GATKCommandLine=<ID=HaplotypeCaller,CommandLine="HaplotypeCaller  --annotation-group StandardAnnotation --annotation-group AS_StandardAnnotation --emit-ref-confidence GVCF --annotate-with-num-discovered-alleles true --output varall-gatkh-rdsbwa-NA19240chr2122.gvcf.gz.temp.gz --intervals /tmp/tempExtral.10600-UMRUbdTebwJFtLcIEVGB/_Extral_temp_1.bed --input /data/genomecomb.testdata/tmp/exomes_yri_chr2122/samples/NA19240chr2122/map-rdsbwa-NA19240chr2122.bam --reference /complgen/refseq/hg19/genome_hg19.fa  --disable-tool-default-annotations false --gvcf-gq-bands 1 --gvcf-gq-bands 2 --gvcf-gq-bands 3 --gvcf-gq-bands 4 --gvcf-gq-bands 5 --gvcf-gq-bands 6 --gvcf-gq-bands 7 --gvcf-gq-bands 8 --gvcf-gq-bands 9 --gvcf-gq-bands 10 --gvcf-gq-bands 11 --gvcf-gq-bands 12 --gvcf-gq-bands 13 --gvcf-gq-bands 14 --gvcf-gq-bands 15 --gvcf-gq-bands 16 --gvcf-gq-bands 17 --gvcf-gq-bands 18 --gvcf-gq-bands 19 --gvcf-gq-bands 20 --gvcf-gq-bands 21 --gvcf-gq-bands 22 --gvcf-gq-bands 23 --gvcf-gq-bands 24 --gvcf-gq-bands 25 --gvcf-gq-bands 26 --gvcf-gq-bands 27 --gvcf-gq-bands 28 --gvcf-gq-bands 29 --gvcf-gq-bands 30 --gvcf-gq-bands 31 --gvcf-gq-bands 32 --gvcf-gq-bands 33 --gvcf-gq-bands 34 --gvcf-gq-bands 35 --gvcf-gq-bands 36 --gvcf-gq-bands 37 --gvcf-gq-bands 38 --gvcf-gq-bands 39 --gvcf-gq-bands 40 --gvcf-gq-bands 41 --gvcf-gq-bands 42 --gvcf-gq-bands 43 --gvcf-gq-bands 44 --gvcf-gq-bands 45 --gvcf-gq-bands 46 --gvcf-gq-bands 47 --gvcf-gq-bands 48 --gvcf-gq-bands 49 --gvcf-gq-bands 50 --gvcf-gq-bands 51 --gvcf-gq-bands 52 --gvcf-gq-bands 53 --gvcf-gq-bands 54 --gvcf-gq-bands 55 --gvcf-gq-bands 56 --gvcf-gq-bands 57 --gvcf-gq-bands 58 --gvcf-gq-bands 59 --gvcf-gq-bands 60 --gvcf-gq-bands 70 --gvcf-gq-bands 80 --gvcf-gq-bands 90 --gvcf-gq-bands 99 --indel-size-to-eliminate-in-ref-model 10 --use-alleles-trigger false --disable-optimizations false --just-determine-active-regions false --dont-genotype false --dont-trim-active-regions false --max-disc-ar-extension 25 --max-gga-ar-extension 300 --padding-around-indels 150 --padding-around-snps 20 --kmer-size 10 --kmer-size 25 --dont-increase-kmer-sizes-for-cycles false --allow-non-unique-kmers-in-ref false --num-pruning-samples 1 --recover-dangling-heads false --do-not-recover-dangling-branches false --min-dangling-branch-length 4 --consensus false --max-num-haplotypes-in-population 128 --error-correct-kmers false --min-pruning 2 --debug-graph-transformations false --kmer-length-for-read-error-correction 25 --min-observations-for-kmer-to-be-solid 20 --likelihood-calculation-engine PairHMM --base-quality-score-threshold 18 --pair-hmm-gap-continuation-penalty 10 --pair-hmm-implementation FASTEST_AVAILABLE --pcr-indel-model CONSERVATIVE --phred-scaled-global-read-mismapping-rate 45 --native-pair-hmm-threads 4 --native-pair-hmm-use-double-precision false --debug false --use-filtered-reads-for-annotations false --bam-writer-type CALLED_HAPLOTYPES --dont-use-soft-clipped-bases false --capture-assembly-failure-bam false --error-correct-reads false --do-not-run-physical-phasing false --min-base-quality-score 10 --smith-waterman JAVA --use-new-qual-calculator false --heterozygosity 0.001 --indel-heterozygosity 1.25E-4 --heterozygosity-stdev 0.01 --standard-min-confidence-threshold-for-calling 10.0 --max-alternate-alleles 6 --max-genotype-count 1024 --sample-ploidy 2 --genotyping-mode DISCOVERY --contamination-fraction-to-filter 0.0 --output-mode EMIT_VARIANTS_ONLY --all-site-pls false --min-assembly-region-size 50 --max-assembly-region-size 300 --assembly-region-padding 100 --max-reads-per-alignment-start 50 --active-probability-threshold 0.002 --max-prob-propagation-distance 50 --interval-set-rule UNION --interval-padding 0 --interval-exclusion-padding 0 --interval-merging-rule ALL --read-validation-stringency SILENT --seconds-between-progress-updates 10.0 --disable-sequence-dictionary-validation false --create-output-bam-index true --create-output-bam-md5 false --create-output-variant-index true --create-output-variant-md5 false --lenient false --add-output-sam-program-record true --add-output-vcf-command-line true --cloud-prefetch-buffer 40 --cloud-index-prefetch-buffer -1 --disable-bam-index-caching false --help false --version false --showHidden false --verbosity INFO --QUIET false --use-jdk-deflater false --use-jdk-inflater false --gcs-max-retries 20 --disable-tool-default-read-filters false --minimum-mapping-quality 20",Version=4.0.3.0,Date="April 19, 2018 5:09:59 PM CEST">
		##GVCFBlock0-1=minGQ=0(inclusive),maxGQ=1(exclusive)
		##GVCFBlock1-2=minGQ=1(inclusive),maxGQ=2(exclusive)
		##GVCFBlock10-11=minGQ=10(inclusive),maxGQ=11(exclusive)
		##GVCFBlock11-12=minGQ=11(inclusive),maxGQ=12(exclusive)
		##GVCFBlock12-13=minGQ=12(inclusive),maxGQ=13(exclusive)
		##GVCFBlock13-14=minGQ=13(inclusive),maxGQ=14(exclusive)
		##GVCFBlock14-15=minGQ=14(inclusive),maxGQ=15(exclusive)
		##GVCFBlock15-16=minGQ=15(inclusive),maxGQ=16(exclusive)
		##GVCFBlock16-17=minGQ=16(inclusive),maxGQ=17(exclusive)
		##GVCFBlock17-18=minGQ=17(inclusive),maxGQ=18(exclusive)
		##GVCFBlock18-19=minGQ=18(inclusive),maxGQ=19(exclusive)
		##GVCFBlock19-20=minGQ=19(inclusive),maxGQ=20(exclusive)
		##GVCFBlock2-3=minGQ=2(inclusive),maxGQ=3(exclusive)
		##GVCFBlock20-21=minGQ=20(inclusive),maxGQ=21(exclusive)
		##GVCFBlock21-22=minGQ=21(inclusive),maxGQ=22(exclusive)
		##GVCFBlock22-23=minGQ=22(inclusive),maxGQ=23(exclusive)
		##GVCFBlock23-24=minGQ=23(inclusive),maxGQ=24(exclusive)
		##GVCFBlock24-25=minGQ=24(inclusive),maxGQ=25(exclusive)
		##GVCFBlock25-26=minGQ=25(inclusive),maxGQ=26(exclusive)
		##GVCFBlock26-27=minGQ=26(inclusive),maxGQ=27(exclusive)
		##GVCFBlock27-28=minGQ=27(inclusive),maxGQ=28(exclusive)
		##GVCFBlock28-29=minGQ=28(inclusive),maxGQ=29(exclusive)
		##GVCFBlock29-30=minGQ=29(inclusive),maxGQ=30(exclusive)
		##GVCFBlock3-4=minGQ=3(inclusive),maxGQ=4(exclusive)
		##GVCFBlock30-31=minGQ=30(inclusive),maxGQ=31(exclusive)
		##GVCFBlock31-32=minGQ=31(inclusive),maxGQ=32(exclusive)
		##GVCFBlock32-33=minGQ=32(inclusive),maxGQ=33(exclusive)
		##GVCFBlock33-34=minGQ=33(inclusive),maxGQ=34(exclusive)
		##GVCFBlock34-35=minGQ=34(inclusive),maxGQ=35(exclusive)
		##GVCFBlock35-36=minGQ=35(inclusive),maxGQ=36(exclusive)
		##GVCFBlock36-37=minGQ=36(inclusive),maxGQ=37(exclusive)
		##GVCFBlock37-38=minGQ=37(inclusive),maxGQ=38(exclusive)
		##GVCFBlock38-39=minGQ=38(inclusive),maxGQ=39(exclusive)
		##GVCFBlock39-40=minGQ=39(inclusive),maxGQ=40(exclusive)
		##GVCFBlock4-5=minGQ=4(inclusive),maxGQ=5(exclusive)
		##GVCFBlock40-41=minGQ=40(inclusive),maxGQ=41(exclusive)
		##GVCFBlock41-42=minGQ=41(inclusive),maxGQ=42(exclusive)
		##GVCFBlock42-43=minGQ=42(inclusive),maxGQ=43(exclusive)
		##GVCFBlock43-44=minGQ=43(inclusive),maxGQ=44(exclusive)
		##GVCFBlock44-45=minGQ=44(inclusive),maxGQ=45(exclusive)
		##GVCFBlock45-46=minGQ=45(inclusive),maxGQ=46(exclusive)
		##GVCFBlock46-47=minGQ=46(inclusive),maxGQ=47(exclusive)
		##GVCFBlock47-48=minGQ=47(inclusive),maxGQ=48(exclusive)
		##GVCFBlock48-49=minGQ=48(inclusive),maxGQ=49(exclusive)
		##GVCFBlock49-50=minGQ=49(inclusive),maxGQ=50(exclusive)
		##GVCFBlock5-6=minGQ=5(inclusive),maxGQ=6(exclusive)
		##GVCFBlock50-51=minGQ=50(inclusive),maxGQ=51(exclusive)
		##GVCFBlock51-52=minGQ=51(inclusive),maxGQ=52(exclusive)
		##GVCFBlock52-53=minGQ=52(inclusive),maxGQ=53(exclusive)
		##GVCFBlock53-54=minGQ=53(inclusive),maxGQ=54(exclusive)
		##GVCFBlock54-55=minGQ=54(inclusive),maxGQ=55(exclusive)
		##GVCFBlock55-56=minGQ=55(inclusive),maxGQ=56(exclusive)
		##GVCFBlock56-57=minGQ=56(inclusive),maxGQ=57(exclusive)
		##GVCFBlock57-58=minGQ=57(inclusive),maxGQ=58(exclusive)
		##GVCFBlock58-59=minGQ=58(inclusive),maxGQ=59(exclusive)
		##GVCFBlock59-60=minGQ=59(inclusive),maxGQ=60(exclusive)
		##GVCFBlock6-7=minGQ=6(inclusive),maxGQ=7(exclusive)
		##GVCFBlock60-70=minGQ=60(inclusive),maxGQ=70(exclusive)
		##GVCFBlock7-8=minGQ=7(inclusive),maxGQ=8(exclusive)
		##GVCFBlock70-80=minGQ=70(inclusive),maxGQ=80(exclusive)
		##GVCFBlock8-9=minGQ=8(inclusive),maxGQ=9(exclusive)
		##GVCFBlock80-90=minGQ=80(inclusive),maxGQ=90(exclusive)
		##GVCFBlock9-10=minGQ=9(inclusive),maxGQ=10(exclusive)
		##GVCFBlock90-99=minGQ=90(inclusive),maxGQ=99(exclusive)
		##GVCFBlock99-100=minGQ=99(inclusive),maxGQ=100(exclusive)
		##INFO=<ID=AS_InbreedingCoeff,Number=A,Type=Float,Description="allele specific heterozygosity as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation; relate to inbreeding coefficient">
		##INFO=<ID=AS_QD,Number=A,Type=Float,Description="Allele-specific Variant Confidence/Quality by Depth">
		##INFO=<ID=AS_RAW_BaseQRankSum,Number=1,Type=String,Description="raw data for allele specific rank sum test of base qualities">
		##INFO=<ID=AS_RAW_MQ,Number=1,Type=String,Description="Allele-specfic raw data for RMS Mapping Quality">
		##INFO=<ID=AS_RAW_MQRankSum,Number=1,Type=String,Description="Allele-specfic raw data for Mapping Quality Rank Sum">
		##INFO=<ID=AS_RAW_ReadPosRankSum,Number=1,Type=String,Description="allele specific raw data for rank sum test of read position bias">
		##INFO=<ID=AS_SB_TABLE,Number=1,Type=String,Description="Allele-specific forward/reverse read counts for strand bias tests">
		##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
		##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
		##INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
		##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
		##INFO=<ID=ExcessHet,Number=1,Type=Float,Description="Phred-scaled p-value for exact test of excess heterozygosity">
		##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
		##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
		##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
		##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
		##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
		##INFO=<ID=NDA,Number=1,Type=Integer,Description="Number of alternate alleles discovered (but not necessarily genotyped) at this site">
		##INFO=<ID=RAW_MQ,Number=1,Type=Float,Description="Raw data for RMS Mapping Quality">
		##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
		##contig=<ID=chr1,length=249250621>
		##contig=<ID=chr2,length=243199373>
		##contig=<ID=chr3,length=198022430>
		##contig=<ID=chr4,length=191154276>
		##contig=<ID=chr5,length=180915260>
		##contig=<ID=chr6,length=171115067>
		##contig=<ID=chr7,length=159138663>
		##contig=<ID=chr8,length=146364022>
		##contig=<ID=chr9,length=141213431>
		##contig=<ID=chr10,length=135534747>
		##contig=<ID=chr11,length=135006516>
		##contig=<ID=chr12,length=133851895>
		##contig=<ID=chr13,length=115169878>
		##contig=<ID=chr14,length=107349540>
		##contig=<ID=chr15,length=102531392>
		##contig=<ID=chr16,length=90354753>
		##contig=<ID=chr17,length=81195210>
		##contig=<ID=chr18,length=78077248>
		##contig=<ID=chr19,length=59128983>
		##contig=<ID=chr20,length=63025520>
		##contig=<ID=chr21,length=48129895>
		##contig=<ID=chr22,length=51304566>
		##contig=<ID=chrM,length=16571>
		##contig=<ID=chrX,length=155270560>
		##contig=<ID=chrY,length=59373566>
		##source=HaplotypeCaller
		#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA19240chr2122
		chr21	9452821	.	T	<NON_REF>	.	.	END=9452821	GT:DP:GQ:MIN_DP:PL	0/0:7:12:7:0,12,180
		chr21	9452822	.	T	<NON_REF>	.	.	END=9452826	GT:DP:GQ:MIN_DP:PL	0/0:7:9:7:0,9,135
		chr21	9452827	.	A	<NON_REF>	.	.	END=9452831	GT:DP:GQ:MIN_DP:PL	0/0:8:12:8:0,12,180
		chr21	9452832	.	TATC	T,<NON_REF>	76.73	.	AS_RAW_BaseQRankSum=|1.5,1|NaN;AS_RAW_MQ=11797.00|5840.00|0.00;AS_RAW_MQRankSum=|-1.4,1|NaN;AS_RAW_ReadPosRankSum=|1.0,1|NaN;AS_SB_TABLE=0,4|0,3|0,0;BaseQRankSum=1.579;DP=8;ExcessHet=3.0103;MLEAC=1,0;MLEAF=0.500,0.00;MQRankSum=-1.368;NDA=2;RAW_MQ=20038.00;ReadPosRankSum=1.006	GT:AD:GQ:PL:SB	0/1:4,3,0:99:114,0,159,126,168,294:0,4,0,3
		chr21	9452836	.	A	<NON_REF>	.	.	END=9452836	GT:DP:GQ:MIN_DP:PL	0/0:7:0:7:0,0,28
		chr21	9452837	.	T	<NON_REF>	.	.	END=9452838	GT:DP:GQ:MIN_DP:PL	0/0:7:21:7:0,21,296
	}]\n
	file_write tmp/expected.tsv [deindent {
		chromosome	begin	end	type	ref	alt	name	quality	filter	alleleSeq1	alleleSeq2	zyg	phased	genotypes	alleledepth_ref	alleledepth	coverage	genoqual	MIN_DP	PGT	PID	PL	SB	AS_InbreedingCoeff	AS_QD	AS_RAW_BaseQRankSum	AS_RAW_MQ	AS_RAW_MQRankSum	AS_RAW_ReadPosRankSum	AS_SB_TABLE	BaseQRankSum	totalcoverage	DS	ExcessHet	InbreedingCoeff	MLEAC	MLEAF	MQ	MQRankSum	NDA	RAW_MQ	ReadPosRankSum
		chr21	9452820	9452821	snp	T	.	.	.	.	T	T	r	0	0;0			7	12	7			0,12,180																				
		chr21	9452821	9452826	ref	5	.	.	.	.	5	5	r	0	0;0			7	9	7			0,9,135																				
		chr21	9452826	9452831	ref	5	.	.	.	.	5	5	r	0	0;0			8	12	8			0,12,180																				
		chr21	9452832	9452835	del	ATC		.	76.73	.	ATC		t	0	0;1	4	3		99				114,0,159,126,168,294	0,4,0,3			|1.5,1|NaN	11797.00|5840.00|0.00	|-1.4,1|NaN	|1.0,1|NaN	0,4|0,3|0,0	1.579	8		3.0103		1	0.500		-1.368	2	20038.00	1.006
		chr21	9452835	9452836	snp	A	.	.	.	.	A	A	r	0	0;0			7	0	7			0,0,28																				
		chr21	9452836	9452838	ref	2	.	.	.	.	2	2	r	0	0;0			7	21	7			0,21,296																				
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

test vcf2tsv {gvcf BP_RESOLUTION} {
	cg gatk_gatk_genotypevcfs -dbdir $::refseqdir/hg19 data/varall-gatkh-bwa-sample1.gvcf tmp/test.vcf
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

test vcf2tsv {vcf2tsv -sort 0 all vcf files} {
	foreach file [glob data/*vcf] {
		cg vcf2tsv -split 1 -sort 0 $file tmp/test.tsv
		cg vcf2tsv -split 1 $file tmp/expected.tsv
		exec diff tmp/test.tsv tmp/expected.tsv
	}
} {}

test vcfcat {vcfcat basic} {
	write_vcf tmp/temp1.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
		20	14370	rs6054257	g	a	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
	}
	write_vcf tmp/temp2.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
		21	14380	x	g	a	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
	}
	write_vcf tmp/expected.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
		20	14370	rs6054257	g	a	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
		21	14380	x	g	a	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
	}
	cg vcfcat tmp/temp1.vcf tmp/temp2.vcf > tmp/temp3.vcf
	exec diff tmp/temp3.vcf tmp/expected.vcf
} {}

test vcfcat {vcfcat -o} {
	write_vcf tmp/temp1.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
		20	14370	rs6054257	g	a	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
	}
	write_vcf tmp/temp2.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
		21	14380	x	g	a	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
	}
	write_vcf tmp/expected.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
		20	14370	rs6054257	g	a	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
		21	14380	x	g	a	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
	}
	cg vcfcat -o tmp/temp3.vcf.gz tmp/temp1.vcf tmp/temp2.vcf
	exec gunzip tmp/temp3.vcf.gz
	exec diff tmp/temp3.vcf tmp/expected.vcf
} {}

testsummarize
