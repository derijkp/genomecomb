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

test vcf2tsv {vcf2tsv} {
	exec cg vcf2tsv data/test.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test.vcf2tsv
} {}

test vcf2tsv {vcf2tsv ins and del} {
	exec cg vcf2tsv data/test2.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test2.vcf2tsv
} {}

test vcf2tsv {vcf2tsv 1000glow} {
	exec cg vcf2tsv data/test1000glow.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test1000glow.vcf2tsv
} {}

test vcf2tsv {vcf2tsv -split ori} {
	exec cg vcf2tsv -split ori data/test.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-testori.vcf2tsv
} {}

test vcf2tsv {vcf2tsv ins and del -split ori} {
	exec cg vcf2tsv -split ori data/test2.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test2ori.vcf2tsv
} {}

test vcf2tsv {vcf2tsv -split 1} {
	exec cg vcf2tsv -s 1 data/test.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-tests.vcf2tsv
} {}

test vcf2tsv {vcf2tsv ins and del -split 1} {
	exec cg vcf2tsv -s 1 data/test2.vcf tmp/temp.tsv
	exec diff tmp/temp.tsv data/expected-test2s.vcf2tsv
} {}

test vcf2tsv {vcf2tsv ins and del -split 1 -typelist .} {
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

test vcf2tsv {vcf2tsv space in name} {
	file copy -force data/test.vcf "tmp/test it.vcf"
	exec cg vcf2tsv "tmp/test it.vcf" "tmp/test it.tsv"
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

test vcf2tsv {svtest.vcf} {
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

testsummarize
