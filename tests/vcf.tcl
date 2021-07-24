#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

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

test vcfcat {vcfcat first empy} {
	file_write tmp/temp1.vcf {}
	write_vcf tmp/temp2.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
		21	14380	x	g	a	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
	}
	write_vcf tmp/expected.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
		21	14380	x	g	a	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
	}
	cg vcfcat tmp/temp1.vcf tmp/temp2.vcf > tmp/temp3.vcf
	exec diff tmp/temp3.vcf tmp/expected.vcf
} {}

test vcfcat {vcfcat -s 1 -o} {
	write_vcf tmp/temp1.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
		21	14370	rs6054257	g	a	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
	}
	write_vcf tmp/temp2.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
		20	14380	x	g	a	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
		20	14080	x	g	b	30	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
	}
	write_vcf tmp/expected.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
		20	14080	x	g	b	30	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
		20	14380	x	g	a	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
		21	14370	rs6054257	g	a	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
	}
	cg vcfcat -s 1 -o tmp/temp3.vcf.gz tmp/temp1.vcf tmp/temp2.vcf
	file delete tmp/temp3.vcf
	exec gunzip tmp/temp3.vcf.gz
	exec diff tmp/temp3.vcf tmp/expected.vcf
} {}

test vcfcat {vcfcat -sample} {
	write_vcf tmp/temp1.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      test
		20	14370	rs6054257	g	a	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51
	}
	write_vcf tmp/temp2.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      test
		21	14380	x	g	a	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51
	}
	write_vcf tmp/expected.vcf {
		CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      renamed
		20	14370	rs6054257	g	a	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51
		21	14380	x	g	a	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51
	}
	cg vcfcat -sample renamed tmp/temp1.vcf tmp/temp2.vcf > tmp/temp3.vcf
	exec diff tmp/temp3.vcf tmp/expected.vcf
} {}

test sortvcf {sortvcf basic} {
	file_write tmp/temp.vcf [deindent {
		##fileformat=VCFv4.0
		##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
		##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
		##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
		##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
		##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
		##contig=<ID=chrM,length=16571>
		##contig=<ID=chr1,length=249250621>
		##contig=<ID=chr2,length=243199373>
		##contig=<ID=chr3,length=198022430>
		#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001
		M	17330	.	T	A	3	q10	DP=11	GT:GQ:DP	0|0:49:3
		1	1110696	rs6040355	A	G,T	67	PASS	DP=10	GT:GQ:DP	1|2:21:6
		2	1230237	.	T	.	47	PASS	DP=13	GT:GQ:DP	0|0:54:7
	}]\n
	file_write tmp/expected.vcf [deindent {
		##fileformat=VCFv4.0
		##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
		##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
		##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
		##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
		##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
		##contig=<ID=chr1,length=249250621>
		##contig=<ID=chr2,length=243199373>
		##contig=<ID=chr3,length=198022430>
		##contig=<ID=chrM,length=16571>
		#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001
		1	1110696	rs6040355	A	G,T	67	PASS	DP=10	GT:GQ:DP	1|2:21:6
		2	1230237	.	T	.	47	PASS	DP=13	GT:GQ:DP	0|0:54:7
		M	17330	.	T	A	3	q10	DP=11	GT:GQ:DP	0|0:49:3
	}]\n
	exec cg sortvcf -stack 1 tmp/temp.vcf tmp/sorted.vcf
	exec diff tmp/sorted.vcf tmp/expected.vcf
} {}

test sortvcf {sortvcf compressed} {
	file_write tmp/temp.vcf [deindent {
		##fileformat=VCFv4.0
		##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
		##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
		##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
		##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
		##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
		##contig=<ID=chrM,length=16571>
		##contig=<ID=chr1,length=249250621>
		##contig=<ID=chr2,length=243199373>
		##contig=<ID=chr3,length=198022430>
		#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001
		M	17330	.	T	A	3	q10	DP=11	GT:GQ:DP	0|0:49:3
		1	1110696	rs6040355	A	G,T	67	PASS	DP=10	GT:GQ:DP	1|2:21:6
		2	1230237	.	T	.	47	PASS	DP=13	GT:GQ:DP	0|0:54:7
	}]\n
	file_write tmp/expected.vcf [deindent {
		##fileformat=VCFv4.0
		##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
		##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
		##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
		##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
		##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
		##contig=<ID=chr1,length=249250621>
		##contig=<ID=chr2,length=243199373>
		##contig=<ID=chr3,length=198022430>
		##contig=<ID=chrM,length=16571>
		#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001
		1	1110696	rs6040355	A	G,T	67	PASS	DP=10	GT:GQ:DP	1|2:21:6
		2	1230237	.	T	.	47	PASS	DP=13	GT:GQ:DP	0|0:54:7
		M	17330	.	T	A	3	q10	DP=11	GT:GQ:DP	0|0:49:3
	}]\n
	exec cg sortvcf -stack 1 tmp/temp.vcf tmp/sorted.vcf.gz
	exec gunzip tmp/sorted.vcf.gz
	exec diff tmp/sorted.vcf tmp/expected.vcf
} {}

testsummarize
