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

testsummarize
