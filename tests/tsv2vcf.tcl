#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test tsv2vcf {basic} {
	exec cg vcf2tsv data/test3.vcf tmp/temp.tsv
	exec cg tsv2vcf -dbdir $::refseqdir/hg19 tmp/temp.tsv tmp/temp.vcf
	exec diff tmp/temp.vcf data/test3.vcf
} {20,21d19
< ##split=1
< ##info=tsv converted from vcf
child process exited abnormally} error

test tsv2vcf {no metadata in comments} {
	exec cg vcf2tsv data/test3.vcf tmp/temp.tsv.temp
	cg select -rc 1 tmp/temp.tsv.temp tmp/temp.tsv
	exec cg tsv2vcf -dbdir $::refseqdir/hg19 tmp/temp.tsv tmp/temp.vcf
	exec diff tmp/temp.vcf data/test3.vcf
} {1a2,6
> ##fileDate=20090805
> ##reference=hg19
> ##phasing=partial
> ##FILTER=<ID=q10,Description="Quality below 10">
> ##FILTER=<ID=s50,Description="Less than 50% of samples have data">
4c9
< ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth (filtered)">
---
> ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
child process exited abnormally} error

test tsv2vcf {gatkh results} {
	exec cg vcf2tsv data/test4.vcf tmp/temp.tsv
	exec cg tsv2vcf -dbdir $::refseqdir/hg19 tmp/temp.tsv tmp/temp.vcf
	exec diff tmp/temp.vcf data/test4.vcf
} {47,48d46
< ##split=1
< ##info=tsv converted from vcf
74d71
< ##samplename=NA19240m
child process exited abnormally} error

test tsv2vcf {gatk results} {
	exec cg vcf2tsv data/var-gatk-bwa-test.vcf tmp/temp.tsv
	exec cg tsv2vcf -dbdir $::refseqdir/hg19 tmp/temp.tsv | grep -v {##GATKCommandLine.UnifiedGenotyper\|##samplename\|##split\|##info} > tmp/temp.vcf
	exec grep -v {##GATKCommandLine.UnifiedGenotyper} data/var-gatk-bwa-test.vcf > tmp/expected.vcf
	exec diff tmp/temp.vcf tmp/expected.vcf
} {27,28c27
< ##INFO=<ID=RPA_ref,Number=1,Type=Integer,Description="reference only value of: Number of times tandem repeat unit is repeated, for each allele (including reference)">
< ##INFO=<ID=RPA,Number=A,Type=Integer,Description="alleles only values of: Number of times tandem repeat unit is repeated, for each allele (including reference)">
---
> ##INFO=<ID=RPA,Number=.,Type=Integer,Description="Number of times tandem repeat unit is repeated, for each allele (including reference)">
child process exited abnormally} error

test tsv2vcf {gatkh results, no metadata} {
	exec cg vcf2tsv data/test4.vcf tmp/temp.tsv.temp
	cg select -rc 1 tmp/temp.tsv.temp tmp/temp.tsv
	exec cg tsv2vcf -sample NA19240m -dbdir $::refseqdir/hg19 tmp/temp.tsv tmp/temp.vcf
	exec grep -v {##contig=\|##ALT=\|##source=\|##FILTER=} data/test4.vcf > tmp/expected.vcf
	exec diff tmp/temp.vcf tmp/expected.vcf
} {4c4
< ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth (filtered)">
---
> ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
29c29
< ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
---
> ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
child process exited abnormally} error

test tsv2vcf {sam results} {
	exec cg vcf2tsv data/var-sam-bwa-test.vcf tmp/temp.tsv
	exec cg tsv2vcf -dbdir $::refseqdir/hg19 tmp/temp.tsv | grep -v {##split\|##info=\|##samplename\|##bcftools_callCommand} > tmp/temp.vcf
	exec grep -v {##bcftools_callCommand} data/var-sam-bwa-test.vcf > tmp/expected.vcf
	exec diff tmp/temp.vcf tmp/expected.vcf
} {15c15
< ##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)">
---
> ##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">
child process exited abnormally} error

test tsv2vcf {sam results,no metadata} {
	exec cg vcf2tsv data/var-sam-bwa-test.vcf tmp/temp.tsv.temp
	cg select -overwrite 1 -rc 1 tmp/temp.tsv.temp tmp/temp.tsv
	exec cg tsv2vcf -sample NA19240m -dbdir $::refseqdir/hg19 tmp/temp.tsv tmp/temp.vcf
	exec grep -v {##reference\|##samtools\|##bcftools\|##contig=\|##ALT=\|##source=\|##FILTER=} data/var-sam-bwa-test.vcf > tmp/expected.vcf
	exec diff tmp/temp.vcf tmp/expected.vcf
} {3,4c3,4
< ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
< ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth (filtered)">
---
> ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
> ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of high-quality bases">
9,10c9,10
< ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
< ##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)">
---
> ##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
> ##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">
21c21
< ##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
---
> ##INFO=<ID=MQ,Number=1,Type=Integer,Description="Root-mean-square mapping quality of covering reads">
child process exited abnormally} error

test tsv2vcf {gatkh results} {
	exec cg vcf2tsv data/var-gatkh-bwa-test.vcf tmp/temp.tsv
	exec cg tsv2vcf -dbdir $::refseqdir/hg19 tmp/temp.tsv | grep -v {##GATKCommandLine\|##samplename\|##source=\|##split\|##info} > tmp/temp.vcf
	exec grep -v {##GATKCommandLine\|##source=} data/var-gatkh-bwa-test.vcf > tmp/expected.vcf
	exec diff tmp/temp.vcf tmp/expected.vcf
} {1a2
> ##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
44d44
< ##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
child process exited abnormally} error

test tsv2vcf {gatkh results,no metadata} {
	exec cg vcf2tsv data/var-gatkh-bwa-test.vcf tmp/temp.tsv.temp
	cg select -overwrite 1 -rc 1 tmp/temp.tsv.temp tmp/temp.tsv
	exec cg tsv2vcf -sample NA19240m -dbdir $::refseqdir/hg19 tmp/temp.tsv | grep -v {##GATKCommandLine\|##samplename\|##source=\|##split\|##info} > tmp/temp.vcf
	exec grep -v {##GATKCommandLine\|##source=\|##contig=\|##ALT=\|##source=\|##FILTER=} data/var-gatkh-bwa-test.vcf > tmp/expected.vcf
	exec diff tmp/temp.vcf tmp/expected.vcf
} {4c4
< ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth (filtered)">
---
> ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
9c9
< ##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase Set">
---
> ##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phasing set (typically the position of the first variant in the set)">
29c29
< ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
---
> ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
40c40
< ##INFO=<ID=RAW_MQandDP,Number=.,Type=String,Description="RAW_MQandDP, no further description found">
---
> ##INFO=<ID=RAW_MQandDP,Number=2,Type=Integer,Description="Raw data (sum of squared MQ and total depth) for improved RMS Mapping Quality calculation. Incompatible with deprecated RAW_MQ formulation.">
child process exited abnormally} error

test tsv2vcf {strelka results} {
	exec cg vcf2tsv data/var-strelka-bwa-test.vcf tmp/temp.tsv
	exec cg tsv2vcf -dbdir $::refseqdir/hg19 tmp/temp.tsv | grep -v {##GATKCommandLine\|##samplename\|##source=\|##split\|##info} > tmp/temp.vcf
	exec grep -v {##GATKCommandLine\|##source=} data/var-strelka-bwa-test.vcf > tmp/expected.vcf
	exec diff tmp/temp.vcf tmp/expected.vcf
} {1c1
< ##fileformat=VCFv4.2
---
> ##fileformat=VCFv4.1
48c48
< ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed. For indels this value only includes reads which confidently support each allele (posterior prob 0.51 or higher that read contains indicated allele vs all other intersecting indel alleles)">
---
> ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed. For indels this value only includes reads which confidently support each allele (posterior prob 0.51 or higher that read contains indicated allele vs all other intersecting indel alleles)">
55a56
> ##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the region described in this record">
child process exited abnormally} error

test tsv2vcf {strelka results,no metadata} {
	exec cg vcf2tsv data/var-strelka-bwa-test.vcf tmp/temp.tsv.temp
	cg select -rc 1 tmp/temp.tsv.temp tmp/temp.tsv
	exec cg tsv2vcf -sample NA19240m -dbdir $::refseqdir/hg19 tmp/temp.tsv tmp/temp.vcf
	exec grep -v {##content\|##startTime\|##fileDate=\|##reference=\|##cmdline=\|#Depth\|##contig=\|##ALT=\|##source\|##FILTER=} data/var-strelka-bwa-test.vcf > tmp/expected.vcf
	exec diff tmp/temp.vcf tmp/expected.vcf
} {1c1
< ##fileformat=VCFv4.2
---
> ##fileformat=VCFv4.1
5c5
< ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth (filtered)">
---
> ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Filtered basecall depth used for site genotyping. In a non-variant multi-site block this value represents the average of all sites in the block.">
8c8
< ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
---
> ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed. For indels this value only includes reads which confidently support each allele (posterior prob 0.51 or higher that read contains indicated allele vs all other intersecting indel alleles)">
12,13c12,13
< ##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase Set">
< ##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.">
---
> ##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">
> ##FORMAT=<ID=SB,Number=1,Type=Float,Description="Sample site strand bias">
15a16
> ##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the region described in this record">
22c23
< ##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
---
> ##INFO=<ID=MQ,Number=1,Type=Integer,Description="RMS of mapping quality">
child process exited abnormally} error

test tsv2vcf {old conversion} {
	file copy -force data/convertedvcf_old.tsv tmp/temp.tsv
	exec cg tsv2vcf -dbdir $::refseqdir/hg19 tmp/temp.tsv tmp/temp.vcf
	exec diff tmp/temp.vcf data/convertedvcf_old.vcf
} {}

test tsv2vcf {compressed} {
	exec cg lz4 < data/convertedvcf_old.tsv > tmp/temp.tsv.lz4
	exec cg tsv2vcf -dbdir $::refseqdir/hg19 tmp/temp.tsv.lz4 tmp/temp.vcf
	exec diff tmp/temp.vcf data/convertedvcf_old.vcf
} {}

test tsv2vcf {same variant twice in tsv (e.g. after liftover)} {
	file_write tmp/temp.tsv [deindent {
		chromosome	begin	end	type	ref	alt	name	quality	filter	alleleSeq1	alleleSeq2	zyg	phased	genotypes	alleledepth_ref	alleledepth	coverage	genoqual	PGT	PID	PL	RGQ	SB	allelecount	frequency	totalallelecount	AS_BaseQRankSum	AS_FS	AS_InbreedingCoeff	AS_MQ	AS_MQRankSum	AS_QD	AS_RAW_BaseQRankSum	AS_RAW_MQ	AS_RAW_MQRankSum	AS_RAW_ReadPosRankSum	AS_ReadPosRankSum	AS_SB_TABLE	AS_SOR	BaseQRankSum	ClippingRankSum	totalcoverage	DS	ExcessHet	FS	InbreedingCoeff	MLEAC	MLEAF	MQ	MQRankSum	NDA	QD	RAW_MQ	ReadPosRankSum	SOR
		chr21	42735718	42735719	snp	G	A	.	125.60	.	G	A	t	0	0;1	4	4	8	99			133,0,104			1	0.500	2	1.100	0.000		60.00	0.000	15.70					-0.400		0.693	1.10	0.00	8		3.0103	0.000		1	0.500	60.00	0.00	2	15.70		-3.660e-01	0.693
		chr21	42735718	42735719	snp	G	A	.	100.00	.	A	A	m	0	1;1	4	4	8	80			99,0,98			1	0.500	2	1.200	0.000		60.00	0.000	15.80					-0.500		0.793	1.20	0.10	8		3.1103	0.100		1	0.500	70.00	0.00	2	15.80		-3.760e-01	0.793
		chr21	42735718	42735719	snp	G	A	.	101.00	.	A	G	t	0	1;0	4	4	8	80			99,0,98			1	0.500	2	1.200	0.000		60.00	0.000	15.80					-0.500		0.793	1.20	0.10	8		3.1103	0.100		1	0.500	70.00	0.00	2	15.80		-3.760e-01	0.793
		chr21	42735718	42735719	snp	G	C	.	99.60	.	G	C	t	0	0;1	4	4	8	99			133,0,104			1	0.500	2	1.100	0.000		60.00	0.000	15.70					-0.400		0.693	1.10	0.00	8		3.0103	0.000		1	0.500	60.00	0.00	2	15.70		-3.660e-01	0.693
	}]\n
	exec cg tsv2vcf -dbdir $::refseqdir/hg19 tmp/temp.tsv tmp/temp.vcf
	# catchstderr_exec gatk ValidateVariants -R $::refseqdir/hg19/genome_hg19.fa -V tmp/temp.vcf --validation-type-to-exclude ALL
	exec tail -4 tmp/temp.vcf
} {#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	temp
chr21	42735719	.	G	A,C	125.60	.	AC=1,1;AF=0.500,0.500;AN=2;AS_BaseQRankSum=1.100,1.100;AS_FS=0.000,0.000;AS_MQ=60.00,60.00;AS_MQRankSum=0.000,0.000;AS_QD=15.70,15.70;AS_ReadPosRankSum=-0.400,-0.400;AS_SOR=0.693,0.693;BaseQRankSum=1.10;ClippingRankSum=0.00;DP=8;ExcessHet=3.0103;FS=0.000;MLEAC=1,1;MLEAF=0.500,0.500;MQ=60.00;MQRankSum=0.00;NDA=2;QD=15.70;ReadPosRankSum=-3.660e-01;SOR=0.693	GT:AD:DP:GQ:PL	0/2:4,4,4:8:99:133,0,104
chr21	42735719	.	G	A	101.00	.	AC=1;AF=0.500;AN=2;AS_BaseQRankSum=1.200;AS_FS=0.000;AS_MQ=60.00;AS_MQRankSum=0.000;AS_QD=15.80;AS_ReadPosRankSum=-0.500;AS_SOR=0.793;BaseQRankSum=1.20;ClippingRankSum=0.10;DP=8;ExcessHet=3.1103;FS=0.100;MLEAC=1;MLEAF=0.500;MQ=70.00;MQRankSum=0.00;NDA=2;QD=15.80;ReadPosRankSum=-3.760e-01;SOR=0.793	GT:AD:DP:GQ:PL	1/0:4,4:8:80:99,0,98
chr21	42735719	.	G	A	100.00	.	AC=1;AF=0.500;AN=2;AS_BaseQRankSum=1.200;AS_FS=0.000;AS_MQ=60.00;AS_MQRankSum=0.000;AS_QD=15.80;AS_ReadPosRankSum=-0.500;AS_SOR=0.793;BaseQRankSum=1.20;ClippingRankSum=0.10;DP=8;ExcessHet=3.1103;FS=0.100;MLEAC=1;MLEAF=0.500;MQ=70.00;MQRankSum=0.00;NDA=2;QD=15.80;ReadPosRankSum=-3.760e-01;SOR=0.793	GT:AD:DP:GQ:PL	1/1:4,4:8:80:99,0,98}


testsummarize
