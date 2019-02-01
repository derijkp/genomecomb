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

testsummarize
