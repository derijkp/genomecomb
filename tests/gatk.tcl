#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test gatkexec {version} {
	regexp {^4\.[0-9]\.[0-9]\.[0-9]$} [gatkexec version]
} 1

test gatkexec {error} {
	gatkexec {-XX:ParallelGCThreads=1 -d64 -Xms512m -Xmx4g} HaplotypeCaller -R bla -I bla -O bla
} {*The specified fasta file * does not exist.*} error match

test gatkexec {error -finishedpattern} {
	catch {
		gatkexec -finishedpattern {does not exist} {-XX:ParallelGCThreads=1 -d64 -Xms512m -Xmx4g} HaplotypeCaller -R bla -I bla -O bla
	} msg
} 0

test gatk {gatk_compar} {
	file copy data/varall-gatkh-bwa-sample1.gvcf tmp
	file copy data/varall-gatkh-bwa-sample2.gvcf tmp
	exec cg gatk_compar --stack 1 -dbdir $::refseqdir/hg19 -maxmem 4 -vqsr {} tmp/result.vcf.gz {*}[glob tmp/*.gvcf.gz] >@ stdout 2>@ stderr
	cg tsvdiff tmp/result.vcf.gz data/gatk_compar_expected.vcf.gz
} {}

# differs from using genomicsdbimport ...
#test gatk {gatk_compar -usecombinegvcfs 1} {
#	file copy data/varall-gatkh-bwa-sample1.gvcf tmp
#	file copy data/varall-gatkh-bwa-sample2.gvcf tmp
#	exec cg gatk_compar --stack 1 -usecombinegvcfs 1 -dbdir $::refseqdir/hg19 -maxmem 4 -vqsr {} tmp/result.vcf.gz {*}[glob tmp/*.gvcf.gz] >@ stdout 2>@ stderr
#	cg tsvdiff tmp/result.vcf.gz data/gatk_compar_expected.vcf.gz
#} {}

testsummarize

