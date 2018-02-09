#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test primercheck {basic} {
	exec cg primercheck data/primers.tsv $::refseqdir/hg19 tmp/primerinfo.tsv
	exec diff tmp/primerinfo.tsv data/primercheck-results.tsv
} {} 

test primercheck {basic} {
	write_tab tmp/primers.tsv {
		name	primer1	primer2	info
		disc1	CCCTCAACTTGTCACTTAAAGAAATCAC	AGTCAAACGAAAACTTCATAGGCTTCT	DISC1
		disc1extra	CCCTCAACTTGTCACTTAAAGAAATCAC	CCAGTCAAACGAAAACTTCATAGGCTTCT	DISC1extra
	}
	exec cg primercheck -i * tmp/primers.tsv $::refseqdir/hg19 tmp/primerinfo.tsv
	write_tab tmp/expected.tsv {
		chromosome	begin	end	name	outer_begin	outer_end	strand	primer1	primer2	numamplicons	amplicons	primer1_hits	primer1_snpsmaxfreq	primer1_snps	primer2_hits	primer2_snpsmaxfreq	primer2_snps	amplicon_fts	remarks	info
		1	231931035	231931071	disc1	231931007	231931098	+	CCCTCAACTTGTCACTTAAAGAAATCAC	AGTCAAACGAAAACTTCATAGGCTTCT	1	1:231931035-231931071	18	0.001826	rs76230451(snp@chr1:231931033-231931034;freq=0.001826;valid=by-cluster,by-frequency,by-1000genomes;submitterCount=2)	12	-	{}	{}	{}	DISC1
		1	231931035	231931071	disc1extra	?	?	+	CCCTCAACTTGTCACTTAAAGAAATCAC	CCAGTCAAACGAAAACTTCATAGGCTTCT	1	1:231931035-231931071	18	?	?	12	?	?	?	{Main target not found (no amplicon with full primer match)}	DISC1extra
	}
	exec diff tmp/primerinfo.tsv tmp/expected.tsv
} {} 

test primercheck {freq} {
	cg cplinked $::refseqdir/hg19 tmp
	file delete {*}[glob tmp/hg19/var_hg19_snp135.*]
	file copy data/var_hg19_partofsnp135.tsv.lz4 tmp/hg19/var_hg19_snp135.tsv.lz4
	cg lz4index tmp/hg19/var_hg19_snp135.tsv.lz4
	cg maketabix tmp/hg19/var_hg19_snp135.tsv.lz4
	exec cg primercheck data/primers.tsv tmp/hg19 tmp/primerinfo.tsv
	exec diff tmp/primerinfo.tsv data/primercheck-results.tsv
} {} 

test primercheck {freqp} {
	cg cplinked $::refseqdir/hg19 tmp
	file delete {*}[glob tmp/hg19/var_hg19_snp135.*]
	cg select -f {chrom start end type ref alt name {freqp=catch($freq * 100.0,$freq)} avHetSE strand molType valid func weight exceptions submitterCount submitters alleleFreqCount alleles alleleNs alleleFreqs bitfields} data/var_hg19_partofsnp135.tsv.lz4 tmp/hg19/var_hg19_snp135.tsv.lz4
	cg lz4index tmp/hg19/var_hg19_snp135.tsv.lz4
	cg maketabix tmp/hg19/var_hg19_snp135.tsv.lz4
	exec cg primercheck data/primers.tsv tmp/hg19 tmp/primerinfo.tsv
	exec diff tmp/primerinfo.tsv data/primercheck-results.tsv
} {} 

file delete -force tmp/temp.sft

set ::env(PATH) $keeppath

testsummarize
