#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test genome_seq {basic} {
	exec cg genome_seq -i name data/reg_genome_seq.tsv $::refseqdir/hg19 > tmp/temp.fas
	exec diff tmp/temp.fas data/expected-reg_genome_seq.fas
} {} 

test genome_seq {basic -f -1 -n -1} {
	exec cg genome_seq -i name -f -1 -n -1 data/reg_genome_seq.tsv $::refseqdir/hg19 > tmp/temp.fas
	exec diff tmp/temp.fas data/expected-reg_genome_seq_fullmasked.fas
} {} 

test genome_seq {basic -p {snp 1000g}} {
	exec cg genome_seq -i name -f -1 -n -1 -p {snp 1000g} data/reg_genome_seq.tsv $::refseqdir/hg19 > tmp/temp.fas
	exec diff tmp/temp.fas data/expected-reg_genome_seq_fullmasked.fas
} {} 

test genome_seq {basic -p {1000g}} {
	exec cg genome_seq -i name -f -1 -n -1 -p {1000g} data/reg_genome_seq2.tsv $::refseqdir/hg19 > tmp/temp.fas
	exec diff tmp/temp.fas data/expected-reg_genome_seq2.fas
} {} 

test genome_seq {outfile} {
	exec cg genome_seq -i name data/reg_genome_seq.tsv $::refseqdir/hg19 tmp/temp.fas
	exec diff tmp/temp.fas data/expected-reg_genome_seq.fas
} {} 

test genome_seq {-f 0.002 -n 0.4 freqp} {
	cg cplinked $::refseqdir/hg19 tmp
	file delete {*}[glob tmp/hg19/var_hg19_snp135.*]
	file copy data/var_hg19_partofsnp135.tsv.lz4 tmp/hg19/var_hg19_snp135.tsv.lz4
	cg lz4index tmp/hg19/var_hg19_snp135.tsv.lz4
	cg maketabix tmp/hg19/var_hg19_snp135.tsv.lz4
	exec cg genome_seq -i name -f 0.002 -n 0.4 data/reg_genome_seq.tsv tmp/hg19 > tmp/temp.fas
	exec diff tmp/temp.fas data/expected-reg_genome_seq_f0.002n0.4.fas
} {} 

test genome_seq {-f 0.002 -n 0.4 freqp} {
	cg cplinked $::refseqdir/hg19 tmp
	file delete {*}[glob tmp/hg19/var_hg19_snp135.*]
	cg select -f {chrom start end type ref alt name {freqp=catch($freq * 100.0,$freq)}} data/var_hg19_partofsnp135.tsv.lz4 tmp/hg19/var_hg19_snp135.tsv.lz4
	cg lz4index tmp/hg19/var_hg19_snp135.tsv.lz4
	cg maketabix tmp/hg19/var_hg19_snp135.tsv.lz4
	exec cg genome_seq -i name -f 0.002 -n 0.4 data/reg_genome_seq.tsv tmp/hg19 > tmp/temp.fas
	exec diff tmp/temp.fas data/expected-reg_genome_seq_f0.002n0.4.fas
} {} 

test genome_seq {no freq or freqp error} {
	cg cplinked $::refseqdir/hg19 tmp
	file delete {*}[glob tmp/hg19/var_hg19_snp135.*]
	cg select -f {chrom start end type ref alt name} data/var_hg19_partofsnp135.tsv.lz4 tmp/hg19/var_hg19_snp135.tsv.lz4
	cg lz4index tmp/hg19/var_hg19_snp135.tsv.lz4
	cg maketabix tmp/hg19/var_hg19_snp135.tsv.lz4
	exec cg genome_seq -i name -f 0.002 -n 0.4 data/reg_genome_seq.tsv tmp/hg19 > tmp/temp.fas
	exec diff tmp/temp.fas data/expected-reg_genome_seq_f0.002n0.4.fas
} {file *tmp/hg19/var_hg19_snp135.tsv.gz has no freq or freqp field} match error

test genome_seq {concat and makemap} {
	exec cg genome_seq -c "" -m tmp/temp.map -i name data/reg_genome_seq.tsv $::refseqdir/hg19 > tmp/temp.fas
	exec diff tmp/temp.fas data/expected-reg_genome_seq.ifas
	exec diff tmp/temp.map data/expected-reg_genome_seq.map
} {} 

test genome_seq {gcsplit} {
	exec cg genome_seq -g 50 -gs 60 -i name data/reg_genome_seq.tsv $::refseqdir/hg19 tmp/temp.fas
	exec diff tmp/temp-lowgc.fas data/expected-reg_genome_seq-lowgc.fas
	exec diff tmp/temp-highgc.fas data/expected-reg_genome_seq-highgc.fas
} {} 

test genome_seq {-g 0} {
	exec cg genome_seq -g 0 data/reg_genome_seq.tsv $::refseqdir/hg19 tmp/temp.fas
	exec diff tmp/temp.fas data/expected-reg_genome_seq-noid-g0.fas
} {}

test genome_seq {-g 0 -gd 0} {
	exec cg genome_seq -g 0 -gd 0 data/reg_genome_seq.tsv $::refseqdir/hg19 tmp/temp.fas
	exec diff tmp/temp.fas data/expected-reg_genome_seq-noid.fas
} {}

test genome_seq {gcsplit and -l} {
	exec cg genome_seq -g 50 -gs 60 -l _ -i name data/reg_genome_seq.tsv $::refseqdir/hg19 tmp/temp.fas
	exec diff tmp/temp-highgc.fas data/expected-reg_genome_seq-highgc-l.fas
} {} 

test genome_seq {split} {
	exec cg genome_seq -g 50 -s 1 -i name data/reg_genome_seq.tsv $::refseqdir/hg19 tmp/
	exec diff tmp/r1.fas data/expected-reg_genome_seq-1.fas
	exec diff tmp/r9.fas data/expected-reg_genome_seq-9.fas
	lsort -dict [glob tmp/*.fas]
} {tmp/r1.fas tmp/r2.fas tmp/r3.fas tmp/r4.fas tmp/r5.fas tmp/r6.fas tmp/r7.fas tmp/r8.fas tmp/r9.fas tmp/r10.fas tmp/r11.fas tmp/r12.fas tmp/r13.fas tmp/r14.fas tmp/r15.fas tmp/r16.fas tmp/r17.fas} 

test genome_seq {gcsplit split} {
	exec cg genome_seq -g 50 -s 1 -gs 60 -i name data/reg_genome_seq.tsv $::refseqdir/hg19 tmp/temp.fas
	exec diff tmp/temp-lowgc/r1.fas data/expected-reg_genome_seq-1.fas
	exec diff tmp/temp-highgc/r9.fas data/expected-reg_genome_seq-9.fas
	list [lsort -dict [glob tmp/temp-lowgc/*.fas]] [lsort -dict [glob tmp/temp-highgc/*.fas]]
} {{tmp/temp-lowgc/r1.fas tmp/temp-lowgc/r2.fas tmp/temp-lowgc/r4.fas tmp/temp-lowgc/r5.fas tmp/temp-lowgc/r6.fas tmp/temp-lowgc/r7.fas tmp/temp-lowgc/r8.fas tmp/temp-lowgc/r10.fas tmp/temp-lowgc/r11.fas tmp/temp-lowgc/r15.fas tmp/temp-lowgc/r16.fas} {tmp/temp-highgc/r3.fas tmp/temp-highgc/r9.fas tmp/temp-highgc/r12.fas tmp/temp-highgc/r13.fas tmp/temp-highgc/r14.fas tmp/temp-highgc/r17.fas}} 

test genome_seq {bugfix gcval error: -n 0.01 -p Common -r "N" -i "name" -gs 65 -gd 0} {
	exec cg genome_seq -n 0.01 -p Common -r "N" -i "name" -gs 60 -gd 0 data/reg_genome_seq.tsv $::refseqdir/hg19 tmp/temp.fas
} {}

test genome_seq {regionlist} {
	exec cg genome_seq {chr1:69328-69576 chr18:103887-104141} $::refseqdir/hg19 tmp/result.fas
	file_write tmp/expected.fas {>chr1-69328-69576
GCCAGCGCAAAGTCATCTCTTTCAAGGGCTGCCTTGTTCAGATATTTCTCCTTCACTTCTTTGGTGGGAGTGAGATGGTGATNCTCATAGCCATGGGCTNTGACAGATATATAGCAATATGCAANCCCCTACACTACACTACAATTAtGTGTGGCAACGCATGTGTCgGCATTATGGCTGTCACATGGGGAATTGGCTTTCTCCATTCGGTGAGCCAGTTGGCNTTTGCCGTGCACTTACTCTTCTGT
>chr18-103887-104141
ctccagaggcNgaggcaggagaatggtgtgaacctgggaggaggagcttgcagtgagccgggatcatgccactgcattccagcccgggcaacagagcaagactccatctcaaaaaaaaaaaaaaaaaaaaaGGATCTCTGTTTAGAATGCTACCTATTGCCTTCTGGATAGAATCACAACTCTTTACCACAATCAACACAGCTTCAgccctgcttctatatccagcctcatctatttctgctcctcctccttat
}
	exec diff tmp/result.fas tmp/expected.fas
} {}

test genome_seq {regionlist of one} {
	exec cg genome_seq {chr18:103887-104141} $::refseqdir/hg19 tmp/result.fas
	file_write tmp/expected.fas {>chr18-103887-104141
ctccagaggcNgaggcaggagaatggtgtgaacctgggaggaggagcttgcagtgagccgggatcatgccactgcattccagcccgggcaacagagcaagactccatctcaaaaaaaaaaaaaaaaaaaaaGGATCTCTGTTTAGAATGCTACCTATTGCCTTCTGGATAGAATCACAACTCTTTACCACAATCAACACAGCTTCAgccctgcttctatatccagcctcatctatttctgctcctcctccttat
}
	exec diff tmp/result.fas tmp/expected.fas
} {}

file delete -force tmp/temp.sft

set ::env(PATH) $keeppath

testsummarize
