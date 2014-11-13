#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test genome_seq {basic} {
	exec cg genome_seq -i name data/reg_genome_seq.tsv /complgen/refseq/hg19_test > tmp/temp.fas  2> /dev/null
	exec diff tmp/temp.fas data/expected-reg_genome_seq.fas
} {} 

test genome_seq {basic -f -1 -n -1} {
	exec cg genome_seq -i name -f -1 -n -1 data/reg_genome_seq.tsv /complgen/refseq/hg19_test > tmp/temp.fas  2> /dev/null
	exec diff tmp/temp.fas data/expected-reg_genome_seq_fullmasked.fas
} {} 

test genome_seq {basic -p {snp 1000g}} {
	exec cg genome_seq -i name -f -1 -n -1 -p {snp 1000g} data/reg_genome_seq.tsv /complgen/refseq/hg19_test > tmp/temp.fas  2> /dev/null
	exec diff tmp/temp.fas data/expected-reg_genome_seq_fullmasked.fas
} {} 

test genome_seq {basic -p {1000g}} {
	exec cg genome_seq -i name -f -1 -n -1 -p {1000g} data/reg_genome_seq2.tsv /complgen/refseq/hg19_test > tmp/temp.fas  2> /dev/null
	exec diff tmp/temp.fas data/expected-reg_genome_seq2.fas
} {} 

test genome_seq {outfile} {
	exec cg genome_seq -i name data/reg_genome_seq.tsv /complgen/refseq/hg19_test tmp/temp.fas  2> /dev/null
	exec diff tmp/temp.fas data/expected-reg_genome_seq.fas
} {} 

test genome_seq {concat and makemap} {
	exec cg genome_seq -c "" -m tmp/temp.map -i name data/reg_genome_seq.tsv /complgen/refseq/hg19_test > tmp/temp.fas  2> /dev/null
	exec diff tmp/temp.fas data/expected-reg_genome_seq.ifas
	exec diff tmp/temp.map data/expected-reg_genome_seq.map
} {} 

test genome_seq {gcsplit} {
	exec cg genome_seq -g 50 -gs 60 -i name data/reg_genome_seq.tsv /complgen/refseq/hg19_test tmp/temp.fas  2> /dev/null
	exec diff tmp/temp-lowgc.fas data/expected-reg_genome_seq-lowgc.fas
	exec diff tmp/temp-highgc.fas data/expected-reg_genome_seq-highgc.fas
} {} 

test genome_seq {-g 0} {
	exec cg genome_seq -g 0 data/reg_genome_seq.tsv /complgen/refseq/hg19_test tmp/temp.fas  2> /dev/null
	exec diff tmp/temp.fas data/expected-reg_genome_seq-noid.fas
} {1c1
< >chr1-69328-69576 GC:48.0
---
> >chr1-69328-69576
3c3
< >chr1-69386-69634 GC:47.2
---
> >chr1-69386-69634
5c5
< >chr1-231830367-231830615 GC:59.7
---
> >chr1-231830367-231830615
7c7
< >chr2-145425276-145425524 GC:46.4
---
> >chr2-145425276-145425524
9c9
< >chr2-145427578-145427826 GC:37.1
---
> >chr2-145427578-145427826
11c11
< >chr2-145671229-145671490 GC:41.4
---
> >chr2-145671229-145671490
13c13
< >chr2-177497692-177497940 GC:31.5
---
> >chr2-177497692-177497940
15c15
< >chr2-178087040-178087288 GC:41.1
---
> >chr2-178087040-178087288
17c17
< >chr4-2069607-2069855 GC:54.0
---
> >chr4-2069607-2069855
19c19
< >chr4-5526991-5527238 GC:38.9
---
> >chr4-5526991-5527238
21c21
< >chr4-24981564-24981812 GC:36.3
---
> >chr4-24981564-24981812
23c23
< >chr7-26416288-26416536 GC:69.8
---
> >chr7-26416288-26416536
25c25
< >chr8-1949001-1949249 GC:56.5
---
> >chr8-1949001-1949249
27c27
< >chr18-103887-104141 GC:47.2
---
> >chr18-103887-104141
29c29
< >chr22-16108815-16109063 GC:41.1
---
> >chr22-16108815-16109063
31c31
< >chrM-16047-16295 GC:46.8
---
> >chrM-16047-16295
33c33
< >chrX-302322-302570 GC:65.3
---
> >chrX-302322-302570
}

test genome_seq {-g 0 -gd 0} {
	exec cg genome_seq -g 0 -gd 0 data/reg_genome_seq.tsv /complgen/refseq/hg19_test tmp/temp.fas  2> /dev/null
	exec diff tmp/temp.fas data/expected-reg_genome_seq-noid.fas
} {}

test genome_seq {gcsplit and -l} {
	exec cg genome_seq -g 50 -gs 60 -l _ -i name data/reg_genome_seq.tsv /complgen/refseq/hg19_test tmp/temp.fas  2> /dev/null
	exec diff tmp/temp-highgc.fas data/expected-reg_genome_seq-highgc-l.fas
} {} 

test genome_seq {split} {
	file delete {*}[glob tmp/*.fas]
	exec cg genome_seq -g 50 -s 1 -i name data/reg_genome_seq.tsv /complgen/refseq/hg19_test tmp/  2> /dev/null
	exec diff tmp/r1.fas data/expected-reg_genome_seq-1.fas
	exec diff tmp/r9.fas data/expected-reg_genome_seq-9.fas
	lsort -dict [glob tmp/*.fas]
} {tmp/r1.fas tmp/r2.fas tmp/r3.fas tmp/r4.fas tmp/r5.fas tmp/r6.fas tmp/r7.fas tmp/r8.fas tmp/r9.fas tmp/r10.fas tmp/r11.fas tmp/r12.fas tmp/r13.fas tmp/r14.fas tmp/r15.fas tmp/r16.fas tmp/r17.fas} 

test genome_seq {gcsplit split} {
	file delete -force {*}[glob tmp/*]
	exec cg genome_seq -g 50 -s 1 -gs 60 -i name data/reg_genome_seq.tsv /complgen/refseq/hg19_test tmp/temp.fas  2> /dev/null
	exec diff tmp/temp-lowgc/r1.fas data/expected-reg_genome_seq-1.fas
	exec diff tmp/temp-highgc/r9.fas data/expected-reg_genome_seq-9.fas
	list [lsort -dict [glob tmp/temp-lowgc/*.fas]] [lsort -dict [glob tmp/temp-highgc/*.fas]]
} {{tmp/temp-lowgc/r1.fas tmp/temp-lowgc/r2.fas tmp/temp-lowgc/r4.fas tmp/temp-lowgc/r5.fas tmp/temp-lowgc/r6.fas tmp/temp-lowgc/r7.fas tmp/temp-lowgc/r8.fas tmp/temp-lowgc/r10.fas tmp/temp-lowgc/r11.fas tmp/temp-lowgc/r15.fas tmp/temp-lowgc/r16.fas} {tmp/temp-highgc/r3.fas tmp/temp-highgc/r9.fas tmp/temp-highgc/r12.fas tmp/temp-highgc/r13.fas tmp/temp-highgc/r14.fas tmp/temp-highgc/r17.fas}} 

test primercheck {basic} {
	exec cg primercheck data/primers.tsv /complgen/refseq/hg19_test tmp/primerinfo.tsv 2> /dev/null 2> /dev/null
	exec diff tmp/primerinfo.tsv data/primercheck-results.tsv
} {} 

file delete -force tmp/temp.sft

set ::env(PATH) $keeppath

testsummarize
