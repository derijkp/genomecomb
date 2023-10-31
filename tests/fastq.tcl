#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

proc checksplitoutput {numseq {maxparts 1000000}} {
	set part 1
	for {set pos 0} {$pos < 100} {incr pos $numseq} {
		putslog "checking $pos"
		exec zcat tmp/split/p${part}_seq_R1.fq.gz > tmp/test.fastq
		set pipe "zcat tmp/seq_R1.fq.gz"
		if {$pos > 0} {append pipe  " | tail -n +[expr {4*$pos+1}]"}
		if {$part < $maxparts} {
			append pipe " | head -[expr {4*$numseq}]"
		}
		append pipe " > tmp/expected.fastq"
		exec {*}$pipe
		exec diff tmp/test.fastq tmp/expected.fastq
		if {$part >= $maxparts} break
		incr part
	}
}

test fastq_split {fastq_split basic} {
	file copy {*}[glob data/seq_R*.fq.gz] tmp/
	set numseq 50
	cg fastq_split -numseq $numseq tmp/seq_R1.fq.gz tmp/split/seq_R1.fq.gz
	checksplitoutput $numseq
	bsort [glob tmp/split/*]
} {tmp/split/p1_seq_R1.fq.gz tmp/split/p2_seq_R1.fq.gz}

test fastq_split {fastq_split remainder -numseq 40} {
	file copy {*}[glob data/seq_R*.fq.gz] tmp/
	set numseq 40
	cg fastq_split -numseq $numseq tmp/seq_R1.fq.gz tmp/split/seq_R1.fq.gz
	checksplitoutput $numseq
	bsort [glob tmp/split/*]
} {tmp/split/p1_seq_R1.fq.gz tmp/split/p2_seq_R1.fq.gz tmp/split/p3_seq_R1.fq.gz}

test fastq_split {fastq_split remainder -numseq 33} {
	file copy {*}[glob data/seq_R*.fq.gz] tmp/
	set numseq 40
	cg fastq_split -numseq $numseq tmp/seq_R1.fq.gz tmp/split/seq_R1.fq.gz
	checksplitoutput $numseq
	bsort [glob tmp/split/*]
} {tmp/split/p1_seq_R1.fq.gz tmp/split/p2_seq_R1.fq.gz tmp/split/p3_seq_R1.fq.gz}

test fastq_split {fastq_split -numseq 99} {
	file copy {*}[glob data/seq_R*.fq.gz] tmp/
	set numseq 99
	cg fastq_split -numseq $numseq tmp/seq_R1.fq.gz tmp/split/seq_R1.fq.gz
	checksplitoutput $numseq
	bsort [glob tmp/split/*]
} {tmp/split/p1_seq_R1.fq.gz tmp/split/p2_seq_R1.fq.gz}

test fastq_split {fastq_split -numseq 25} {
	file copy {*}[glob data/seq_R*.fq.gz] tmp/
	set numseq 25
	cg fastq_split -numseq $numseq tmp/seq_R1.fq.gz tmp/split/seq_R1.fq.gz
	checksplitoutput $numseq
	bsort [glob tmp/split/*]
} {tmp/split/p1_seq_R1.fq.gz tmp/split/p2_seq_R1.fq.gz tmp/split/p3_seq_R1.fq.gz tmp/split/p4_seq_R1.fq.gz}

test fastq_split {fastq_split -maxparts} {
	file copy {*}[glob data/seq_R*.fq.gz] tmp/
	set numseq 25
	cg fastq_split -numseq $numseq -maxparts 2 tmp/seq_R1.fq.gz tmp/split/seq_R1.fq.gz
	checksplitoutput 25 2
	bsort [glob tmp/split/*]
} {tmp/split/p1_seq_R1.fq.gz tmp/split/p2_seq_R1.fq.gz}

test fastq_split {fastq_split -parts 4} {
	file copy {*}[glob data/seq_R*.fq.gz] tmp/
	set numseq 25
	cg fastq_split -parts 4 tmp/seq_R1.fq.gz tmp/split/seq_R1.fq.gz
	checksplitoutput $numseq
	bsort [glob tmp/split/*]
} {tmp/split/p1_seq_R1.fq.gz tmp/split/p2_seq_R1.fq.gz tmp/split/p3_seq_R1.fq.gz tmp/split/p4_seq_R1.fq.gz}

test fastq_split {fastq_split -parts 3} {
	file copy {*}[glob data/seq_R*.fq.gz] tmp/
	exec cg fastq_split -parts 3 tmp/seq_R1.fq.gz tmp/split/seq_R1.fq.gz >@ stdout
	checksplitoutput 33 3
	bsort [glob tmp/split/*]
} {tmp/split/p1_seq_R1.fq.gz tmp/split/p2_seq_R1.fq.gz tmp/split/p3_seq_R1.fq.gz}

test fastq2tsv {fastq2tsv} {
	cg fastq2tsv data/seq_R1.fq.gz tmp/test.tsv
} {}

test fastq2tsv {fastq2tsv 2} {
	file delete tmp/test.tsv
	cg fastq2tsv data/seq_R1.fq.gz tmp/test.tsv
	file delete tmp/test.fq.gz
	cg tsv2fastq tmp/test.tsv tmp/test.fq.gz
	cg zcat data/seq_R1.fq.gz > tmp/seq_R1.fq
	cg zcat tmp/test.fq.gz > tmp/test.fq
	exec diff tmp/seq_R1.fq tmp/test.fq
} {}

test fastq2tsv {fastq2tsv comments, pipeing} {
	file_write tmp/test.fastq [deindent {
		@SRR792091.9203/1	RG:Z:test
		CCAGTCCGAGAGCTCTGCCGGGTCCTGGGCACCAGGGGCCACTCTGCTGTGTCCCTTTGCTCCAGTCCCAATGTCACTCGCATGGTGCAGGGTCATGGCT
		+
		BCCFFFFFHHHHHJJJJJJIIJGHIJJJJIJJJJJIJJIIJJJJJJJHHHHEHFFFFFFEEEEEDCDDDDDDDDEDDDDDDDDDBCCCDDDD?CDDDDDD
		@SRR792091.86979/1	RG:Z:test
		TAATTTCTTTCTTGGTTGCCTCTGCCAAGCTCAGCCTGTTTGTGATCTCCTGCTGGCCCCGGCACTTCACAATCATGTAGCCCTTCTTGAGGGGGTACGT
		+
		@@@FDBDEHHHGHICIJEGGHIJIIJJIFCHGEGIGEG@GHIG@FBBF@DGHGGIGECDHHBBEFFFFFEECEEDDDEDDCCCDDDDDDDD<@9@99A>?
	}]\n
	file_write tmp/expected_test.tsv [deindent {
		name	sequence	quality	comments
		@SRR792091.9203/1	CCAGTCCGAGAGCTCTGCCGGGTCCTGGGCACCAGGGGCCACTCTGCTGTGTCCCTTTGCTCCAGTCCCAATGTCACTCGCATGGTGCAGGGTCATGGCT	BCCFFFFFHHHHHJJJJJJIIJGHIJJJJIJJJJJIJJIIJJJJJJJHHHHEHFFFFFFEEEEEDCDDDDDDDDEDDDDDDDDDBCCCDDDD?CDDDDDD	RG:Z:test
		@SRR792091.86979/1	TAATTTCTTTCTTGGTTGCCTCTGCCAAGCTCAGCCTGTTTGTGATCTCCTGCTGGCCCCGGCACTTCACAATCATGTAGCCCTTCTTGAGGGGGTACGT	@@@FDBDEHHHGHICIJEGGHIJIIJJIFCHGEGIGEG@GHIG@FBBF@DGHGGIGECDHHBBEFFFFFEECEEDDDEDDCCCDDDDDDDD<@9@99A>?	RG:Z:test
	}]\n
	file delete tmp/test.tsv
	cg fastq2tsv tmp/test.fastq tmp/test.tsv
	exec diff tmp/test.tsv tmp/expected_test.tsv
	# pipes
	exec cg fastq2tsv < tmp/test.fastq > tmp/test.tsv
	exec diff tmp/test.tsv tmp/expected_test.tsv
	# 2 input
	cg cat tmp/expected_test.tsv tmp/expected_test.tsv > tmp/expected_test2.tsv
	file delete tmp/test.tsv
	exec cg fastq2tsv tmp/test.fastq tmp/test.fastq tmp/test.tsv
	exec diff tmp/test.tsv tmp/expected_test2.tsv
} {}

test fastq2tsv {fastq2tsv comments and -fields} {
	file_write tmp/test.fastq [deindent {
		@SRR792091.9203/1	RG:Z:test
		CCAGTCCGAGAGCTCTGCCGGGTCCTGGGCACCAGGGGCCACTCTGCTGTGTCCCTTTGCTCCAGTCCCAATGTCACTCGCATGGTGCAGGGTCATGGCT
		+
		BCCFFFFFHHHHHJJJJJJIIJGHIJJJJIJJJJJIJJIIJJJJJJJHHHHEHFFFFFFEEEEEDCDDDDDDDDEDDDDDDDDDBCCCDDDD?CDDDDDD
		@SRR792091.86979/1	RG:Z:test	MM:Z:test2
		TAATTTCTTTCTTGGTTGCCTCTGCCAAGCTCAGCCTGTTTGTGATCTCCTGCTGGCCCCGGCACTTCACAATCATGTAGCCCTTCTTGAGGGGGTACGT
		+
		@@@FDBDEHHHGHICIJEGGHIJIIJJIFCHGEGIGEG@GHIG@FBBF@DGHGGIGECDHHBBEFFFFFEECEEDDDEDDCCCDDDDDDDD<@9@99A>?
	}]\n
	file_write tmp/expected_test.tsv [deindent {
		name	sequence	RG	comments
		@SRR792091.9203/1	CCAGTCCGAGAGCTCTGCCGGGTCCTGGGCACCAGGGGCCACTCTGCTGTGTCCCTTTGCTCCAGTCCCAATGTCACTCGCATGGTGCAGGGTCATGGCT	test	
		@SRR792091.86979/1	TAATTTCTTTCTTGGTTGCCTCTGCCAAGCTCAGCCTGTTTGTGATCTCCTGCTGGCCCCGGCACTTCACAATCATGTAGCCCTTCTTGAGGGGGTACGT	test	MM:Z:test2
	}]\n
	file delete tmp/test.tsv
	cg fastq2tsv -fields {name sequence RG comments} tmp/test.fastq tmp/test.tsv
	exec diff tmp/test.tsv tmp/expected_test.tsv
} {}

testsummarize
