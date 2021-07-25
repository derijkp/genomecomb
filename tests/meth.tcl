#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl
set keepdir [pwd]

set test_cleantmp 0

# tests
# =====

lappend dopts -threads 1
set runopts {-stack 1}

# download testdata in process_ont.tcl

test meth {meth_nanopolish} {
	cd $::smalltestdir
	file delete -force tmp/meth_nanopolish
	file mkdir tmp/meth_nanopolish/fast5
	file mkdir tmp/meth_nanopolish/fastq
	foreach file [glob -nocomplain ori/ont_f5c_chr22_meth_example/fast5/*] {
		mklink $file tmp/meth_nanopolish/fast5/[file tail $file]
	}
	foreach file [glob -nocomplain ori/ont_f5c_chr22_meth_example/fastq/*] {
		mklink $file tmp/meth_nanopolish/fastq/[file tail $file]
	}
	mklink  ori/ont_f5c_chr22_meth_example/map-hlongshot-sminimap2-methtest.bam tmp/meth_nanopolish/map-hlongshot-sminimap2-methtest.bam
	mklink  ori/ont_f5c_chr22_meth_example/map-hlongshot-sminimap2-methtest.bam.bai tmp/meth_nanopolish/map-hlongshot-sminimap2-methtest.bam.bai

	puts time:[time {
		exec cg meth_nanopolish {*}$::dopts \
			-distrreg chr \
			-threads 6 \
			-refseq $::refseqdir/hg19 \
			tmp/meth_nanopolish/fast5 tmp/meth_nanopolish/fastq tmp/meth_nanopolish/map-hlongshot-sminimap2-methtest.bam \
			tmp/meth_nanopolish/meth-methtest.tsv.zst \
			>& tmp/meth_nanopolish/meth_nanopolish.log
	}]

	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 -x *.bam -x *.bai -x fastqc_report.html \
		-x colinfo -x meth.html -x *.zsti -x *.lz4i -x *.finished -x info_analysis.tsv \
		-x *.analysisinfo -x *.png -x *.submitting \
		-x *log_jobs -x *.index -x *.log \
		tmp/meth_nanopolish expected/meth_nanopolish]
	join [list_remove $result {}] \n
} {}

test meth {meth_nanopolish with gz results} {
	cd $::smalltestdir
	file delete -force tmp/methgz
	file mkdir tmp/methgz/samples/methtest/fast5
	file mkdir tmp/methgz/samples/methtest/fastq
	foreach file [glob -nocomplain ori/ont_f5c_chr22_meth_example/fast5/*] {
		mklink $file tmp/methgz/samples/methtest/fast5/[file tail $file]
	}
	foreach file [glob -nocomplain ori/ont_f5c_chr22_meth_example/fastq/*] {
		mklink $file tmp/methgz/samples/methtest/fastq/[file tail $file]
	}
	file mkdir tmp/methgz/samples/methtest2/fast5
	file mkdir tmp/methgz/samples/methtest2/fastq
	foreach file {ori/ont_f5c_chr22_meth_example/fast5/batch_0.fast5 ori/ont_f5c_chr22_meth_example/fast5/batch_1.fast5} {
		mklink $file tmp/methgz/samples/methtest2/fast5/[file tail $file]
	}
	foreach file {ori/ont_f5c_chr22_meth_example/fastq/batch_0.fastq.gz ori/ont_f5c_chr22_meth_example/fastq/batch_1.fastq.gz} {
		mklink $file tmp/methgz/samples/methtest2/fastq/[file tail $file]
	}
	# set ::env(CG_FAST_SHAREDSTORAGE) /buffer/tmp
	puts time:[time {
		exec cg process_project {*}$::runopts {*}$::dopts -split 1 \
			-paired 0 -clip 0 \
			-meth_nanopolish-compression gz \
			-maxfastqdistr 250 \
			-aligner {minimap2} \
			-removeduplicates 0 \
			-realign 0 \
			-distrreg chr \
			-svcallers {sniffles npinv} \
			-varcallers longshot \
			-methcallers nanopolish \
			-hap_bam 1 \
			-threads 6 \
			-reports {-fastqc predictgender} \
			tmp/methgz $::refseqdir/hg19 >& tmp/methgz.log
	}]
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 -x *.bam -x *.bai -x fastqc_report.html \
		-x colinfo -x meth.html -x *.zsti -x *.lz4i -x *.finished -x info_analysis.tsv \
		-x *.gz.tbi -x *sniffles*.vcf \
		-x *.analysisinfo -x *.png -x *.submitting \
		-x *log_jobs -x *.index \
		tmp/methgz expected/methgz]
	join [list_remove $result {}] \n
} {}

test meth {meth_nanopolish preset gpc} {
	cd $::smalltestdir
	file delete -force tmp/methgpc
	file mkdir tmp/methgpc/samples/methtest/fast5
	file mkdir tmp/methgpc/samples/methtest/fastq
	foreach file [glob -nocomplain ori/ont_f5c_chr22_meth_example/fast5/*] {
		mklink $file tmp/methgpc/samples/methtest/fast5/[file tail $file]
	}
	foreach file [glob -nocomplain ori/ont_f5c_chr22_meth_example/fastq/*] {
		mklink $file tmp/methgpc/samples/methtest/fastq/[file tail $file]
	}
	file mkdir tmp/methgpc/samples/methtest2/fast5
	file mkdir tmp/methgpc/samples/methtest2/fastq
	foreach file {ori/ont_f5c_chr22_meth_example/fast5/batch_0.fast5 ori/ont_f5c_chr22_meth_example/fast5/batch_1.fast5} {
		mklink $file tmp/methgpc/samples/methtest2/fast5/[file tail $file]
	}
	foreach file {ori/ont_f5c_chr22_meth_example/fastq/batch_0.fastq.gz ori/ont_f5c_chr22_meth_example/fastq/batch_1.fastq.gz} {
		mklink $file tmp/methgpc/samples/methtest2/fastq/[file tail $file]
	}
	# set ::env(CG_FAST_SHAREDSTORAGE) /buffer/tmp
	puts time:[time {
		exec cg process_project {*}$::runopts {*}$::dopts -split 1 \
			-paired 0 -clip 0 \
			-maxfastqdistr 250 \
			-aligner {minimap2} \
			-removeduplicates 0 \
			-realign 0 \
			-distrreg chr \
			-svcallers {} \
			-varcallers {} \
			-methcallers nanopolish_gpc \
			-hap_bam 1 \
			-threads 6 \
			-reports {-fastqc predictgender} \
			tmp/methgpc $::refseqdir/hg19 >& tmp/methgpc.log
	}]
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 -x *.bam -x *.bai -x fastqc_report.html \
		-x colinfo -x meth.html -x *.zsti -x *.lz4i -x *.finished -x info_analysis.tsv \
		-x *.analysisinfo -x *.png -x *.submitting \
		-x *log_jobs -x *.index \
		tmp/methgpc expected/methgpc]
	join [list_remove $result {}] \n
} {}

testsummarize


