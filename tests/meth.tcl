#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl
set keepdir [pwd]

set test_cleantmp 0

# tests
# =====

# gatk is not deterministic with threads, so no threads for testing and comparing to previous runs
lappend dopts -threads 1
set runopts {-stack 1}

if 0 {

	# download testdata
	cd $::smalltestdir/ori
	wget https://f5c.page.link/f5c_na12878_test
	mv f5c_na12878_test f5c_na12878_test.tar.gz
	tar xvzf f5c_na12878_test.tar.gz
	mv chr22_meth_example ont_f5c_chr22_meth_example
	cd ont_f5c_chr22_meth_example
	mkdir fast5
	mv fast5_files single_fast5_files
	~/build/bin-x86_64/ont-fast5-api/bin/python3.8 ~/build/bin-x86_64/ont-fast5-api/bin/single_to_multi_fast5 -i single_fast5_files -s fast5 -n 4000
	~/build/bin-x86_64/ont-fast5-api/bin/python3.8 ~/build/bin-x86_64/ont-fast5-api/bin/single_to_multi_fast5 -i single_fast5_files -s fast5 -n 4000
	#
	unset -nocomplain fast5file2batcha
	unset -nocomplain readid2batcha
	set f [open fast5/filename_mapping.txt]
	while {[gets $f line] != -1} {
		foreach {file batch} [split $line \t] break
		set fast5file2batcha($file) $batch
	}
	close $f
	set f [open reads.fastq.index.readdb]
	while {[gets $f line] != -1} {
		foreach {readid file} [split $line \t] break
		if {$file eq ""} continue
		set file [file tail $file]
		set readid2batcha($readid) $fast5file2batcha($file)
	}
	close $f
	file mkdir fastq
	unset -nocomplain fa
	set f [open reads.fastq]
	while {[gets $f line] != -1} {
		set readid [string range [lindex $line 0] 1 end]
		putsvars readid
		if {![info exists readid2batcha($readid)]} {
			puts "skipping $readid: not found"
			gets $f
			gets $f
			gets $f
			continue
		}
		set batch $readid2batcha($readid)
		if {![info exists fa($batch)]} {
			set fastqfile [file root [file tail $batch]].fastq
			puts "Creating $fastqfile"
			set fa($batch) [open fastq/$fastqfile w]
		}
		puts $fa($batch) $line
		puts $fa($batch) [gets $f]
		puts $fa($batch) [gets $f]
		puts $fa($batch) [gets $f]
	}
	close $f
	foreach batch [array names fa] {
		close $fa($batch)
	}
	foreach file [glob fastq/*.fastq] {
		exec bgzip $file
	}
}

test meth {meth_nanopolish} {
	cd $::smalltestdir
	file delete -force tmp/meth
	file mkdir tmp/meth/samples/methtest/fast5
	file mkdir tmp/meth/samples/methtest/fastq
	foreach file [glob -nocomplain ori/ont_f5c_chr22_meth_example/fast5/*] {
		mklink $file tmp/meth/samples/methtest/fast5/[file tail $file]
	}
	foreach file [glob -nocomplain ori/ont_f5c_chr22_meth_example/fastq/*] {
		mklink $file tmp/meth/samples/methtest/fastq/[file tail $file]
	}
	file mkdir tmp/meth/samples/methtest2/fast5
	file mkdir tmp/meth/samples/methtest2/fastq
	foreach file {ori/ont_f5c_chr22_meth_example/fast5/batch_0.fast5 ori/ont_f5c_chr22_meth_example/fast5/batch_1.fast5} {
		mklink $file tmp/meth/samples/methtest2/fast5/[file tail $file]
	}
	foreach file {ori/ont_f5c_chr22_meth_example/fastq/batch_0.fastq.gz ori/ont_f5c_chr22_meth_example/fastq/batch_1.fastq.gz} {
		mklink $file tmp/meth/samples/methtest2/fastq/[file tail $file]
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
			-svcallers {sniffles npinv} \
			-varcallers longshot \
			-methcallers nanopolish \
			-hap_bam 1 \
			-threads 6 \
			-reports {-fastqc predictgender} \
			tmp/meth $::refseqdir/hg19 >& tmp/meth.log
	}]
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 -x *.bam -x *.bai -x fastqc_report.html \
		-x colinfo -x meth.html -x *.zsti -x *.lz4i -x *.finished -x info_analysis.tsv \
		-x *.analysisinfo -x *.png -x *.submitting \
		-x *log_jobs -x *.index \
		tmp/meth expected/meth]
} {}


testsummarize


