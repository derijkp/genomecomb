#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test reports {hsmetrics} {
	test_cleantmp
	cg select -q {$chromosome in "chr21 chr22"} $::refseqdir/hg19/extra/reg_hg19_exome_SeqCap_EZ_v3.tsv.lz4 tmp/regfile.tsv
	set bamfile genomecomb.testdata/ori/test-map-rdsbwa-NA19240chr2122.bam
	set regionfile tmp/regfile.tsv
	set resultfile tmp/result.hsmetrics
	cg hsmetrics $bamfile $regionfile $resultfile 2> /dev/null
	cg select -rc 1 tmp/result.hsmetrics tmp/result.hsmetrics.nocomments
	exec diff tmp/result.hsmetrics.nocomments genomecomb.testdata/expected/bam_histo-NA19240chr2122.hsmetrics
} {}

test reports {coverage_report} {
	test_cleantmp
	cg select -f {chromosome=chr_clip($chromosome) begin end info} data/reg_hg19_smallpartexome.tsv tmp/regfile.tsv
	set regionfile tmp/regfile.tsv
	mklink genomecomb.testdata/ori/test-map-rdsbwa-NA19240chr2122.bam tmp/test.bam
	mklink genomecomb.testdata/ori/test-map-rdsbwa-NA19240chr2122.bam.bai tmp/test.bam.bai
	set bamfile tmp/test.bam
	set intervals {1 5 10 20 50 100 200 500 1000}
	cg coverage_report $regionfile $bamfile
	exec diff tmp/test.histo genomecomb.testdata/expected/bam_histo-NA19240_smallpartchr2122.tsv
} {}

test reports {report_vars} {
	mklink data/annot_compar-exomes_yri_parts.tsv tmp/vars.tsv
	cg report_vars -stack 1 -v 2 -sample gatk-rdsbwa-NA19238chr2122 \
		-targetfile $::refseqdir/hg19/extra/reg_hg19_exome_SeqCap_EZ_v3.tsv.lz4 \
		-refcodingfile $::refseqdir/hg19/extra/reg_hg19_refcoding.tsv.lz4 \
		tmp/vars.tsv tmp/report_vars.tsv
	exec diff tmp/report_vars.tsv data/report_vars.tsv
} {}

test reports {process_reports} {
	cd $::smalltestdir
	file delete -force tmp/test_reports
	file mkdir tmp/test_reports/NA19240mx2
	file copy {*}[glob expected/exomes_yri_mx2/samples/NA19240mx2/*] tmp/test_reports/NA19240mx2
	file delete -force tmp/test_reports/NA19240mx2/reports
	mklink refseqtest/hg19/extra/reg_hg19_exome_SeqCap_EZ_v3.tsv.lz4 tmp/test_reports/NA19240mx2/reg_hg19_targets.tsv.lz4
	cg process_reports -stack 1 -v 0 tmp/test_reports/NA19240mx2 refseqtest/hg19 2>@ stderr >@ stdout
	cg tsvdiff -q 1 -x fastqc_report.html -x *.png tmp/test_reports/NA19240mx2/reports expected/exomes_yri_mx2/samples/NA19240mx2/reports
	set e [checkdiff -y --suppress-common-lines tmp/test_reports/NA19240mx2/reports/fastqc_fw-NA19240mx2.fastqc/fastqc_report.html expected/exomes_yri_mx2/samples/NA19240mx2/reports/fastqc_fw-NA19240mx2.fastqc/fastqc_report.html]
	if {[llength [split [string trim $e] \n]] > 2} {error "too many differences in fastqc_report.html"}
	set e [checkdiff -y --suppress-common-lines tmp/test_reports/NA19240mx2/reports/fastqc_rev-NA19240mx2.fastqc/fastqc_report.html expected/exomes_yri_mx2/samples/NA19240mx2/reports/fastqc_rev-NA19240mx2.fastqc/fastqc_report.html]
	if {[llength [split [string trim $e] \n]] > 2} {error "too many differences in fastqc_report.html"}
} {}

test reports {process_reports no targetfile} {
	cd $::smalltestdir
	file delete -force tmp/test_reportsnotarget
	file mkdir tmp/test_reportsnotarget/NA19240mx2
	file copy {*}[glob expected/exomes_yri_mx2/samples/NA19240mx2/*] tmp/test_reportsnotarget/NA19240mx2
	file delete -force tmp/test_reportsnotarget/NA19240mx2/reports tmp/test_reportsnotarget/NA19240mx2/reg_hg19_targets.tsv.lz4
	cg process_reports -stack 1 -v 0 tmp/test_reportsnotarget/NA19240mx2 refseqtest/hg19 2>@ stderr >@ stdout
	file_regsub > >\n tmp/test_reportsnotarget/NA19240mx2/reports/fastqc_fw-NA19240mx2.fastqc/fastqc_report.html tmp/test_reportsnotarget/NA19240mx2/reports/fastqc_fw-NA19240mx2.fastqc/fastqc_report.html.temp
	file rename -force tmp/test_reportsnotarget/NA19240mx2/reports/fastqc_fw-NA19240mx2.fastqc/fastqc_report.html.temp tmp/test_reportsnotarget/NA19240mx2/reports/fastqc_fw-NA19240mx2.fastqc/fastqc_report.html
	file_regsub > >\n tmp/test_reportsnotarget/NA19240mx2/reports/fastqc_rev-NA19240mx2.fastqc/fastqc_report.html tmp/test_reportsnotarget/NA19240mx2/reports/fastqc_rev-NA19240mx2.fastqc/fastqc_report.html.temp
	file rename -force tmp/test_reportsnotarget/NA19240mx2/reports/fastqc_rev-NA19240mx2.fastqc/fastqc_report.html.temp tmp/test_reportsnotarget/NA19240mx2/reports/fastqc_rev-NA19240mx2.fastqc/fastqc_report.html
	cg tsvdiff -q 1 -x fastqc_report.html tmp/test_reportsnotarget/NA19240mx2/reports expected/test_reportsnotarget/NA19240mx2/reports
	set e [checkdiff -y --suppress-common-lines tmp/test_reportsnotarget/NA19240mx2/reports/fastqc_fw-NA19240mx2.fastqc/fastqc_report.html expected/test_reportsnotarget/NA19240mx2/reports/fastqc_fw-NA19240mx2.fastqc/fastqc_report.html]
	if {[llength [split [string trim $e] \n]] > 2} {error "too many differences in fastqc_report.html"}
	set e [checkdiff -y --suppress-common-lines tmp/test_reportsnotarget/NA19240mx2/reports/fastqc_rev-NA19240mx2.fastqc/fastqc_report.html expected/test_reportsnotarget/NA19240mx2/reports/fastqc_rev-NA19240mx2.fastqc/fastqc_report.html]
	if {[llength [split [string trim $e] \n]] > 2} {error "too many differences in fastqc_report.html"}
} {}

test reports {process_reportscombine} {
	mkdir tmp/samples/ERR194147_30x_NA12878
	file copy data/reports tmp/samples/ERR194147_30x_NA12878
	mkdir tmp/samples/test/reports
	foreach file {
		report_fastq_fw-ERR194147_30x_NA12878.tsv
		report_fastq_rev-ERR194147_30x_NA12878.tsv		
	} {
		regsub ERR194147_30x_NA12878 $file test newfile
		set c [file_read data/reports/$file]
		regsub -all ERR194147_30x_NA12878 $c test c
		file_write tmp/samples/test/reports/$newfile $c
	}
	cg select -overwrite 1 -f {depth ontarget {offtarget=int(0.9*$offtarget)}} data/reports/histodepth-rdsbwa-ERR194147_30x_NA12878.tsv tmp/samples/test/reports/histodepth-rdsbwa-test.tsv
	cg process_reportscombine {*}$::dopts -overwrite 1 -dbdir $::refseqdir/hg19 tmp/combinereports tmp/samples/ERR194147_30x_NA12878/reports tmp/samples/test/reports
	if {[catch {exec diff -r tmp/combinereports data/expected-combinereports} msg]} {
		if {![string match [deindent {
			diff -r tmp/combinereports/report-combinereports.html data/expected-combinereports/report-combinereports.html
			305c305
			< The full "sequenced genome" region in */reg_hg19_sequencedgenome.tsv.lz4 was used as target region for calulations in the table
			---
			> The full "sequenced genome" region in */reg_hg19_sequencedgenome.tsv.lz4 was used as target region for calulations in the table
			child process exited abnormally
			}] $msg]} {
				error $msg
		}
	}
} {}

test reports {process_reportscombine 2} {
	cd $::smalltestdir
	file delete -force tmp/combinereports
	cg process_reportscombine {*}$::dopts tmp/combinereports {*}[glob expected/exomes_yri_mx2/samples/* tmp/genomes_yri_mx2/samples/NA19240ilmx2/reports] expected/exomes_yri_mx2/samples/NA19240mx2/reports
	cg tsvdiff tmp/combinereports expected/combinereports
} {}

testsummarize
