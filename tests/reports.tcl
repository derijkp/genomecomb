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
	cd $::bigtestdir
	file delete -force tmp/test_reports
	file mkdir tmp/test_reports/NA19240mx2
	file copy {*}[glob expected/exomes_yri_mx2/samples/NA19240mx2/*] tmp/test_reports/NA19240mx2
	file delete -force tmp/test_reports/NA19240mx2/reports
	mklink refseqtest/hg19/extra/reg_hg19_exome_SeqCap_EZ_v3.tsv.lz4 tmp/test_reports/NA19240mx2/reg_hg19_targets.tsv.lz4
	cg process_reports -stack 1 -v 0 tmp/test_reports/NA19240mx2 refseqtest/hg19 2>@ stderr >@ stdout
	file_regsub > >\n tmp/test_reportsnotarget/NA19240mx2/reports/fastqc_fw-NA19240mx2.fastqc/fastqc_report.html tmp/test_reportsnotarget/NA19240mx2/reports/fastqc_fw-NA19240mx2.fastqc/fastqc_report.html.temp
	file rename -force tmp/test_reportsnotarget/NA19240mx2/reports/fastqc_fw-NA19240mx2.fastqc/fastqc_report.html.temp tmp/test_reportsnotarget/NA19240mx2/reports/fastqc_fw-NA19240mx2.fastqc/fastqc_report.html
	file_regsub > >\n tmp/test_reportsnotarget/NA19240mx2/reports/fastqc_rev-NA19240mx2.fastqc/fastqc_report.html tmp/test_reportsnotarget/NA19240mx2/reports/fastqc_rev-NA19240mx2.fastqc/fastqc_report.html.temp
	file rename -force tmp/test_reportsnotarget/NA19240mx2/reports/fastqc_rev-NA19240mx2.fastqc/fastqc_report.html.temp tmp/test_reportsnotarget/NA19240mx2/reports/fastqc_rev-NA19240mx2.fastqc/fastqc_report.html
	cg tsvdiff -q 1 -x fastqc_report.html -x *.png tmp/test_reports/NA19240mx2/reports expected/exomes_yri_mx2/samples/NA19240mx2/reports
	catch {checkdiff -y --suppress-common-lines tmp/test_reports/NA19240mx2/reports/fastqc_fw-NA19240mx2.fastqc/fastqc_report.html expected/exomes_yri_mx2/samples/NA19240mx2/reports/fastqc_fw-NA19240mx2.fastqc/fastqc_report.html} e
	if {[llength [split [string trim $e] \n]] > 2} {error "too many differences in fastqc_report.html"}
	catch {checkdiff -y --suppress-common-lines tmp/test_reports/NA19240mx2/reports/fastqc_rev-NA19240mx2.fastqc/fastqc_report.html expected/exomes_yri_mx2/samples/NA19240mx2/reports/fastqc_rev-NA19240mx2.fastqc/fastqc_report.html} e
	if {[llength [split [string trim $e] \n]] > 2} {error "too many differences in fastqc_report.html"}
} {}

test reports {process_reports no targetfile} {
	cd $::bigtestdir
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
	catch {checkdiff -y --suppress-common-lines tmp/test_reportsnotarget/NA19240mx2/reports/fastqc_fw-NA19240mx2.fastqc/fastqc_report.html expected/test_reportsnotarget/NA19240mx2/reports/fastqc_fw-NA19240mx2.fastqc/fastqc_report.html} e
	if {[llength [split [string trim $e] \n]] > 2} {error "too many differences in fastqc_report.html"}
	catch {checkdiff -y --suppress-common-lines tmp/test_reportsnotarget/NA19240mx2/reports/fastqc_rev-NA19240mx2.fastqc/fastqc_report.html expected/test_reportsnotarget/NA19240mx2/reports/fastqc_rev-NA19240mx2.fastqc/fastqc_report.html} e
	if {[llength [split [string trim $e] \n]] > 2} {error "too many differences in fastqc_report.html"}
} {}

testsummarize
