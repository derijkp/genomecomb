#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl
set keepdir [pwd]

set test_cleantmp 0

# set optx {-x *.analysisinfo -x *flagstat*}
set ::optx {}

# tests
# =====

# gatk is not deterministic with threads, so no threads for testing and comparing to previous runs
lappend dopts -threads 1
set runopts {-stack 1}

# set dopts {--stack 1 --verbose 2 -threads 1 -d sge -map-mem 10G}

test process_small {process_project mastr_mx2} {
	cd $::smalltestdir
	set basename mastr_mx2
	file delete -force tmp/${basename}
	file mkdir tmp/${basename}/samples
	foreach sample [glob ori/${basename}.start/samples/*] {
		set sdest tmp/${basename}/samples/[string range [file tail $sample] 0 end-3]/fastq
		file mkdir $sdest
		foreach file [glob -nocomplain $sample/ori/*] {
			file copy $file $sdest
		}
	}
	file copy -force ori/wgs2.mastr/samplicons-wgs2.tsv tmp/${basename}/samplicons-wgs2.tsv
	# file copy ori/${basename}/demultiplex_stats.tsv tmp/${basename}
	# if you want to see output while running
	exec cg process_project {*}$::runopts {*}$::dopts -split 1 \
		-varcallers {gatk sam} \
		-reports -predictgender \
		-minfastqreads 10 \
		-amplicons tmp/${basename}/samplicons-wgs2.tsv \
		-extra_reports_mastr 1 \
		tmp/${basename} $::refseqdir/hg19 >& tmp/${basename}.log
	# check vs expected
	# filter out blanco hsmetrics results (are not deterministic)
	cg select -q {$sample ne "crsbwa-blanco2_8485"} tmp/${basename}/reports/report_hsmetrics-${basename}.tsv tmp/${basename}/reports/filtered_report_hsmetrics-${basename}.tsv
	cg select -q {!($sample eq "crsbwa-blanco2_8485" && $source eq "hsmetrics")} tmp/${basename}/reports/report_stats-${basename}.tsv tmp/${basename}/reports/filtered_report_stats-${basename}.tsv
	set result {}
	lappend result [tsvdiff -q 1 \
		-x *.bam -x *.bai -x *.tbi -x *.png -x *.zsti -x *.lz4i -x *.index \
		-x *.finished -x *.submitting -x *log_jobs -x info_analysis.tsv -x projectinfo.tsv \
		-x fastqc_report.html -x ${basename}.html -x report_stats-${basename}.tsv \
		-x colinfo -x *.stats.zst \
		-x ${basename}_hsmetrics_report.tsv -x report_hsmetrics-${basename}.tsv -x hsmetrics-crsbwa-blanco2_8485.hsmetrics \
		-x *-blanco2_8485* \
		-x *.html \
		-ignorefields {
			clipping_cg_version sammerge_version bamclean_version clipamplicons_version
			varcaller_version removeduplicates_version samstats_version hsmetrics_version flagstat_version regextrac_samtools
			regextract_version varcaller_cg_version annotate_cg_version varcaller_version histodepth_version histo_version
			report_vars_version predictgender_version report_covered_version svmulticompar_version regextrac_samtools
		} \
		tmp/${basename} expected/${basename}]
	lappend result [diffhtmlreport tmp/${basename}/reports/report-$basename.html expected/${basename}/reports/report-$basename.html]
	lappend result [diffanalysisinfo tmp/${basename}/compar/annot_compar-${basename}.tsv.analysisinfo expected/${basename}/compar/annot_compar-${basename}.tsv.analysisinfo]
	lappend result [checkdiff -I HistogramID -I htmlwidget -I {^<!} -I {^<h2>20} -I {meta charset} -I {script.src *=} \
		tmp/${basename}/${basename}.html expected/${basename}/${basename}.html]
	join [list_remove $result {}] \n
} {}

test process_small {process_project mastr_mx2_gatkh} {
	cd $::smalltestdir
	set basename mastr_mx2_gatkh
	file delete -force tmp/${basename}
	file mkdir tmp/${basename}/samples
	foreach sample [glob ori/mastr_mx2.start/samples/*] {
		set sdest tmp/${basename}/samples/[string range [file tail $sample] 0 end-3]/fastq
		file mkdir $sdest
		foreach file [glob -nocomplain $sample/ori/*] {
			file copy $file $sdest
		}
	}
	file copy -force ori/wgs2.mastr/samplicons-wgs2.tsv tmp/${basename}/samplicons-wgs2.tsv
	# file copy ori/${basename}/demultiplex_stats.tsv tmp/${basename}
	# if you want to see output while running
	cg process_project {*}$::runopts {*}$::dopts -split 1 -reports -predictgender \
		-minfastqreads 10 -amplicons tmp/${basename}/samplicons-wgs2.tsv -extra_reports_mastr 1 \
		-varcallers {gatkh} -extra_reports_mastr gatkh \
		tmp/${basename} $::refseqdir/hg19 >& tmp/${basename}.log
	# check vs expected
	# filter out blanco hsmetrics results (are not deterministic)
	cg select -q {$sample ne "crsbwa-blanco2_8485"} tmp/${basename}/reports/report_hsmetrics-${basename}.tsv tmp/${basename}/reports/filtered_report_hsmetrics-${basename}.tsv
	cg select -q {!($sample eq "crsbwa-blanco2_8485" && $source eq "hsmetrics")} tmp/${basename}/reports/report_stats-${basename}.tsv tmp/${basename}/reports/filtered_report_stats-${basename}.tsv
	set result {}
	lappend result [tsvdiff -q 1 \
		-x *.bam -x *.bai -x *.tbi -x *.png -x *.zsti -x *.lz4i -x *.index \
		-x *.finished -x *.submitting -x *log_jobs -x info_analysis.tsv -x projectinfo.tsv \
		-x fastqc_report.html -x ${basename}.html -x report_stats-${basename}.tsv \
		-x colinfo -x *.stats.zst \
		-x ${basename}_hsmetrics_report.tsv -x report_hsmetrics-${basename}.tsv -x hsmetrics-crsbwa-blanco2_8485.hsmetrics \
		-x *-blanco2_8485* \
		{*}[get ::optx {}] \
		-x *.html \
		-ignorefields {
			clipping_cg_version sammerge_version bamclean_version clipamplicons_version alignedsamstats_version unalignedsamstats_version
			varcaller_version removeduplicates_version samstats_version hsmetrics_version flagstat_version regextrac_samtools
			regextract_version varcaller_cg_version annotate_cg_version histodepth_version histo_version
			report_vars_version predictgender_version report_covered_version svmulticompar_version
		} \
		tmp/${basename} expected/${basename}]
	lappend result [diffhtmlreport tmp/${basename}/reports/report-$basename.html expected/${basename}/reports/report-$basename.html]
	lappend result [diffanalysisinfo tmp/${basename}/compar/annot_compar-${basename}.tsv.analysisinfo expected/${basename}/compar/annot_compar-${basename}.tsv.analysisinfo]
	lappend result [checkdiff -I version_os -I param_dbfiles -I param_dbdir -I command -I version_genomecomb -I maxopenfiles -I maxopenfiles \
		tmp/${basename}/compar/info_analysis.tsv expected/${basename}/compar/info_analysis.tsv]
	lappend result [checkdiff -I HistogramID -I htmlwidget -I {^<!} -I {^<h2>20} -I {meta charset} -I {script.src *=} \
		tmp/${basename}/${basename}.html expected/${basename}/${basename}.html]
	join [list_remove $result {}] \n
} {}

test process_small {process_project mastr_mx2 cram gatkh and strelka} {
	cd $::smalltestdir
	set basename mastr_mx2_cram
	file delete -force tmp/${basename}
	file mkdir tmp/${basename}/samples
	foreach sample [glob ori/mastr_mx2.start/samples/*] {
		set sdest tmp/${basename}/samples/[string range [file tail $sample] 0 end-3]/fastq
		file mkdir $sdest
		foreach file [glob -nocomplain $sample/ori/*] {
			file copy $file $sdest
		}
	}
	file copy -force ori/wgs2.mastr/samplicons-wgs2.tsv tmp/${basename}/samplicons-wgs2.tsv
	# file copy ori/${basename}/demultiplex_stats.tsv tmp/${basename}
	# if you want to see output while running
	cg process_project {*}$::runopts {*}$::dopts -split 1 -reports -predictgender \
		-minfastqreads 10 -amplicons tmp/${basename}/samplicons-wgs2.tsv \
		-extra_reports_mastr 1 \
		-varcallers {gatkh strelka} \
		-aliformat cram \
		-extra_reports_mastr gatkh \
		tmp/${basename} $::refseqdir/hg19 >& tmp/${basename}.log
	# check vs expected
	# filter out blanco hsmetrics results (are not deterministic)
	cg select -q {$sample ne "crsbwa-blanco2_8485"} tmp/${basename}/reports/report_hsmetrics-${basename}.tsv tmp/${basename}/reports/filtered_report_hsmetrics-${basename}.tsv
	cg select -q {!($sample eq "crsbwa-blanco2_8485" && $source eq "hsmetrics")} tmp/${basename}/reports/report_stats-${basename}.tsv tmp/${basename}/reports/filtered_report_stats-${basename}.tsv
	set result {}
	lappend result [tsvdiff -q 1 -x *log_jobs -x *.index \
		-x *.zsti -x *.lz4i -x *.tbi -x *.png -x *.cram -x *.crai -x *.bai \
		-x *.finished -x *.submitting -x info_analysis.tsv \
		-x projectinfo.tsv -x *.vcf.gz\
		-x fastqc_report.html -x colinfo -x *.stats.zst -x ${basename}.html -x hsmetrics-crsbwa-blanco2_8485.hsmetrics \
		-x ${basename}_hsmetrics_report.tsv -x report_hsmetrics-${basename}.tsv -x report_stats-${basename}.tsv \
		{*}[get ::optx {}] \
		-x *.html \
		-ignorefields {
			clipping_cg_version sammerge_version bamclean_version clipamplicons_version
			varcaller_version removeduplicates_version samstats_version hsmetrics_version flagstat_version regextrac_samtools
			regextract_version varcaller_cg_version annotate_cg_version histodepth_version histo_version
			report_vars_version predictgender_version report_covered_version svmulticompar_version
		} \
		tmp/${basename} expected/${basename}]
	lappend result [diffhtmlreport tmp/${basename}/reports/report-$basename.html expected/${basename}/reports/report-$basename.html]
	lappend result [diffanalysisinfo tmp/${basename}/compar/annot_compar-${basename}.tsv.analysisinfo expected/${basename}/compar/annot_compar-${basename}.tsv.analysisinfo]
	lappend result [checkdiff -I HistogramID -I htmlwidget -I {^<!} -I {^<h2>20} -I {meta charset} -I {script.src *=} \
		tmp/${basename}/${basename}.html expected/${basename}/${basename}.html]
	lappend result [checkdiff -I version_os -I param_dbfiles -I command -I param_dbdir -I version_genomecomb \
		tmp/${basename}/compar/info_analysis.tsv expected/${basename}/compar/info_analysis.tsv]
	join [list_remove $result {}] \n
} {}

#test process_small {process_project mastr_mx2 space in name} {
#	cd $::smalltestdir
#	set mastrdir {tmp/mastr_mx2 space}
#	set expected {expected/mastr_mx2 space}
#	file delete -force $mastrdir
#	file mkdir $mastrdir/samples
#	foreach sample [glob ori/mastr_mx2.start/samples/*] {
#		set sdest $mastrdir/samples/[string range [file tail $sample] 0 end-3]/fastq
#		file mkdir $sdest
#		foreach file [glob -nocomplain $sample/ori/*] {
#			file copy $file $sdest
#		}
#	}
#	file copy -force ori/wgs2.mastr/samplicons-wgs2.tsv $mastrdir/samplicons-wgs2.tsv
#	# file copy ori/mastr_mx2/demultiplex_stats.tsv $mastrdir
#	# if you want to see output while running
#	cg process_project {*}$::runopts {*}$::dopts -split 1 \
#		-minfastqreads 10 -amplicons $mastrdir/samplicons-wgs2.tsv -extra_reports_mastr 1 \
#		$mastrdir $::refseqdir/hg19 >& tmp/mastr_mx2_space.log
#	# check vs expected
#	set result {}
#	lappend result [tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x *.index -x fastqc_report.html \
#		-x colinfo -x *.stats.zst -x "mastr_mx2 space.html" -x *.zsti -x *.lz4i -x *.finished -x info_analysis.tsv -x *.png -x *.submitting \
#		-x *.analysisinfo \
#		$mastrdir $expected
#	lappend result [diffanalysisinfo $mastrdir/compar/annot_compar-*.tsv.analysisinfo $expected/compar/annot_compar-*.tsv.analysisinfo]
#	lappend result [checkdiff -y --suppress-common-lines "$mastrdir/mastr_mx2 space.html" "$expected/mastr_mx2 space.html" | grep -v -E {HistogramID|htmlwidget-|^<!|^<h2>20|meta charset|script.src *=}]
#	lappend result [checkdiff -y --suppress-common-lines $mastrdir/compar/info_analysis.tsv $expected/compar/info_analysis.tsv | grep -v -E {version_os|param_adapterfile|param_targetvarsfile|param_dbfiles|command|version_genomecomb}]
#	join [list_remove $result {}] \n
#} {}

test process_small {process_project -jobsample 1 mastr_mx2_js1} {
	cd $::smalltestdir
	set basename mastr_mx2_js1
	file delete -force tmp/${basename}
	file mkdir tmp/${basename}/samples
	foreach sample [glob ori/mastr_mx2.start/samples/*] {
		set sdest tmp/${basename}/samples/[string range [file tail $sample] 0 end-3]/fastq
		file mkdir $sdest
		foreach file [glob -nocomplain $sample/ori/*] {
			file copy $file $sdest
		}
	}
	file copy -force ori/wgs2.mastr/samplicons-wgs2.tsv tmp/${basename}/samplicons-wgs2.tsv
	cg process_project {*}$::runopts {*}$::dopts -jobsample 1 -split 1 \
		-varcallers {gatk sam} \
		-minfastqreads 10 \
		-amplicons tmp/${basename}/samplicons-wgs2.tsv \
		-extra_reports_mastr 1 \
		tmp/${basename} $::refseqdir/hg19 >& tmp/${basename}.log
	# check vs expected
	# filter out blanco hsmetrics results (are not deterministic)
	cg select -q {$sample ne "crsbwa-blanco2_8485"} tmp/${basename}/reports/report_hsmetrics-${basename}.tsv tmp/${basename}/reports/filtered_report_hsmetrics-${basename}.tsv
	cg select -q {!($sample eq "crsbwa-blanco2_8485" && $source eq "hsmetrics")} tmp/${basename}/reports/report_stats-${basename}.tsv tmp/${basename}/reports/filtered_report_stats-${basename}.tsv
	set result {}
	lappend result [tsvdiff -q 1 \
		-x *.bam -x *.bai -x *.tbi -x *.png -x *.zsti -x *.lz4i -x *.index \
		-x *.finished -x *.submitting -x *log_jobs -x info_analysis.tsv -x projectinfo.tsv \
		-x fastqc_report.html -x ${basename}.html -x report_stats-${basename}.tsv \
		-x colinfo -x *.stats.zst \
		-x ${basename}_hsmetrics_report.tsv -x report_hsmetrics-${basename}.tsv -x hsmetrics-crsbwa-blanco2_8485.hsmetrics \
		-x *-blanco2_8485* \
		{*}[get ::optx {}] \
		-x *.html \
		-ignorefields {
			clipping_cg_version sammerge_version bamclean_version clipamplicons_version
			varcaller_version removeduplicates_version samstats_version hsmetrics_version flagstat_version regextrac_samtools
			regextract_version varcaller_cg_version annotate_cg_version histodepth_version histo_version
			report_vars_version predictgender_version report_covered_version svmulticompar_version
		} \
		tmp/${basename} expected/${basename}]
	lappend result [diffhtmlreport tmp/${basename}/reports/report-$basename.html expected/${basename}/reports/report-$basename.html]
	lappend result [diffanalysisinfo tmp/${basename}/compar/annot_compar-*.tsv.analysisinfo expected/${basename}/compar/annot_compar-*.tsv.analysisinfo]
	lappend result [checkdiff -I version_os -I param_dbfiles -I param_dbdir -I command -I version_genomecomb -I maxopenfiles \
		tmp/${basename}/compar/info_analysis.tsv expected/${basename}/compar/info_analysis.tsv]
	lappend result [checkdiff -I HistogramID -I htmlwidget -I {^<!} -I {^<h2>20} -I {meta charset} -I {script.src *=} \
		tmp/${basename}/${basename}.html expected/${basename}/${basename}.html]
	join [list_remove $result {}] \n
} {}

test process_small {process_sample one_exome_yri_mx2} {
	cd $::smalltestdir
	set basename one_exome_yri_mx2
	file delete -force tmp/${basename}
	file mkdir tmp/${basename}/samples
	foreach sample {
		 NA19240mx2
	} {
		file mkdir tmp/${basename}/samples/$sample/fastq
		file copy {*}[glob ori/exomes_yri_mx2.start/samples/$sample/ori/*.fq.gz] tmp/${basename}/samples/$sample/fastq
		mklink $::refseqdir/hg19/extra/reg_hg19_exome_SeqCap_EZ_v3.tsv.lz4 tmp/${basename}/samples/$sample/reg_hg19_targets.tsv.lz4
	}
	cg process_sample {*}$::runopts {*}$::dopts -split 1 \
		-dbdir $::refseqdir/hg19 \
		tmp/${basename}/samples/NA19240mx2 >& tmp/${basename}.log
	# cg process_sample --stack 1 --verbose 2 -d status -split 1 -dbdir $::refseqdir/hg19 tmp/${basename}/samples/NA19240mx2 | less -S
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 \
		-x *.bam -x *.bai -x *.tbi -x *.png -x *.zsti -x *.lz4i -x *.index \
		-x *.finished -x *.submitting -x *log_jobs -x info_analysis.tsv -x projectinfo.tsv \
		-x fastqc_report.html -x ${basename}.html -x report_stats-${basename}.tsv \
		-x colinfo -x *.stats.zst \
		-x ${basename}_hsmetrics_report.tsv -x report_hsmetrics-${basename}.tsv -x hsmetrics-crsbwa-blanco2_8485.hsmetrics \
		-x *bam.dupmetrics \
		{*}[get ::optx {}] \
		-ignorefields {
			clipping_cg_version sammerge_version bamclean_version clipamplicons_version
			varcaller_version removeduplicates_version samstats_version hsmetrics_version flagstat_version regextrac_samtools
			regextract_version varcaller_cg_version annotate_cg_version histodepth_version histo_version
			report_vars_version predictgender_version report_covered_version svmulticompar_version
		} \
		tmp/${basename}/samples/NA19240mx2 expected/${basename}/samples/NA19240mx2]
	lappend result [checkdiff -I version_os -I param_dbfiles -I param_dbdir -I command -I version_genomecomb -I maxopenfiles \
		tmp/${basename}/samples/NA19240mx2/info_analysis.tsv expected/${basename}/samples/NA19240mx2/info_analysis.tsv]
	join [list_remove $result {}] \n
} {}

test process_small {process_sample one_d_exome_yri_mx2 distrreg} {
	cd $::smalltestdir
	set basename one_d_exome_yri_mx2
	file delete -force tmp/${basename}
	file mkdir tmp/${basename}/samples
	foreach sample {
		 NA19240mx2
	} {
		file mkdir tmp/${basename}/samples/$sample/fastq
		file copy {*}[glob ori/exomes_yri_mx2.start/samples/$sample/ori/*.fq.gz] tmp/${basename}/samples/$sample/fastq
		mklink $::refseqdir/hg19/extra/reg_hg19_exome_SeqCap_EZ_v3.tsv.lz4 tmp/${basename}/samples/$sample/reg_hg19_targets.tsv.lz4
	}
	cg process_sample {*}$::runopts {*}$::dopts -split 1 \
		-dbdir $::refseqdir/hg19 -distrreg 1 \
		tmp/${basename}/samples/NA19240mx2 >& tmp/${basename}.log
	# cg process_sample --stack 1 --verbose 2 -d status -split 1 -dbdir $::refseqdir/hg19 tmp/${basename}/samples/NA19240mx2 | less -S
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 \
		-x *.bam -x *.bai -x *.tbi -x *.png -x *.zsti -x *.lz4i -x *.index \
		-x *.finished -x *.submitting -x *log_jobs -x info_analysis.tsv -x projectinfo.tsv \
		-x fastqc_report.html -x ${basename}.html -x report_stats-${basename}.tsv \
		-x colinfo -x *.stats.zst \
		-x ${basename}_hsmetrics_report.tsv -x report_hsmetrics-${basename}.tsv -x hsmetrics-crsbwa-blanco2_8485.hsmetrics \
		-ignorefields {
			clipping_cg_version sammerge_version bamclean_version clipamplicons_version
			varcaller_version removeduplicates_version samstats_version hsmetrics_version flagstat_version regextrac_samtools
			regextract_version varcaller_cg_version annotate_cg_version histodepth_version histo_version
			report_vars_version predictgender_version report_covered_version svmulticompar_version
		} \
		{*}[get ::optx {}] \
		tmp/${basename}/samples/NA19240mx2 expected/${basename}/samples/NA19240mx2]
	lappend result [checkdiff -I version_os -I param_dbfiles -I param_dbdir -I command -I version_genomecomb -I maxopenfiles \
		tmp/${basename}/samples/NA19240mx2/info_analysis.tsv expected/${basename}/samples/NA19240mx2/info_analysis.tsv]
	join [list_remove $result {}] \n
} {}

test process_small {process_project exomes_yri_mx2} {
	cd $::smalltestdir
	set basename exomes_yri_mx2
	file delete -force tmp/${basename}
	file mkdir tmp/${basename}/samples
	foreach sample {
		NA19238mx2  NA19239mx2  NA19240mx2
	} {
		file mkdir tmp/${basename}/samples/$sample/fastq
		file copy {*}[glob ori/${basename}.start/samples/$sample/ori/*.fq.gz] tmp/${basename}/samples/$sample/fastq
		mklink $::refseqdir/hg19/extra/reg_hg19_exome_SeqCap_EZ_v3.tsv.lz4 tmp/${basename}/samples/$sample/reg_hg19_targets.tsv.lz4
	}
	cg process_project {*}$::runopts {*}$::dopts -split 1 \
		-varcallers {gatk sam} \
		-dbdir $::refseqdir/hg19 tmp/${basename} >& tmp/${basename}.log
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 \
		-x *.bam -x *.bai -x *.tbi -x *.png -x *.zsti -x *.lz4i -x *.index \
		-x *.finished -x *.submitting -x *log_jobs -x info_analysis.tsv -x projectinfo.tsv \
		-x fastqc_report.html -x ${basename}.html -x report_stats-${basename}.tsv \
		-x colinfo -x *.stats.zst \
		-x ${basename}_hsmetrics_report.tsv -x report_hsmetrics-${basename}.tsv -x hsmetrics-crsbwa-blanco2_8485.hsmetrics \
		{*}[get ::optx {}] \
		-x *.html \
		-ignorefields {
			clipping_cg_version sammerge_version bamclean_version clipamplicons_version
			varcaller_version removeduplicates_version samstats_version hsmetrics_version flagstat_version regextrac_samtools
			regextract_version varcaller_cg_version annotate_cg_version histodepth_version histo_version
			report_vars_version predictgender_version report_covered_version svmulticompar_version
		} \
		tmp/${basename} expected/${basename}]
	lappend result [diffhtmlreport tmp/${basename}/reports/report-$basename.html expected/${basename}/reports/report-$basename.html]
	lappend result [diffanalysisinfo tmp/${basename}/compar/annot_compar-*.tsv.analysisinfo expected/${basename}/compar/annot_compar-*.tsv.analysisinfo]
	# lappend result [checkdiff -y --suppress-common-lines tmp/${basename}/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics expected/${basename}/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics | grep -v "Started on" | grep -v bammarkduplicates2]
	foreach file1 [glob tmp/${basename}/compar/info_analysis.tsv tmp/${basename}/samples/*/info_analysis.tsv] {
		regsub ^tmp $file1 expected file2
		lappend result [checkdiff -I version_os -I param_dbfiles -I param_dbdir -I command -I version_genomecomb -I maxopenfiles \
			-I param_adapterfile \
			$file1 $file2]
	}
	join [list_remove $result {}] \n
} {}

test process_small {process_project exomesfb_yri_mx2 (freebayes)} {
	cd $::smalltestdir
	set basename exomesfb_yri_mx2
	file delete -force tmp/${basename}
	file mkdir tmp/${basename}/samples
	foreach sample {
		NA19238mx2  NA19239mx2  NA19240mx2
	} {
		file mkdir tmp/${basename}/samples/$sample/fastq
		file copy {*}[glob ori/exomes_yri_mx2.start/samples/$sample/ori/*.fq.gz] tmp/${basename}/samples/$sample/fastq
		mklink $::refseqdir/hg19/extra/reg_hg19_exome_SeqCap_EZ_v3.tsv.lz4 tmp/${basename}/samples/$sample/reg_hg19_targets.tsv.lz4
	}
	cg process_project {*}$::runopts {*}$::dopts -split 1 \
		-varcallers {gatk freebayes} \
		-dbdir $::refseqdir/hg19 tmp/${basename} >& tmp/${basename}.log
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 \
		-x *.bam -x *.bai -x *.tbi -x *.png -x *.zsti -x *.lz4i -x *.index \
		-x *.finished -x *.submitting -x *log_jobs -x info_analysis.tsv -x projectinfo.tsv \
		-x fastqc_report.html -x ${basename}.html -x report_stats-${basename}.tsv \
		-x colinfo -x *.stats.zst \
		-x ${basename}_hsmetrics_report.tsv -x report_hsmetrics-${basename}.tsv -x hsmetrics-crsbwa-blanco2_8485.hsmetrics \
		{*}[get ::optx {}] \
		-x *.html \
		-ignorefields {
			clipping_cg_version sammerge_version bamclean_version clipamplicons_version
			varcaller_version removeduplicates_version samstats_version hsmetrics_version flagstat_version regextrac_samtools
			regextract_version varcaller_cg_version annotate_cg_version histodepth_version histo_version
			report_vars_version predictgender_version report_covered_version svmulticompar_version
		} \
		tmp/${basename} expected/${basename}]
	lappend result [diffhtmlreport tmp/${basename}/reports/report-$basename.html expected/${basename}/reports/report-$basename.html]
	lappend result [diffanalysisinfo tmp/${basename}/compar/annot_compar-*.tsv.analysisinfo expected/${basename}/compar/annot_compar-*.tsv.analysisinfo]
	# lappend result [checkdiff -y --suppress-common-lines tmp/${basename}/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics expected/${basename}/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics | grep -v "Started on" | grep -v bammarkduplicates2]
	foreach file1 [glob tmp/${basename}/compar/info_analysis.tsv tmp/${basename}/samples/*/info_analysis.tsv] {
		regsub ^tmp $file1 expected file2
		lappend result [checkdiff -I version_os -I param_dbfiles -I param_dbdir -I command -I version_genomecomb -I maxopenfiles \
			-I param_adapterfile \
			$file1 $file2]
	}
	join [list_remove $result {}] \n
} {}

test process_small {process_project exomes_gatkh_yri_mx2 (haplotypecaller)} {
	cd $::smalltestdir
	set basename exomes_gatkh_yri_mx2
	file delete -force tmp/${basename}
	file mkdir tmp/${basename}/samples
	foreach sample {
		NA19238mx2  NA19239mx2  NA19240mx2
	} {
		file mkdir tmp/${basename}/samples/$sample/fastq
		file copy {*}[glob ori/exomes_yri_mx2.start/samples/$sample/ori/*.fq.gz] tmp/${basename}/samples/$sample/fastq
		mklink $::refseqdir/hg19/extra/reg_hg19_exome_SeqCap_EZ_v3.tsv.lz4 tmp/${basename}/samples/$sample/reg_hg19_targets.tsv.lz4
	}
	cg process_project {*}$::runopts {*}$::dopts -split 1 \
		-varcallers {gatkh freebayes} \
		-dbdir $::refseqdir/hg19 tmp/${basename} >& tmp/${basename}.log
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 \
		-x *.bam -x *.bai -x *.tbi -x *.png -x *.zsti -x *.lz4i -x *.index \
		-x *log_jobs -x *.finished -x *.submitting \
		-x *bam.dupmetrics -x info_analysis.tsv -x projectinfo.tsv -x *.analysisinfo \
		-x colinfo -x *.stats.zst -x fastqc_report.html -x *.vcf.gz \
		{*}[get ::optx {}] \
		-x *.html \
		-ignorefields {
			clipping_cg_version sammerge_version bamclean_version clipamplicons_version
			varcaller_version removeduplicates_version samstats_version hsmetrics_version flagstat_version regextrac_samtools
			regextract_version varcaller_cg_version annotate_cg_version histodepth_version histo_version
			report_vars_version predictgender_version report_covered_version svmulticompar_version
		} \
		tmp/${basename} expected/${basename}]
	lappend result [diffhtmlreport tmp/${basename}/reports/report-$basename.html expected/${basename}/reports/report-$basename.html]
	lappend result [diffanalysisinfo tmp/${basename}/compar/annot_compar-*.tsv.analysisinfo expected/${basename}/compar/annot_compar-*.tsv.analysisinfo]
	# lappend result [checkdiff -y --suppress-common-lines tmp/${basename}/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics expected/${basename}/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics | grep -v "Started on" | grep -v bammarkduplicates2]
	foreach file1 [glob tmp/${basename}/compar/info_analysis.tsv tmp/${basename}/samples/*/info_analysis.tsv] {
		regsub ^tmp $file1 expected file2
		lappend result [checkdiff -I version_os -I param_dbfiles -I param_dbdir -I command -I version_genomecomb -I maxopenfiles \
			-I param_adapterfile \
			$file1 $file2]
	}
	join [list_remove $result {}] \n
} {}

#test process_small {process_project exomes yri mx2} {
#	cd $::smalltestdir
#	file delete -force tmp/exomes_yri_mx2
#	file mkdir tmp/exomes_yri_mx2/samples
#	foreach sample {
#		NA19238mx2  NA19239mx2  NA19240mx2
#	} {
#		file mkdir tmp/exomes_yri_mx2/samples/$sample/fastq
#		file copy {*}[glob ori/exomes_yri_mx2.start/samples/$sample/ori/*.fq.gz] tmp/exomes_yri_mx2/samples/$sample/fastq
#	}
#	cg process_project --stack 1 --verbose 2 {*}$::runopts {*}$::dopts -split 1 -dbdir $::refseqdir/hg19 tmp/exomes_yri_mx2 >& tmp/exomes_yri_mx2.log
#	# check vs expected
#	checkdiff -y --suppress-common-lines tmp/exomes_yri_mx2/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics expected/exomes_yri_mx2/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics | grep -v "Started on" | grep -v bammarkduplicates2
#	checkdiff -qr -x *log_jobs -x *.bam -x *.bai -x colinfo -x *.stats.zst -x *_fastqc -x *bam.dupmetrics tmp/exomes_yri_mx2 expected/exomes_yri_mx2
#	# could have used this, but previous is faster
#	# cg tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x colinfo -x *.stats.zst -x *_fastqc -x *bam.dupmetrics tmp/exomes_yri_mx2 expected/exomes_yri_mx2 > temp
#} {}

test process_small {process_sample one_genome_yri_mx2} {
	cd $::smalltestdir
	set basename one_genome_yri_mx2
	set ref $::refseqdir/hg19
	file delete -force tmp/${basename}
	file mkdir tmp/${basename}/samples/NA19240cgmx2
	mklink ori/genomes_yri_mx2.start/samples/NA19240cgmx2/ori tmp/${basename}/samples/NA19240cgmx2/ori
	cg process_sample {*}$::runopts {*}$::dopts -split 1 \
		-dbdir $ref tmp/${basename}/samples/NA19240cgmx2 >& tmp/${basename}.log
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 \
		-x *.bam -x *.bai -x *.tbi -x *.png -x *.zsti -x *.lz4i -x *.index \
		-x *.finished -x *.submitting -x *log_jobs -x info_analysis.tsv -x projectinfo.tsv \
		-x summary-NA19240cgmx2.txt \
		{*}[get ::optx {}] \
		-ignorefields {
			clipping_cg_version sammerge_version bamclean_version clipamplicons_version
			varcaller_version removeduplicates_version samstats_version hsmetrics_version flagstat_version regextrac_samtools
			regextract_version varcaller_cg_version annotate_cg_version histodepth_version histo_version
			report_vars_version predictgender_version report_covered_version svmulticompar_version
		} \
		tmp/${basename}/samples/NA19240cgmx2 expected/genomes_yri_mx2/samples/NA19240cgmx2]
	join [list_remove $result {}] \n
} {}

test process_small {process_project genomes_yri_mx2} {
	cd $::smalltestdir
	set basename genomes_yri_mx2
	set dest tmp/${basename}
	file delete -force tmp/${basename}
	file mkdir tmp/${basename}
	foreach sample {
		NA19238cgmx2 NA19239cgmx2 NA19240cgmx2
		NA19240ilmx2
	} {
		file mkdir tmp/${basename}/samples/$sample
		mklink ori/${basename}.start/samples/$sample/ori tmp/${basename}/samples/$sample/ori
	}
	# cg process_project --stack 1 --verbose 2 -d 2 -split 1 -dbdir /complgen/refseq/testdb2/hg19 tmp/${basename}
	cg process_project {*}$::runopts {*}$::dopts -split 1 \
		-varcallers {gatk sam} \
		-dbdir $::refseqdir/hg19 tmp/${basename} >& tmp/${basename}.log
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 \
		-x *.bam -x *.bai -x *.tbi -x *.png -x *.zsti -x *.lz4i -x *.index \
		-x *.finished -x *.submitting -x *log_jobs -x info_analysis.tsv -x projectinfo.tsv \
		-x fastqc_report.html -x ${basename}.html -x report_stats-${basename}.tsv \
		-x colinfo -x *.stats.zst \
		-x ${basename}_hsmetrics_report.tsv -x report_hsmetrics-${basename}.tsv -x hsmetrics-crsbwa-blanco2_8485.hsmetrics \
		-x summary-*.txt \
		{*}[get ::optx {}] \
		-x *.html \
		-ignorefields {
			clipping_cg_version sammerge_version bamclean_version clipamplicons_version
			varcaller_version removeduplicates_version samstats_version hsmetrics_version flagstat_version regextrac_samtools
			regextract_version varcaller_cg_version annotate_cg_version histodepth_version histo_version
			report_vars_version predictgender_version report_covered_version svmulticompar_version
		} \
		tmp/${basename} expected/${basename}]
	lappend result [diffhtmlreport tmp/${basename}/reports/report-$basename.html expected/${basename}/reports/report-$basename.html]
	lappend result [diffanalysisinfo tmp/${basename}/compar/annot_compar-*.tsv.analysisinfo expected/${basename}/compar/annot_compar-*.tsv.analysisinfo]
	foreach cgsample {NA19238cgmx2 NA19239cgmx2 NA19240cgmx2} {
		lappend result [checkdiff -I finished \
			tmp/${basename}/samples/$cgsample/summary-$cgsample.txt expected/${basename}/samples/$cgsample/summary-$cgsample.txt]
	}
	foreach file1 [glob tmp/${basename}/compar/info_analysis.tsv tmp/${basename}/samples/*/info_analysis.tsv] {
		regsub ^tmp $file1 expected file2
		lappend result [checkdiff -I version_os -I param_dbfiles -I param_dbdir -I command -I version_genomecomb -I maxopenfiles \
			-I param_adapterfile \
			$file1 $file2]
	}
	join [list_remove $result {}] \n
} {}

test process_small {process_project cg_mx2} {
	cd $::smalltestdir
	set basename cg_mx2
	set dest tmp/${basename}
	file delete -force tmp/${basename}
	file mkdir tmp/${basename}
	cg project_addsample -transfer rel tmp/${basename} cgNA19240mx2 ori/mixed_yri_mx2/cgNA19240mx2
	cg process_project {*}$::runopts {*}$::dopts -split 1 \
		-varcallers {gatk sam} \
		-dbdir $::refseqdir/hg19 \
		-dbfile [gzfile /complgen/refseq/hg19/extra/var_hg19_dbnsfp.tsv] \
		tmp/${basename} >& tmp/${basename}.log
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 \
		-x *.bam -x *.bai -x *.tbi -x *.png -x *.zsti -x *.lz4i -x *.index \
		-x *.finished -x *.submitting -x *log_jobs -x info_analysis.tsv -x projectinfo.tsv \
		-x fastqc_report.html -x ${basename}.html -x report_stats-${basename}.tsv \
		-x colinfo -x *.stats.zst -x summary-*.txt \
		-x ${basename}_hsmetrics_report.tsv -x report_hsmetrics-${basename}.tsv -x hsmetrics-crsbwa-blanco2_8485.hsmetrics \
		-x *dupmetrics \
		{*}[get ::optx {}] \
		-x *.html \
		-ignorefields {
			clipping_cg_version sammerge_version bamclean_version clipamplicons_version
			varcaller_version removeduplicates_version samstats_version hsmetrics_version flagstat_version regextrac_samtools
			regextract_version varcaller_cg_version annotate_cg_version histodepth_version histo_version
			report_vars_version predictgender_version report_covered_version svmulticompar_version
		} \
		tmp/${basename} expected/${basename}]
	foreach cgsample {cgNA19240mx2} {
		lappend result [checkdiff -I finished \
			tmp/$basename/samples/$cgsample/summary-$cgsample.txt expected/$basename/samples/$cgsample/summary-$cgsample.txt]
	}
	lappend result [diffanalysisinfo tmp/${basename}/compar/annot_compar-*.tsv.analysisinfo expected/${basename}/compar/annot_compar-*.tsv.analysisinfo]
	# lappend result [checkdiff -y --suppress-common-lines tmp/${basename}/samples/gilNA19240mx2/map-dsbwa-gilNA19240mx2.bam.dupmetrics expected/${basename}/samples/gilNA19240mx2/map-dsbwa-gilNA19240mx2.bam.dupmetrics | grep -v "Started on" | grep -v bammarkduplicates2]
	foreach file1 [glob tmp/$basename/compar/info_analysis.tsv tmp/$basename/samples/*/info_analysis.tsv] {
		regsub ^tmp $file1 expected file2
		lappend result [checkdiff -I version_os -I param_dbfiles -I param_dbdir -I command -I version_genomecomb -I maxopenfiles \
			-I param_adapterfile \
			$file1 $file2]
	}
	join [list_remove $result {}] \n
} {}

test process_small {process_project mixed_yri_mx2} {
	cd $::smalltestdir
	set basename mixed_yri_mx2
	set dest tmp/${basename}
	file delete -force tmp/${basename}
	file mkdir tmp/${basename}
	cg project_addsample -transfer rel tmp/${basename} cgNA19240mx2 ori/${basename}/cgNA19240mx2
	cg project_addsample -transfer rel tmp/${basename} gilNA19240mx2 {*}[glob ori/${basename}/gilNA19240mx2/*.fq.gz]
	cg project_addsample -transfer rel -targetfile ori/${basename}/reg_hg19_exome_SeqCap_EZ_v3.tsv.lz4 tmp/${basename} exNA19239mx2 {*}[glob ori/${basename}/exNA19239mx2/*.fq.gz]
	cg project_addsample -transfer rel -targetfile ori/${basename}/reg_hg19_exome_SeqCap_EZ_v3.tsv.lz4 tmp/${basename} exNA19240mx2 ori/${basename}/exNA19240mx2
	cg process_project {*}$::runopts {*}$::dopts -split 1 \
		-varcallers {gatk sam} \
		-dbdir $::refseqdir/hg19 \
		-dbfile [gzfile /complgen/refseq/hg19/extra/var_hg19_dbnsfp.tsv] \
		tmp/${basename} >& tmp/${basename}.log
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 \
		-x *.bam -x *.bai -x *.tbi -x *.png -x *.zsti -x *.lz4i -x *.index \
		-x *.finished -x *.submitting -x *log_jobs -x info_analysis.tsv -x projectinfo.tsv \
		-x fastqc_report.html -x ${basename}.html -x report_stats-${basename}.tsv \
		-x colinfo -x *.stats.zst \
		-x ${basename}_hsmetrics_report.tsv -x report_hsmetrics-${basename}.tsv -x hsmetrics-crsbwa-blanco2_8485.hsmetrics \
		-x *dupmetrics \
		{*}[get ::optx {}] \
		-x *.html \
		-ignorefields {
			clipping_cg_version sammerge_version bamclean_version clipamplicons_version
			varcaller_version removeduplicates_version samstats_version hsmetrics_version flagstat_version regextrac_samtools
			regextract_version varcaller_cg_version annotate_cg_version histodepth_version histo_version
			report_vars_version predictgender_version report_covered_version svmulticompar_version
		} \
		tmp/${basename} expected/${basename}]
	lappend result [diffhtmlreport tmp/${basename}/reports/report-$basename.html expected/${basename}/reports/report-$basename.html]
	lappend result [diffanalysisinfo tmp/${basename}/compar/annot_compar-*.tsv.analysisinfo expected/${basename}/compar/annot_compar-*.tsv.analysisinfo]
	# lappend result [checkdiff -y --suppress-common-lines tmp/${basename}/samples/gilNA19240mx2/map-dsbwa-gilNA19240mx2.bam.dupmetrics expected/${basename}/samples/gilNA19240mx2/map-dsbwa-gilNA19240mx2.bam.dupmetrics | grep -v "Started on" | grep -v bammarkduplicates2]
	foreach file1 [glob tmp/$basename/compar/info_analysis.tsv tmp/$basename/samples/*/info_analysis.tsv] {
		regsub ^tmp $file1 expected file2
		lappend result [checkdiff -I version_os -I param_dbfiles -I param_dbdir -I command -I version_genomecomb -I maxopenfiles \
			-I param_adapterfile \
			$file1 $file2]
	}
	join [list_remove $result {}] \n
} {}

test process_small {process_project -distrreg 1 mixed_yri_mx2_distrreg} {
	cd $::smalltestdir
	set basename mixed_yri_mx2_distrreg
	set dest tmp/${basename}
	file delete -force tmp/${basename}
	file mkdir tmp/${basename}
	cg project_addsample -transfer rel tmp/${basename} cgNA19240mx2 ori/mixed_yri_mx2/cgNA19240mx2
	cg project_addsample -transfer rel tmp/${basename} gilNA19240mx2 {*}[glob ori/mixed_yri_mx2/gilNA19240mx2/*.fq.gz]
	cg project_addsample -transfer rel -targetfile ori/mixed_yri_mx2/reg_hg19_exome_SeqCap_EZ_v3.tsv.lz4 tmp/${basename} exNA19239mx2 {*}[glob ori/mixed_yri_mx2/exNA19239mx2/*.fq.gz]
	cg project_addsample -transfer rel -targetfile ori/mixed_yri_mx2/reg_hg19_exome_SeqCap_EZ_v3.tsv.lz4 tmp/${basename} exNA19240mx2 ori/mixed_yri_mx2/exNA19240mx2
	cg process_project {*}$::runopts {*}$::dopts -distrreg 1 -split 1 \
		-varcallers {gatk sam} \
		-dbdir $::refseqdir/hg19 \
		-dbfile [gzfile /complgen/refseq/hg19/extra/var_hg19_dbnsfp.tsv] \
		tmp/${basename} >& tmp/${basename}.log
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 \
		-x *.bam -x *.bai -x *.tbi -x *.png -x *.zsti -x *.lz4i -x *.index \
		-x *.finished -x *.submitting -x *log_jobs -x info_analysis.tsv -x projectinfo.tsv \
		-x fastqc_report.html -x ${basename}.html -x report_stats-${basename}.tsv \
		-x colinfo -x *.stats.zst -x summary-*.txt \
		-x ${basename}_hsmetrics_report.tsv -x report_hsmetrics-${basename}.tsv -x hsmetrics-crsbwa-blanco2_8485.hsmetrics \
		{*}[get ::optx {}] \
		-x *.html \
		-ignorefields {
			clipping_cg_version sammerge_version bamclean_version clipamplicons_version
			varcaller_version removeduplicates_version samstats_version hsmetrics_version flagstat_version regextrac_samtools
			regextract_version varcaller_cg_version annotate_cg_version histodepth_version histo_version
			report_vars_version predictgender_version report_covered_version svmulticompar_version
		} \
		tmp/${basename} expected/${basename}]
	lappend result [diffhtmlreport tmp/${basename}/reports/report-$basename.html expected/${basename}/reports/report-$basename.html]
	lappend result [diffanalysisinfo tmp/${basename}/compar/annot_compar-*.tsv.analysisinfo expected/${basename}/compar/annot_compar-*.tsv.analysisinfo]
	# lappend result [checkdiff -y --suppress-common-lines tmp/${basename}/samples/gilNA19240mx2/map-dsbwa-gilNA19240mx2.bam.dupmetrics expected/${basename}/samples/gilNA19240mx2/map-dsbwa-gilNA19240mx2.bam.dupmetrics | grep -v "Started on" | grep -v bammarkduplicates2]
	foreach file1 [glob tmp/$basename/compar/info_analysis.tsv tmp/$basename/samples/*/info_analysis.tsv] {
		regsub ^tmp $file1 expected file2
		lappend result [checkdiff -I version_os -I param_dbfiles -I param_dbdir -I command -I version_genomecomb -I maxopenfiles \
			-I param_adapterfile \
			$file1 $file2]
	}
	join [list_remove $result {}] \n
} {}

if 0 {

test process_small {annotate refseqbuild/hg19} {
	cd $::smalltestdir
	file delete tmp/annot_exomes_yri_mx2.tsv
	cg annotate --stack 1 --verbose 2 expected/exomes_yri_mx2/compar/compar-exomes_yri_mx2.tsv tmp/annot_exomes_yri_mx2.tsv /data/genomecomb.testdata/refseqbuild/hg19 /data/genomecomb.testdata/refseqbuild/hg19/extra >& tmp/annot_exomes_yri_mx2.log
	checkdiff tmp/annot_exomes_yri_mx2.tsv expected/annot_exomes_yri_mx2.tsv
} {}

test process_small {annotate refseqbuild/hg38} {
	cd $::smalltestdir
	# compar-exomes_yri_mx2.tsv is actually hg19, but we use it here anyway just to test if the hg38 works
	file delete tmp/annothg38_exomes_yri_mx2.tsv
	cg annotate --stack 1 --verbose 2 expected/exomes_yri_mx2/compar/compar-exomes_yri_mx2.tsv tmp/annothg38_exomes_yri_mx2.tsv /data/genomecomb.testdata/refseqbuild/hg38 /data/genomecomb.testdata/refseqbuild/hg38/extra >& tmp/annothg38_exomes_yri_mx2.log
	chekdiff tmp/annothg38_exomes_yri_mx2.tsv expected/annothg38_exomes_yri_mx2.tsv
}

}

testsummarize


