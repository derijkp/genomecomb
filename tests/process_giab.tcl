#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl
set keepdir [pwd]

if {![info exists argv]} {set argv {}}
set download 0
set test_cleantmp 0
# check or run
set run 0
set forcebenchmark 0
if {[inlist $argv run]} {
	set run 1
}
if {[inlist $argv forcebenchmark]} {
	set forcebenchmark 1
}

proc benchmarkvars {args} {
	set target [lindex $args end]
	if {!$::forcebenchmark && [file exists $target]} {
		puts "skipping $target: already exists and no forcebenchmark"
		return $target
	}
	exec cg benchmarkvars {*}$args
}

# Download database
# =================
if {[get download 0] != 0} {
	# mkdir a directory ~/public/giab and download publically available giab data to it
	if {[isont $download] || $download in "sge slurm"} {set d $download} else {set d 1}
	cg_giab_getdata -d $d precisionfda_v2016_04 ~/public/giab/precisionfda_v2016_04
	cg_giab_gettruth -ref hg38 3.3.2 ~/public/giab/truth/truth_hg38_v3.3.2
	cg_giab_gettruth -ref hg38 4.2.1 ~/public/giab/truth/truth_hg38_v4.2.1

	cg_giab_getdata -d $d platinum_genomes ~/public/platinum_genomes
	cg_giab_gettruth -ref hg38 hybrid ~/public/platinum_genomes/truthset/2017-1.0/hg38/hybrid
}

# tests
# =====

# tests will run in ~/genomecomb_giab_testdata
# this directy should be present, and have a dir expected with the expected results to compare to
# It also expects the giab data to be present in ~/public/giab (as downloaded with previous code)

test process_giab {precisionFDA} {
	mkdir ~/genomecomb_giab_testdata
	cd ~/genomecomb_giab_testdata
	if {$::run} {
		file delete -force tmp/precisionFDA
		file mkdir tmp/precisionFDA/samples
		foreach sample {
			HG001_NA12878 HG002_NA24385_son
		} {
			file mkdir tmp/precisionFDA/samples/$sample/fastq
			foreach file [glob public/giab/fastqs/precisionfda_v2016_04/$sample/split/*fastq.gz] {
				mklink $file tmp/precisionFDA/samples/$sample/fastq/[file tail $file]
			}
		}
		foreach {sample version} {
			HG001_NA12878 3.3.2
			HG002_NA24385_son 3.3.2
			HG002_NA24385_son 4.2.1
		} {
			if {$version eq "3.3.2"} {
				set samplename truth_$sample
			} else {
				set samplename truth[string_change $version {. {}}]_$sample
			}
			set sampledir tmp/precisionFDA/samples/$samplename
			mkdir $sampledir
			set vcf [glob public/giab/truth/truth_hg38_v$version/*/${sample}_hg38/*.vcf.gz]
			set region [glob public/giab/truth/truth_hg38_v$version/*/${sample}_hg38/*.bed.tsv.zst]
			mklink $vcf $sampledir/var-truth-truth-$samplename.vcf.gz 1
			mklink [file root [gzroot $vcf]].tsv.zst $sampledir/var-truth-truth-$samplename.tsv.zst 1
			mklink $region $sampledir/sreg-truth-truth-$samplename.tsv.zst 1
		}
	
		# run
		exec devcg process_project -stack 1 -v 2 -d sge -split 1 \
			-threads 8 -distrreg 30000000 -varcallers {gatkh strelka} -svcallers {manta lumpy} \
			-reports {all} \
			-dbdir /complgen/refseq/hg38 \
			tmp/precisionFDA >& tmp/precisionFDA.log
		puts "precisionFDA run started"
	} else {
		# check vs expected
		if {$forcebenchmark} {
			file delete -force tmp/precisionFDA/benchmarks
		}
		mkdir tmp/precisionFDA/benchmarks
		benchmarkvars -analyses {strelka-rdsbwa-HG001_NA12878 gatkh-rdsbwa-HG001_NA12878} \
			tmp/precisionFDA/compar/annot_compar-precisionFDA.tsv.zst \
			truth-truth-truth_HG001_NA12878 \
			tmp/precisionFDA/benchmarks/benchmark-HG001_NA12878.tsv
		benchmarkvars -analyses {strelka-rdsbwa-HG002_NA24385_son gatkh-rdsbwa-HG002_NA24385_son} \
			tmp/precisionFDA/compar/annot_compar-precisionFDA.tsv.zst \
			truth-truth-truth_HG002_NA24385_son \
			tmp/precisionFDA/benchmarks/benchmark-HG002_NA24385_son.tsv
		benchmarkvars -analyses {strelka-rdsbwa-HG001_NA12878 gatkh-rdsbwa-HG001_NA12878} \
			-regionfile tmp/precisionFDA/samples/truth_HG001_NA12878/sreg-truth-truth-truth_HG001_NA12878.tsv.zst \
			tmp/precisionFDA/compar/annot_compar-precisionFDA.tsv.zst \
			truth-truth-truth_HG001_NA12878 \
			tmp/precisionFDA/benchmarks/benchmark-reg-HG001_NA12878.tsv
		benchmarkvars -analyses {strelka-rdsbwa-HG002_NA24385_son gatkh-rdsbwa-HG002_NA24385_son} \
			-regionfile tmp/precisionFDA/samples/truth_HG002_NA24385_son/sreg-truth-truth-truth_HG002_NA24385_son.tsv.zst \
			tmp/precisionFDA/compar/annot_compar-precisionFDA.tsv.zst \
			truth-truth-truth_HG002_NA24385_son \
			tmp/precisionFDA/benchmarks/benchmark-reg-HG002_NA24385_son.tsv
		benchmarkvars -analyses {strelka-rdsbwa-HG002_NA24385_son gatkh-rdsbwa-HG002_NA24385_son} \
			-regionfile tmp/precisionFDA/samples/truth421_HG002_NA24385_son/sreg-truth-truth-truth421_HG002_NA24385_son.tsv.zst \
			tmp/precisionFDA/compar/annot_compar-precisionFDA.tsv.zst \
			truth-truth-truth421_HG002_NA24385_son \
			tmp/precisionFDA/benchmarks/benchmark-reg421-HG002_NA24385_son.tsv
		if {$::forcebenchmark || ![file exists tmp/precisionFDA/benchmarks/hap.py-gatkh-rdsbwa-HG002_NA24385_son.summary.csv]]} {
			exec devcg giab_benchmark -d sge \
				-refseq /complgen/refseq/hg38 \
				tmp/precisionFDA/hap.py- \
				tmp/precisionFDA/samples/truth_HG002_NA24385_son \
				tmp/precisionFDA/samples/HG002_NA24385_son/var-gatkh-rdsbwa-HG002_NA24385_son.vcf.gz \
				tmp/precisionFDA/samples/HG002_NA24385_son/var-strelka-rdsbwa-HG002_NA24385_son.vcf.gz
		}
		#
		set result {}
		lappend result [tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x colinfo -x fastqc_report.html \
			-x *bam.dupmetrics -x info_analysis.tsv -x *.zsti -x *.lz4i -x *.finished -x *.index -x info_analysis.tsv \
			-x *.analysisinfo -x *.png -x *.submitting \
			tmp/precisionFDA expected/precisionFDA]
		lappend result [diffanalysisinfo tmp/precisionFDA/compar/annot_compar-*.tsv.analysisinfo expected/precisionFDA/compar/annot_compar-*.tsv.analysisinfo]
		# lappend result [checkdiff -y --suppress-common-lines tmp/precisionFDA/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics expected/precisionFDA/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics | grep -v "Started on" | grep -v bammarkduplicates2]
		foreach file1 [glob tmp/precisionFDA/compar/info_analysis.tsv tmp/precisionFDA/samples/*/info_analysis.tsv] {
			regsub ^tmp $file1 expected file2
			lappend result [checkdiff -y --suppress-common-lines $file1 $file2 | grep -v -E {version_os|param_adapterfile|param_targetvarsfile|param_dbfiles|command|version_genomecomb}]
		}
		join [list_remove $result {}] \n
	}
} {}

if 0 {
	# try bqsr
	exec devcg gatk_bqsr -d sge -dbdir /complgen/refseq/hg38 tmp/precisionFDA/samples/HG002_NA24385_son/map-rdsbwa-HG002_NA24385_son.bam tmp/precisionFDA/samples/HG002_NA24385_son/map-bqsr-rdsbwa-HG002_NA24385_son.bam
	cg var -d sge -method gatkh -stack 1 -v 2 -distrreg 1 -datatype exome -mincoverage 5 -mingenoqual 12 \
		tmp/precisionFDA/samples/HG002_NA24385_son/map-bqsr-rdsbwa-HG002_NA24385_son.bam \
		/complgen/refseq/hg38/genome_hg38.ifas > tmp/log 2> tmp/logerror
	cg var -d sge -method strelka -stack 1 -v 2 -distrreg 1 -datatype exome -mincoverage 5 -mingenoqual 12 \
		tmp/precisionFDA/samples/HG002_NA24385_son/map-bqsr-rdsbwa-HG002_NA24385_son.bam \
		/complgen/refseq/hg38/genome_hg38.ifas > tmp/log 2> tmp/logerror
	devcg giab_benchmark -d sge \
		-refseq /complgen/refseq/hg38 \
		tmp/precisionFDA/hap.py- \
		tmp/precisionFDA/samples/truth_HG002_NA24385_son \
		tmp/precisionFDA/samples/HG002_NA24385_son/var-gatkh-bqsr-rdsbwa-HG002_NA24385_son.vcf.gz \
		tmp/precisionFDA/samples/HG002_NA24385_son/var-gatkh-rdsbwa-HG002_NA24385_son.vcf.gz \
		tmp/precisionFDA/samples/HG002_NA24385_son/var-strelka-bqsr-rdsbwa-HG002_NA24385_son.vcf.gz \
		tmp/precisionFDA/samples/HG002_NA24385_son/var-strelka-rdsbwa-HG002_NA24385_son.vcf.gz \
}

proc addtruth {dir truthsets} {
	foreach {sample version} $truthsets {
		if {$version eq "3.3.2"} {
			set samplename truth_$sample
		} else {
			set samplename truth[string_change $version {. {}}]_$sample
		}
		set sampledir $dir/samples/$samplename
		mkdir $sampledir
		set vcf [glob public/giab/truth/truth_hg38_v$version/*/${sample}_hg38/*.vcf.gz]
		set region [glob public/giab/truth/truth_hg38_v$version/*/${sample}_hg38/*.bed.tsv.zst]
		mklink $vcf $sampledir/var-truth-truth-$samplename.vcf.gz 1
		mklink [file root [gzroot $vcf]].tsv.zst $sampledir/var-truth-truth-$samplename.tsv.zst 1
		mklink $region $sampledir/sreg-truth-truth-$samplename.tsv.zst 1
	}
	
}

test process_giab {giab_one} {
	cd ~/genomecomb_giab_testdata
	if {$::run} {
		file delete -force tmp/giab_one
		file mkdir tmp/giab_one/samples
		foreach sample {
			NA12878
		} {
			file mkdir tmp/giab_one/samples/$sample/fastq
			foreach file [glob public/platinum_genomes/*_$sample/fastqsplit/*fastq.gz] {
				mklink $file tmp/giab_one/samples/$sample/fastq/[file tail $file]
			}
		}
		addtruth tmp/giab_one {
			HG001_NA12878 3.3.2
		}
		# run
		exec devcg process_project -stack 1 -v 2 -d sge -split 1 \
			-threads 8 -distrreg 30000000 -varcallers {gatkh strelka} -svcallers {manta lumpy} \
			-reports {all} \
			-dbdir /complgen/refseq/hg38 \
			tmp/giab_one >& tmp/giab_one.log
	} else {
		# check vs expected
		benchmarkvars tmp/giab_one/compar/annot_compar-giab_one.tsv.zst \
			truth-truth-truth_HG001_NA12878 \
			tmp/giab_one/benchmark.tsv
		#
		set result {}
		lappend result [tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x colinfo -x fastqc_report.html \
			-x *bam.dupmetrics -x info_analysis.tsv -x *.zsti -x *.lz4i -x *.finished -x *.index -x info_analysis.tsv \
			-x *.analysisinfo -x *.png -x *.submitting \
			tmp/giab_one expected/giab_one]
		lappend result [diffanalysisinfo tmp/giab_one/compar/annot_compar-*.tsv.analysisinfo expected/giab_one/compar/annot_compar-*.tsv.analysisinfo]
		# lappend result [checkdiff -y --suppress-common-lines tmp/giab_one/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics expected/giab_one/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics | grep -v "Started on" | grep -v bammarkduplicates2]
		foreach file1 [glob tmp/giab_one/compar/info_analysis.tsv tmp/giab_one/samples/*/info_analysis.tsv] {
			regsub ^tmp $file1 expected file2
			lappend result [checkdiff -y --suppress-common-lines $file1 $file2 | grep -v -E {version_os|param_adapterfile|param_targetvarsfile|param_dbfiles|command|version_genomecomb}]
		}
		join [list_remove $result {}] \n
	}
} {}

test process_giab {process_giab} {
	cd ~/genomecomb_giab_testdata
	if {$::run} {
		file delete -force tmp/giab
		file mkdir tmp/giab/samples
		foreach sample {
			NA12878	NA12891	NA12892
		} {
			file mkdir tmp/giab/samples/$sample/fastq
			foreach file [glob public/platinum_genomes/*_$sample/fastqsplit/*fastq.gz] {
				mklink $file tmp/giab/samples/$sample/fastq/[file tail $file]
			}
		}
		addtruth tmp/giab {
			HG001_NA12878 3.3.2
		}
#		mkdir tmp/giab/samples/hybridref_NA12878
#		mklink public/platinum_genomes/truthset/2017-1.0/hg38/hybrid/var_hg38.hybrid.tsv.zst tmp/giab/samples/hybridref_NA12878/var-ref-ref-hybridref_NA12878.tsv.zst
#		mklink public/platinum_genomes/truthset/2017-1.0/hg38/hybrid/reg_hg38.hybrid.tsv.zst tmp/giab/samples/hybridref_NA12878/sreg-ref-ref-hybridref_NA12878.tsv.zst
		# run
		exec devcg process_project -stack 1 -v 2 -d sge -split 1 \
			-threads 8 -distrreg 30000000 -varcallers {gatkh strelka} -svcallers {manta lumpy} \
			-reports {all} \
			-dbdir /complgen/refseq/hg38 \
			tmp/giab >& tmp/giab.log
	} else {
		# check vs expected
		benchmarkvars \
			-analyses {strelka-rdsbwa-NA12878 gatkh-rdsbwa-NA12878} \
			-regionfile tmp/giab/samples/truth_HG001_NA12878/sreg-truth-truth-truth_HG001_NA12878.tsv.zst \
			tmp/giab/compar/annot_compar-giab.tsv.zst \
			truth-truth-truth_HG001_NA12878 \
			tmp/giab/benchmark.tsv
		set result {}
		lappend result [tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x colinfo -x fastqc_report.html \
			-x *bam.dupmetrics -x info_analysis.tsv -x *.zsti -x *.lz4i -x *.finished -x *.index -x info_analysis.tsv \
			-x *.analysisinfo -x *.png -x *.submitting \
			tmp/giab expected/giab]
		lappend result [diffanalysisinfo tmp/giab/compar/annot_compar-*.tsv.analysisinfo expected/giab/compar/annot_compar-*.tsv.analysisinfo]
		# lappend result [checkdiff -y --suppress-common-lines tmp/giab/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics expected/giab/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics | grep -v "Started on" | grep -v bammarkduplicates2]
		foreach file1 [glob tmp/giab/compar/info_analysis.tsv tmp/giab/samples/*/info_analysis.tsv] {
			regsub ^tmp $file1 expected file2
			lappend result [checkdiff -y --suppress-common-lines $file1 $file2 | grep -v -E {version_os|param_adapterfile|param_targetvarsfile|param_dbfiles|command|version_genomecomb}]
		}
		join [list_remove $result {}] \n
	}
} {}

if {$run} {
	puts "All tests started"
} else {
	testsummarize
}




