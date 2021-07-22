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
	exec devcg benchmarkvars {*}$args
}

# Download database
# =================
if {[get download 0]} {
	
	# NA12878: Utah woman, parents are NA12891 and NA12892, genetic disease (CYP2D6 mutation)
	# pilot of giab (https://jimb.stanford.edu/giab)
	# also in platinum genomes
	# - (https://emea.illumina.com/platinumgenomes.html): 
	# - 17 member CEPH pedigree 1463 fully sequenced
	# - Eberle, MA et al. (2017) A reference data set of 5.4 million phased human variants validated by genetic inheritance from sequencing a three-generation 17-member pedigree. Genome Research 27: 157-164. doi:10.1101/gr.210500.116
	mkdir a directory ~/public/giab and download publically available giab data to it
	cg_giab_getfastqs -d sge giab_ont_ultralong ~/public/giab/fastqs/giab_ont_ultralong

	# limited regions (chr1:2500000-3000000 chr2:2000000-2600000 chr6:32000000-33000000 chr10:1500000-2000000)
	set file $::env(HOME)/public/giab/fastqs/giab_ont_ultralong_part12610/HG002_NA24385_son/HG002_NA24385_son_part12610.fastq.gz
	mkdir [file dir $file]
	# extract part (after running test process_giab_ont)
	exec samtools view -b -1 -h tmp/giab_ont/samples/HG002_NA24385_son/map-sminimap2-HG002_NA24385_son.bam chr1:2500000-3000000 chr2:2000000-2600000 chr6:32000000-33000000 chr10:1500000-2000000 > temp.bam
	cg bam2fastq -threads 8 temp.bam $file
	mkdir [file dir $file]/split
	cg fastq_split -d sge -parts 2000 $file [file dir $file]/split/[file tail $file]
	#
	set file $::env(HOME)/public/giab/fastqs/giab_ont_ultralong_part12610/HG003_NA24149_father/HG003_NA24149_father_part12610.fastq.gz
	mkdir [file dir $file]
	exec samtools view -b -1 -h tmp/giab_ont/samples/HG003_NA24149_father/map-sminimap2-HG003_NA24149_father.bam chr1:2500000-3000000 chr2:2000000-2600000 chr6:32000000-33000000 chr10:1500000-2000000 > temp.bam
	cg bam2fastq -threads 8 temp.bam $file
	file delete temp.bam
	mkdir [file dir $file]/split
	cg fastq_split -d sge -parts 2000 $file [file dir $file]/split/[file tail $file]
}

# extra code
# ==========
# go and take from public_genomes.proj/procedure-rungenomes.txt

# tests
# =====

#test process_giab {giab_one_ont} {
#	cd ~/genomecomb_giab_testdata
#	file delete -force tmp/giab_one_ont
#	file mkdir tmp/giab_one_ont/samples
#	foreach sample {
#		NA12878
#	} {
#		file mkdir tmp/giab_one_ont/samples/$sample/fastq
#		foreach file [glob platinum_genomes/ori/*_$sample/fastqsplit/*fastq.gz] {
#			mklink $file tmp/giab_one_ont/samples/$sample/fastq/[file tail $file]
#		}
#	}
#	mkdir tmp/giab_one_ont/samples/truth_NA12878
#	mklink ori/giab/truthset/truth_hg38_giab3.3.2/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.tsv.zst \
#		tmp/giab_one_ont/samples/truth_NA12878/var-truth-truth-truth_NA12878.tsv.zst 1
#	mklink ori/giab/truthset/truth_hg38_giab3.3.2/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed.tsv.zst \
#		tmp/giab_one_ont/samples/truth_NA12878/sreg-truth-truth-truth_NA12878.tsv.zst 1
#	# run
#	exec devcg process_project -stack 1 -v 2 -d sge -split 1 \
#		-threads 8 -distrreg 30000000 -varcallers {gatkh strelka} -svcallers {manta lumpy} \
#		-reports {all} \
#		-dbdir /complgen/refseq/hg38 \
#		tmp/giab_one_ont >& tmp/giab_one_ont.log
#	# check vs expected
#	cg benchmarkvars tmp/giab_one_ont/compar/annot_compar-giab_one_ont.tsv.zst truth-truth-truth_NA12878 tmp/giab_one_ont/banchmark.tsv
#	#
#	set result {}
#	lappend result [tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x colinfo -x fastqc_report.html \
#		-x *bam.dupmetrics -x info_analysis.tsv -x *.zsti -x *.lz4i -x *.finished -x *.index -x info_analysis.tsv \
#		-x *.analysisinfo -x *.png -x *.submitting \
#		tmp/giab_one_ont expected/giab_one_ont]
#	lappend result [diffanalysisinfo tmp/giab_one_ont/compar/annot_compar-*.tsv.analysisinfo expected/giab_one_ont/compar/annot_compar-*.tsv.analysisinfo]
#	# lappend result [checkdiff -y --suppress-common-lines tmp/giab_one_ont/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics expected/giab_one_ont/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics | grep -v "Started on" | grep -v bammarkduplicates2]
#	foreach file1 [glob tmp/giab_one_ont/compar/info_analysis.tsv tmp/giab_one_ont/samples/*/info_analysis.tsv] {
#		regsub ^tmp $file1 expected file2
#		lappend result [checkdiff -y --suppress-common-lines $file1 $file2 | grep -v -E {version_os|param_adapterfile|param_targetvarsfile|param_dbfiles|command|version_genomecomb}]
#	}
#	join [list_remove $result {}] \n
#} {}

test process_giab {process small_giab_ont} {
	cd ~/genomecomb_giab_testdata
	if {$::run} {
		file delete -force tmp/small_giab_ont
		file mkdir tmp/small_giab_ont/samples
		foreach sample {
			HG002_NA24385_son	HG003_NA24149_father
		} {
			file mkdir tmp/small_giab_ont/samples/$sample/fastq
			foreach file [glob public/giab/fastqs/giab_ont_ultralong_part12610/$sample/split/*fastq.gz] {
				mklink $file tmp/small_giab_ont/samples/$sample/fastq/[file tail $file]
			}
			foreach version {3.3.2 4.2.1} {
				if {$version eq "3.3.2"} {
					set samplename truth_$sample
				} else {
					set samplename truth[string_change $version {. {}}]_$sample
				}
				set sampledir tmp/small_giab_ont/samples/$samplename
				mkdir $sampledir
				set vcf [glob public/giab/truth/truth_hg38_v$version/*/${sample}_hg38/*.vcf.gz]
				set region [glob public/giab/truth/truth_hg38_v$version/*/${sample}_hg38/*.bed.tsv.zst]
				mklink $vcf $sampledir/var-truth-truth-$samplename.vcf.gz 1
				mklink [file root [gzroot $vcf]].tsv.zst $sampledir/var-truth-truth-$samplename.tsv.zst 1
				mklink $region $sampledir/sreg-truth-truth-$samplename.tsv.zst 1
			}
		}
		# run
		set varcallers {longshot}
		set svcallers {sniffles cuteSV cuteSV_pacbio npinv}
		exec devcg process_project -stack 1 -v 2 -d sge -split 1 -threads 8 \
			-paired 0 -clip 0 \
			-maxfastqdistr 250 \
			-aligner {minimap2} \
			-removeduplicates 0 \
			-realign 0 \
			-distrreg chr \
			-sniffles-n -1 \
			-svcallers $svcallers \
			-cuteSV-threads 16 \
			-varcallers $varcallers \
			-reports {-fastqc predictgender} \
			-dbdir /complgen/refseq/hg38 \
			tmp/small_giab_ont >& tmp/small_giab_ont.log
	} else {
		# check vs expected
		benchmarkvars \
			-analyses {longshot-sminimap2-HG003_NA24149_father} \
			-regionfile tmp/small_giab_ont/samples/truth421_HG003_NA24149_father/sreg-truth-truth-truth421_HG003_NA24149_father.tsv.zst \
			tmp/small_giab_ont/compar/annot_compar-small_giab_ont.tsv.zst truth-truth-truth421_HG003_NA24149_father tmp/small_giab_ont/benchmark_HG003.tsv
		set result {}
		lappend result [tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x colinfo -x fastqc_report.html \
			-x *bam.dupmetrics -x info_analysis.tsv -x *.zsti -x *.lz4i -x *.finished -x *.index -x info_analysis.tsv \
			-x *.analysisinfo -x *.png -x *.submitting \
			tmp/small_giab_ont expected/giab]
		lappend result [diffanalysisinfo tmp/small_giab_ont/compar/annot_compar-*.tsv.analysisinfo expected/giab/compar/annot_compar-*.tsv.analysisinfo]
		# lappend result [checkdiff -y --suppress-common-lines tmp/small_giab_ont/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics expected/giab/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics | grep -v "Started on" | grep -v bammarkduplicates2]
		foreach file1 [glob tmp/small_giab_ont/compar/info_analysis.tsv tmp/small_giab_ont/samples/*/info_analysis.tsv] {
			regsub ^tmp $file1 expected file2
			lappend result [checkdiff -y --suppress-common-lines $file1 $file2 | grep -v -E {version_os|param_adapterfile|param_targetvarsfile|param_dbfiles|command|version_genomecomb}]
		}
		join [list_remove $result {}] \n
	}
} {}

# downsample 3 times for smaller testset (still has size similar to our normal runs) and thus faster testing
test process_giab {process_giab_ont_ds3} {
	cd ~/genomecomb_giab_testdata
	if {$::run} {
		file delete -force tmp/giab_ont_ds3
		file mkdir tmp/giab_ont_ds3/samples
		foreach sample {
			HG002_NA24385_son	HG003_NA24149_father	HG004_NA24143_mother
		} {
			file mkdir tmp/giab_ont_ds3/samples/$sample/fastq
			# skip 2 fastqs everytime
			foreach {file temp temp} [glob public/giab/fastqs/giab_ont_ultralong/$sample/split/*fastq.gz] {
				mklink $file tmp/giab_ont_ds3/samples/$sample/fastq/[file tail $file]
			}
			foreach version {3.3.2 4.2.1} {
				if {$version eq "3.3.2"} {
					set samplename truth_$sample
				} else {
					set samplename truth[string_change $version {. {}}]_$sample
				}
				set sampledir tmp/giab_ont_ds3/samples/$samplename
				mkdir $sampledir
				set vcf [glob public/giab/truth/truth_hg38_v$version/*/${sample}_hg38/*.vcf.gz]
				set region [glob public/giab/truth/truth_hg38_v$version/*/${sample}_hg38/*.bed.tsv.zst]
				mklink $vcf $sampledir/var-truth-truth-$samplename.vcf.gz 1
				mklink [file root [gzroot $vcf]].tsv.zst $sampledir/var-truth-truth-$samplename.tsv.zst 1
				mklink $region $sampledir/sreg-truth-truth-$samplename.tsv.zst 1
			}
		}
		# run
		puts "Starting process_giab_ont_ds3 run"
		exec devcg process_project -stack 1 -v 2 -d sge -split 1 -threads 6 \
			-paired 0 -clip 0 \
			-maxfastqdistr 250 \
			-aligner {minimap2} \
			-removeduplicates 0 \
			-realign 0 \
			-distrreg chr \
			-sniffles-n -1 \
			-svcallers {sniffles cuteSV cuteSV_pacbio npinv} \
			-varcallers {longshot} \
			-reports {-fastqc predictgender} \
			-dbdir /complgen/refseq/hg38 \
			tmp/giab_ont_ds3 >& tmp/giab_ont_ds3.log
	} else {
		# check vs expected
		foreach sample {HG002_NA24385_son HG003_NA24149_father HG004_NA24143_mother} {
			benchmarkvars \
				-analyses longshot-sminimap2-$sample \
				-regionfile tmp/giab_ont_ds3/samples/truth421_$sample/sreg-truth-truth-truth421_$sample.tsv.zst \
				tmp/giab_ont_ds3/compar/annot_compar-giab_ont_ds3.tsv.zst \
				truth-truth-truth421_$sample \
				tmp/giab_ont_ds3/benchmark_truth421_$sample.tsv
			benchmarkvars \
				-analyses longshot-sminimap2-$sample \
				-regionfile tmp/giab_ont_ds3/samples/truth_$sample/sreg-truth-truth-truth_$sample.tsv.zst \
				tmp/giab_ont_ds3/compar/annot_compar-giab_ont_ds3.tsv.zst \
				truth-truth-truth_$sample \
				tmp/giab_ont_ds3/benchmark_truth_$sample.tsv
			#set truthsample tmp/giab_ont_ds3/samples/truth421_HG003_NA24149_father
			#set sample tmp/giab_ont_ds3/samples/HG003_NA24149_father
			#cg giab_benchmark -refseq /complgen/refseq/hg38 tmp/giab_ont_ds3/giab_benchmark_${sample}_ $truthsample $sample
		}
		set result {}
		lappend result [tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x colinfo -x fastqc_report.html \
			-x *bam.dupmetrics -x info_analysis.tsv -x *.zsti -x *.lz4i -x *.finished -x *.index -x info_analysis.tsv \
			-x *.analysisinfo -x *.png -x *.submitting \
			tmp/giab_ont_ds3 expected/giab]
		lappend result [diffanalysisinfo tmp/giab_ont_ds3/compar/annot_compar-*.tsv.analysisinfo expected/giab/compar/annot_compar-*.tsv.analysisinfo]
		# lappend result [checkdiff -y --suppress-common-lines tmp/giab_ont_ds3/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics expected/giab/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics | grep -v "Started on" | grep -v bammarkduplicates2]
		foreach file1 [glob tmp/giab_ont_ds3/compar/info_analysis.tsv tmp/giab_ont_ds3/samples/*/info_analysis.tsv] {
			regsub ^tmp $file1 expected file2
			lappend result [checkdiff -y --suppress-common-lines $file1 $file2 | grep -v -E {version_os|param_adapterfile|param_targetvarsfile|param_dbfiles|command|version_genomecomb}]
		}
		join [list_remove $result {}] \n
	}
} {}

# full data is too slow for easy testing
#test process_giab {process_giab_ont} {
#	cd ~/genomecomb_giab_testdata
#	if {$::run} {
#		file delete -force tmp/giab_ont
#		file mkdir tmp/giab_ont/samples
#		foreach sample {
#			HG002_NA24385_son	HG003_NA24149_father	HG004_NA24143_mother
#		} {
#			file mkdir tmp/giab_ont/samples/$sample/fastq
#			foreach file [glob public/giab/fastqs/giab_ont_ultralong/$sample/split/*fastq.gz] {
#				mklink $file tmp/giab_ont/samples/$sample/fastq/[file tail $file]
#			}
#			foreach version {3.3.2 4.2.1} {
#				if {$version eq "3.3.2"} {
#					set samplename truth_$sample
#				} else {
#					set samplename truth[string_change $version {. {}}]_$sample
#				}
#				set sampledir tmp/giab_ont/samples/$samplename
#				mkdir $sampledir
#				set vcf [glob public/giab/truth/truth_hg38_v$version/*/${sample}_hg38/*.vcf.gz]
#				set region [glob public/giab/truth/truth_hg38_v$version/*/${sample}_hg38/*.bed.tsv.zst]
#				mklink $vcf $sampledir/var-truth-truth-$samplename.vcf.gz 1
#				mklink [file root [gzroot $vcf]].tsv.zst $sampledir/var-truth-truth-$samplename.tsv.zst 1
#				mklink $region $sampledir/sreg-truth-truth-$samplename.tsv.zst 1
#			}
#		}
#		# run
#		puts "Starting process_giab_ont run"
#		exec devcg process_project -stack 1 -v 2 -d sge -split 1 -threads 6 \
#			-paired 0 -clip 0 \
#			-maxfastqdistr 250 \
#			-aligner {minimap2} \
#			-removeduplicates 0 \
#			-realign 0 \
#			-distrreg chr \
#			-sniffles-n -1 \
#			-svcallers {sniffles cuteSV cuteSV_pacbio npinv} \
#			-varcallers {longshot} \
#			-reports {-fastqc predictgender} \
#			-dbdir /complgen/refseq/hg38 \
#			tmp/giab_ont >& tmp/giab_ont.log
#	} else {
#		# check vs expected
#		benchmarkvars \
#			-analyses {longshot-sminimap2-HG003_NA24149_father} \
#			-regionfile tmp/giab_ont/samples/ref_HG003_NA24149_father/sreg-ref-ref-HG003_NA24149_father.tsv.zst \
#			tmp/giab_ont/compar/annot_compar-giab_ont.tsv.zst ref-ref-HG003_NA24149_father tmp/giab_ont/benchmark_HG003.tsv
#		set result {}
#		lappend result [tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x colinfo -x fastqc_report.html \
#			-x *bam.dupmetrics -x info_analysis.tsv -x *.zsti -x *.lz4i -x *.finished -x *.index -x info_analysis.tsv \
#			-x *.analysisinfo -x *.png -x *.submitting \
#			tmp/giab_ont expected/giab]
#		lappend result [diffanalysisinfo tmp/giab_ont/compar/annot_compar-*.tsv.analysisinfo expected/giab/compar/annot_compar-*.tsv.analysisinfo]
#		# lappend result [checkdiff -y --suppress-common-lines tmp/giab_ont/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics expected/giab/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics | grep -v "Started on" | grep -v bammarkduplicates2]
#		foreach file1 [glob tmp/giab_ont/compar/info_analysis.tsv tmp/giab_ont/samples/*/info_analysis.tsv] {
#			regsub ^tmp $file1 expected file2
#			lappend result [checkdiff -y --suppress-common-lines $file1 $file2 | grep -v -E {version_os|param_adapterfile|param_targetvarsfile|param_dbfiles|command|version_genomecomb}]
#		}
#		join [list_remove $result {}] \n
#	}
#} {}

testsummarize
