#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl
set keepdir [pwd]

set test_cleantmp 0

# tests
# =====

lappend dopts -threads 1
set runopts {-stack 1 -v 2}
set dopts {-d 8 -threads 8}
set dopts {-d 2 -threads 2}

test cg_flames {rna flames} {
	cd ~/genomecomb.smalltestdata
	file delete -force tmp/rna_flames
	file mkdir tmp/rna_flames/samples
	foreach sample {
		HG001_NA12878_cDNA  HG001_NA12878_directRNA  HG001_NA12878_ivtRNA
	} {
		file mkdir tmp/rna_flames/samples/$sample/fastq
		foreach file [glob ori/nanopore-wgs-consortium-rna/$sample/splitfastq/*fastq.gz] {
			mklink $file tmp/rna_flames/samples/$sample/fastq/[file tail $file]
		}
		foreach file [glob ori/nanopore-wgs-consortium-rna/$sample/splitfastq/*.bam*]
			mklink $file tmp/rna_flames/samples/$sample/[file tail $file]
		}
#			foreach version {3.3.2 4.2.1} {
#				if {$version eq "3.3.2"} {
#					set samplename truth_$sample
#				} else {
#					set samplename truth[string_change $version {. {}}]_$sample
#				}
#				set sampledir tmp/rna_flames/samples/$samplename
#				mkdir $sampledir
#				set vcf [glob public/giab/truth/truth_hg38_v$version/*/${sample}_hg38/*.vcf.gz]
#				set region [glob public/giab/truth/truth_hg38_v$version/*/${sample}_hg38/*.bed.tsv.zst]
#				mklink $vcf $sampledir/var-truth-truth-$samplename.vcf.gz 1
#				mklink [file root [gzroot $vcf]].tsv.zst $sampledir/var-truth-truth-$samplename.tsv.zst 1
#				mklink $region $sampledir/sreg-truth-truth-$samplename.tsv.zst 1
#			}
	}
	# run

	set varcallers {longshot}
	exec cg process_project {*}$::runopts {*}$::dopts \
		-split 1 \
		-paired 0 -clip 0 \
		-maxfastqdistr 250 \
		-aligner {minimap2_splice} \
		-removeduplicates 0 \
		-realign 0 \
		-distrreg chr \
		-svcallers {} \
		-flair 1 -flair-compar 1 \
		-varcallers $varcallers \
		-reports {-fastqc predictgender} \
		-dbdir $::refseqdir/hg38 \
		tmp/rna_flames >& tmp/rna_flames.log

	# check vs expected
	benchmarkvars \
		-analyses {longshot-sminimap2-HG003_NA24149_father} \
		-regionfile tmp/rna_flames/samples/truth421_HG003_NA24149_father/sreg-truth-truth-truth421_HG003_NA24149_father.tsv.zst \
		tmp/rna_flames/compar/annot_compar-rna_flames.tsv.zst truth-truth-truth421_HG003_NA24149_father tmp/rna_flames/benchmark_HG003.tsv
	set result {}
	lappend result [tsvdiff -q 1 -x *log_jobs -x *.bam -x *.bai -x colinfo -x fastqc_report.html \
		-x *bam.dupmetrics -x info_analysis.tsv -x *.zsti -x *.lz4i -x *.finished -x *.index -x info_analysis.tsv \
		-x *.analysisinfo -x *.png -x *.submitting \
		tmp/rna_flames expected/giab]
	lappend result [diffanalysisinfo tmp/rna_flames/compar/annot_compar-*.tsv.analysisinfo expected/giab/compar/annot_compar-*.tsv.analysisinfo]
	# lappend result [checkdiff -y --suppress-common-lines tmp/rna_flames/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics expected/giab/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics | grep -v "Started on" | grep -v bammarkduplicates2]
	foreach file1 [glob tmp/rna_flames/compar/info_analysis.tsv tmp/rna_flames/samples/*/info_analysis.tsv] {
		regsub ^tmp $file1 expected file2
		lappend result [checkdiff -y --suppress-common-lines $file1 $file2 | grep -v -E {version_os|param_adapterfile|param_targetvarsfile|param_dbfiles|command|version_genomecomb}]
	}
	join [list_remove $result {}] \n

} {}


testsummarize


