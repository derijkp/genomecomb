#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl
set keepdir [pwd]

set download 0
set test_cleantmp 0

# Download database
# =================
if {[get download 0]} {
	
	# NA12878: Utah woman, parents are NA12891 and NA12892, genetic disease (CYP2D6 mutation)
	# pilot of giab (https://jimb.stanford.edu/giab)
	# also in platinum genomes
	# - (https://emea.illumina.com/platinumgenomes.html): 
	# - 17 member CEPH pedigree 1463 fully sequenced
	# - Eberle, MA et al. (2017) A reference data set of 5.4 million phased human variants validated by genetic inheritance from sequencing a three-generation 17-member pedigree. Genome Research 27: 157-164. doi:10.1101/gr.210500.116
	# 
	# AshkenazimTrio
	# NA24385 (HG002_NA24385_son), parents are HG003_NA24149_father and HG004_NA24143_mother
	# PGP genomes

	# giab data
	# https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/
	# ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/

	cd ~/genomecomb_giab_testdata/public
	set basedir [pwd]/giab
	file mkdir $basedir

	# giab truth set
	file mkdir $basedir/truthset/truth_hg38_giab3.3.2
	cd $basedir/truthset/truth_hg38_giab3.3.2
	# set baseurl ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp
	set baseurl ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab
	foreach url [list \
		$baseurl/release/NA12878_HG001/NISTv3.3.2/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz \
		$baseurl/release/NA12878_HG001/NISTv3.3.2/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz.tbi \
		$baseurl/release/NA12878_HG001/NISTv3.3.2/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed \
		$baseurl/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/README_NISTv3.3.2.txt \
		$baseurl/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh38/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed \
		$baseurl/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh38/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz \
		$baseurl/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh38/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz.tbi \
		$baseurl/release/AshkenazimTrio/HG003_NA24149_father/NISTv3.3.2/GRCh38/HG003_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_highconf.vcf.gz \
		$baseurl/release/AshkenazimTrio/HG003_NA24149_father/NISTv3.3.2/GRCh38/HG003_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_highconf.vcf.gz.tbi \
		$baseurl/release/AshkenazimTrio/HG003_NA24149_father/NISTv3.3.2/GRCh38/HG003_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed \
		$baseurl/release/AshkenazimTrio/HG004_NA24143_mother/NISTv3.3.2/GRCh38/HG004_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_highconf.vcf.gz \
		$baseurl/release/AshkenazimTrio/HG004_NA24143_mother/NISTv3.3.2/GRCh38/HG004_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_highconf.vcf.gz.tbi \
		$baseurl/release/AshkenazimTrio/HG004_NA24143_mother/NISTv3.3.2/GRCh38/HG004_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed
	] {
		set file [file tail $url]
		set ext [file extension [gzroot $file]]
		if {$ext eq ".vcf"} {
			if {[file exists [file root [gzroot $file]].tsv.zst]} {
				puts "Skipping vcf file $file: done"
				continue
			}
			puts "Converting vcf file $file to [file root [gzroot $file]].tsv.zst"
			exec wget -c $url 2>@ stderr >@ stdout
			cg vcf2tsv $file [gzroot $file].temp.tsv.zst
			file rename -force [gzroot $file].temp.tsv.zst [gzroot $file].tsv.zst
		} elseif {$ext eq ".bed"} {
			if {[file exists [gzroot $file].tsv.zst]} {
				puts "Skipping bed file $file: done"
				continue
			}
			exec wget -c $url 2>@ stderr >@ stdout
			puts "Converting bed file $file to [gzroot $file].tsv.zst"
			cg bed2tsv $file [gzroot $file].temp.tsv.zst
			file rename -force [gzroot $file].temp.tsv.zst [gzroot $file].tsv.zst
		}
	}

	# download publically available sequence data
	set giabbaseurl ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data
	#
	# NA12878
	file mkdir $basedir/NA12878/Garvan_NA12878_HG001_HiSeq_Exome
	cd $basedir/NA12878/Garvan_NA12878_HG001_HiSeq_Exome
	exec wget -c $giabbaseurl/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/*1.fastq.gz 2>@ stderr >@ stdout
	#
	# AshkenazimTrio
	set giabsamples {
		HG002_NA24385_son HG003_NA24149_father HG004_NA24143_mother
	}
	foreach sample $giabsamples {
		file mkdir $basedir/AshkenazimTrio/$sample
		cd $basedir/AshkenazimTrio/$sample
		exec wget -c -r -nH --cut-dirs=5 $giabbaseurl/AshkenazimTrio/$sample/OsloUniversityHospital_Exome 2>@ stderr >@ stdout
	}
	# fastq
	foreach sample $giabsamples {
		cd $basedir/AshkenazimTrio/$sample/OsloUniversityHospital_Exome
		set file [glob *.bam]
		set fq1 [file root $file]_R1.fastq.gz
		set fq2 [file root $file]_R2.fastq.gz
		if {![file exists $fq1] || ![file exists $fq2]} {
			puts "Skipping $sample fastqs: already present"
			puts "Making $fq1 and $fq2"
			exec cg bam2fastq $file temp_$fq1 temp_$fq2
			file rename temp_$fq1 $fq1
			file rename temp_$fq2 $fq2
		}
		if {[file exists fastqsplit]} {
			puts "Skipping $sample/fastqsplit: already present"
			continue
		}
		file mkdir fastqsplit.temp
		puts "Splitting $fq1"
		exec cg fastq_split -stack 1 -v 2 -d sge -parts 200 $fq1 fastqsplit.temp/[file tail $fq1] >@ stdout 2>@stderr
		puts "Splitting $fq2"
		exec cg fastq_split -stack 1 -v 2 -d sge -parts 200 $fq2 fastqsplit.temp/[file tail $fq2] >@ stdout 2>@stderr
	}
	foreach sample $giabsamples {
		cd $basedir/AshkenazimTrio/$sample/OsloUniversityHospital_Exome
		file rename fastqsplit.temp fastqsplit
	}
	# experiments from sra
	file mkdir $basedir/sra
	cd $basedir/sra
	foreach {name var runid} {
		NA12878	30x SRR098401
		NA24385 30x ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR296/009/SRR2962669
		NA24631 30x ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR296/003/SRR2962693
	} {
		if {[regexp ^ftp: $runid]} {
			set url $runid
			set runid [lindex [file split $runid] end]
		} else {
			set url ftp://ftp.sra.ebi.ac.uk/vol1/fastq/[string range $runid 0 5]/${runid}
		}
		set dir ${runid}_${var}_${name}_exome
		putsvars dir
		set dest $basedir/sra/$dir
		file mkdir $dest
		cd $dest
		wgetfile $url/${runid}_1.fastq.gz
		wgetfile $url/${runid}_2.fastq.gz
		#
		# split fastqs
		if {[file exists fastqsplit]} {
			puts "Skipping $dest/fastqsplit: already present"
			continue
		}
		file mkdir fastqsplit.temp
		set files [glob *.fastq.gz]
		foreach file $files {
			putsvars file
			exec cg fastq_split -stack 1 -v 2 -d sge -parts 200 $file fastqsplit.temp/[file tail $file] >@ stdout 2>@stderr
		}
		file rename fastqsplit.temp fastqsplit
	}
}

# extra code
# ==========
# go and take from public_genomes.proj/procedure-rungenomes.txt

# tests
# =====

test process_giab {process_giab NA12878} {
	cd ~/genomecomb_giab_testdata
	file delete -force tmp/giabexome
	file mkdir tmp/giabexome/samples
	foreach sample {
		NA12878
	} {
		file mkdir tmp/giabexome/samples/$sample/fastq
		foreach file [glob public/giab/sra/*_${sample}_exome/fastqsplit/*fastq.gz] {
			mklink $file tmp/giabexome/samples/$sample/fastq/[file tail $file] 1
		}
	}
	mkdir tmp/giabexome/samples/truth_NA12878
	mklink public/giab/truthset/truth_hg38_giab3.3.2/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.tsv.zst \
		tmp/giabexome/samples/truth_NA12878/var-truth-truth-truth_NA12878.tsv.zst 1
	mklink public/giab/truthset/truth_hg38_giab3.3.2/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed.tsv.zst \
		tmp/giabexome/samples/truth_NA12878/sreg-truth-truth-truth_NA12878.tsv.zst 1
	# run
	exec cg process_project -stack 1 -v 2 -d sge -split 1 \
		-threads 8 -distrreg chr -varcallers {gatkh strelka} -svcallers {manta lumpy} \
		-reports {all} \
		-cleanup 0 \
		-dbdir /complgen/refseq/hg38 \
		tmp/giabexome >& tmp/giabexome.log
	# check vs expected
	set result {}
	lappend result [tsvdiff -q 1 \
		-x *.bai -x *.zsti -x *.lz4i -x *.png \
		-x *log_jobs -x *.submitting -x *.finished -x info_analysis.tsv -x *.analysisinfo \
		-x *.index -x colinfo \
		-x fastqc_report.html -x *bam.dupmetrics \
		tmp/giabexome expected/giabexome]
	lappend result [diffanalysisinfo tmp/giabexome/compar/annot_compar-*.tsv.analysisinfo expected/giabexome/compar/annot_compar-*.tsv.analysisinfo]
	foreach file1 [glob tmp/giabexome/compar/info_analysis.tsv tmp/giabexome/samples/*/info_analysis.tsv] {
		regsub ^tmp $file1 expected file2
		lappend result [checkdiff \
			-I version_os -I param_adapterfile -I param_targetvarsfile -I param_dbfiles -I command -I version_genomecomb \
			$file1 $file2]
	}
	join [list_remove $result {}] \n
} {}

test process_giab {process_giab trio} {
	cd ~/genomecomb_giab_testdata
	file delete -force tmp/giabexome_trio
	file mkdir tmp/giabexome_trio/samples
	foreach sample {
		HG002_NA24385_son	HG003_NA24149_father	HG004_NA24143_mother
	} {
		file mkdir tmp/giabexome_trio/samples/$sample/fastq
		foreach file [glob public/giab/AshkenazimTrio/${sample}/OsloUniversityHospital_Exome/fastqsplit/*fastq.gz] {
			mklink $file tmp/giabexome_trio/samples/$sample/fastq/[file tail $file] 1
		}
		catch {mkdir tmp/giabexome_trio/samples/truth_$sample}
		mklink [lindex [glob public/giab/truthset/truth_hg38_giab3.3.2/[lindex [split $sample _] 0]_*highconf*.vcf.tsv.zst] 0] \
			tmp/giabexome_trio/samples/truth_$sample/var-ref-ref-truth_${sample}.tsv.zst
		mklink [lindex [glob public/giab/truthset/truth_hg38_giab3.3.2/[lindex [split $sample _] 0]_*highconf*.bed.tsv.zst] 0] \
			tmp/giabexome_trio/samples/truth_$sample/sreg-ref-ref-truth_${sample}.tsv.zst
	}
	# run
	exec cg process_project -stack 1 -v 2 -d sge -split 1 \
		-threads 8 -distrreg chr -varcallers {gatkh strelka} -svcallers {manta lumpy} \
		-reports {all} \
		-cleanup 0 \
		-dbdir /complgen/refseq/hg38 \
		tmp/giabexome_trio >& tmp/giabexome_trio.log
	# check vs expected
	cg benchmarkvars \
		-analyses {strelka-rdsbwa-HG002_NA24385_son gatkh-rdsbwa-HG002_NA24385_son} \
		-regionfile tmp/giabexome_trio/samples/truth_HG002_NA24385_son/sreg-ref-ref-truth_HG002_NA24385_son.tsv.zst \
		tmp/giabexome_trio/compar/annot_compar-giabexome_trio.tsv.zst ref-ref-truth_HG002_NA24385_son tmp/giabexome_trio/benchmark-HG002_NA24385_son.tsv
	set result {}
	lappend result [tsvdiff -q 1 \
		-x *.bai -x *.zsti -x *.lz4i -x *.png \
		-x *log_jobs -x *.submitting -x *.finished -x info_analysis.tsv -x *.analysisinfo \
		-x *.index -x colinfo \
		-x fastqc_report.html -x *bam.dupmetrics \
		tmp/giabexome_trio expected/giabexome_trio]
	lappend result [diffanalysisinfo tmp/giabexome_trio/compar/annot_compar-*.tsv.analysisinfo expected/giabexome_trio/compar/annot_compar-*.tsv.analysisinfo]
	foreach file1 [glob tmp/giabexome_trio/compar/info_analysis.tsv tmp/giabexome_trio/samples/*/info_analysis.tsv] {
		regsub ^tmp $file1 expected file2
		lappend result [checkdiff \
			-I version_os -I param_adapterfile -I param_targetvarsfile -I param_dbfiles -I command -I version_genomecomb \
			$file1 $file2]
	}
	join [list_remove $result {}] \n
} {}

testsummarize

# giabexome on dipsy (288 threads) -> wall time 2h12min -> all cluster max load -> 14.02 min