#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

# this script will download publically available (ont) data
# and prepare it for genomecomb tests
# By default only genomic promethion data for HG002, HG003 and HG004 and cDDNA data for HG001 are downloaded,
# and smaller test sets (limited to certain region) created for these samples.
# you can add "full" as argument to the command to also download other samples from the 
# human-pangenomics NHGRI_UCSC_panel and giab data sets
# You can give the option -d sge to distribute the processing of the data on a grid engine based cluster,
# or -d <nrthreads> to distribute over different cores/threads on a single machine. (Only processing, 
# such as alignment is distributed, not the actual downloading itself)

# The original public data used will be downloaded to the directory
# "public" in your homedir (~/public) by default.
# you can change this by setting the variable publicdir to a different location here, or
# by making another storage area first and making ~/public a link to it

# The genomecomb testdata is made in ~/genomecomb.smalltestdata ($smalltestdir)
# again, here you can use a softlink to use another storage location

# Initialise
# ==========

# see where we are, change to tests dir, and source tools.tcl there
set script [file join [pwd] [info script]]
set scriptname [file tail $script]
while 1 {
	if {[catch {set script [file join [pwd] [file readlink $script]]}]} break
}
cd [file dir $script]

source tools.tcl
set keepdir [pwd]

# default settings
# chr6:32000000-33000000 (includes hla -> too big for small test)
set rnaregions {chr1:2600000-3000000 chr2:1500000-2500000 chr6:32000000-32300000 chr10:1100000-1800000}
set wgsregions {chr1:2543891-2598871 chr2:1386718-1759764 chr6:32145542-32185752 chr10:975157-1170215}

if {![info exists argv]} {set argv {}}

set argv [job_init {*}$argv]
if {$argv eq "full"} {
	set full 1
} else {
	set full 0
}

mkdir $publicdir
cd $publicdir
job_logfile $publicdir/make_testdata_ont $publicdir make_testdata_ont.sh

# Genome in a Bottle
# ==================
# https://www.nist.gov/programs-projects/genome-bottle

# pilot
# NA12878 (HG001): Utah woman, parents are NA12891 and NA12892, genetic disease (CYP2D6 mutation)
# pilot of giab (https://jimb.stanford.edu/giab)
# also in platinum genomes
# - (https://emea.illumina.com/platinumgenomes.html): 
# - 17 member CEPH pedigree 1463 fully sequenced
# - Eberle, MA et al. (2017) A reference data set of 5.4 million phased human variants validated by genetic inheritance from sequencing a three-generation 17-member pedigree. Genome Research 27: 157-164. doi:10.1101/gr.210500.116

# main current set (Ashkenazi trio)
# HG002_NA24385_son PGP Ashkenazi Jewish son HG002 (NA24385 NIST RM8391)
# HG003_NA24149_father father HG003 (NA24149)
# HG004_NA24143_mother mother HG004 (NA24143)

# Download basic public giab data
# ===============================

mkdir $publicdir/giab

# Download truth data
# -------------------
# this is the same for SRS
cg_giab_gettruth -ref hg38 3.3.2 $publicdir/giab/truth/truth_hg38_v3.3.2
cg_giab_gettruth -ref hg38 4.2.1 $publicdir/giab/truth/truth_hg38_v4.2.1
cg_giab_gettruth -ref hg38 sv0.6 $publicdir/giab/truth/truthsv_hg38_v0.6
cg_giab_gettruth -ref hg38 hybrid $publicdir/platinum_genomes/truthset/2017-1.0/hg38/hybrid

# Download ont genomic data
# =========================

if {$full} {
	giab_getdata_job precisionfda_v2016_04 $publicdir/giab/precisionfda_v2016_04
}

# Ashkenazi trio giab_ont_ultralong
# ---------------------------------

if {$full} {
	# download data
	giab_getdata_job giab_ont_ultralong $publicdir/giab/fastqs/giab_ont_ultralong
	
	# align
	cd $publicdir/giab/giab_ont_ultralong
	foreach sample {HG002_NA24385_son HG003_NA24149_father HG004_NA24143_mother} {
		set bam $publicdir/giab/giab_ont_ultralong/$sample/map-sminimap2-${sample}_hg38.bam
		if {[file exists $bam]} continue
		map_job -method minimap2 \
			$bam \
			$::refseqdir/hg38 \
			$sample {*}[glob $publicdir/giab/giab_ont_ultralong/$sample/fastqsplit/*.fastq.gz]
		bam_index_job $publicdir/giab/giab_ont_ultralong/$sample/map-sminimap2-${sample}_hg38.bam
	}
	
	# extract regions
	foreach sample {HG002_NA24385_son HG003_NA24149_father HG004_NA24143_mother} {
		job smallfastq-$sample -deps {
			$publicdir/giab/giab_ont_ultralong/$sample/map-sminimap2-${sample}_hg38.bam
		} -targets {
			$smalltestdir/ori/giab_ont_ultralong_regions/map-sminimap2-regions_${sample}_hg38.bam
		} -vars {
			publicdir smalltestdir sample regions
		} -code {
			set dir $smalltestdir/ori/giab_ont_ultralong_regions
			set regionsbam $dir/map-sminimap2-regions_${sample}_hg38.bam
			set regionsfastq $dir/regions_${sample}_hg38.fastq.gz
			exec samtools view -b -1 -h \
				$publicdir/giab/giab_ont_ultralong/$sample/map-sminimap2-${sample}_hg38.bam \
				{*}$regions \
				> $regionsbam
			cg bam2fastq -threads 8 $regionsbam $regionsfastq
			mkdir $dir/split
			cg fastq_split -d sge -parts 20 $file $dir/split/regions_${sample}_hg38.fastq.gz
		}
	}
}

# nanopore-human-pangenomics
# ==========================
# https://humanpangenome.org/
# https://github.com/human-pangenomics/hpgp-data
# https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=NHGRI_UCSC_panel/

mkdir $publicdir/nanopore-human-pangenomics
cd $publicdir/nanopore-human-pangenomics

# 		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/nanopore/HG002_giab_ULfastqs_guppy3.2.4.fastq.gz

set todo {
	HG002 {
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/nanopore/Guppy_4.2.2/GM24385_1_Guppy_4.2.2_prom.fastq.gz
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/nanopore/Guppy_4.2.2/GM24385_2_Guppy_4.2.2_prom.fastq.gz
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/nanopore/Guppy_4.2.2/GM24385_3_Guppy_4.2.2_prom.fastq.gz
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/nanopore/Guppy_4.2.2/GM24385_4_Guppy_4.2.2_prom.fastq.gz
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/nanopore/Guppy_4.2.2/GM24385_5_Guppy_4.2.2_prom.fastq.gz
	}
	HG003 {
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG003/nanopore/Guppy_4.2.2/GM24149_1_Guppy_4.2.2_prom.fastq.gz
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG003/nanopore/Guppy_4.2.2/GM24149_2_Guppy_4.2.2_prom.fastq.gz
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG003/nanopore/Guppy_4.2.2/GM24149_3_Guppy_4.2.2_prom.fastq.gz
	}
	HG004 {
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG004/nanopore/Guppy_4.2.2/GM24143_1_Guppy_4.2.2_prom.fastq.gz
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG004/nanopore/Guppy_4.2.2/GM24143_2_Guppy_4.2.2_prom.fastq.gz
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG004/nanopore/Guppy_4.2.2/GM24143_3_Guppy_4.2.2_prom.fastq.gz
	}
}
if {$full} {
	lappend todo HG001 {
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG001/nanopore/Guppy_4.2.2/HG001_Circulomics_Guppy_4.2.2.fastq.gz
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG001/nanopore/Guppy_4.2.2/HG001_NBT2018_Guppy_4.2.2.fastq.gz
	}
	HG005 {
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG005/nanopore/Guppy_4.2.2/01_09_20_R941_GM24631_1_Guppy_4.2.2_prom.fastq.gz
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG005/nanopore/Guppy_4.2.2/01_09_20_R941_GM24631_2_Guppy_4.2.2_prom.fastq.gz
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG005/nanopore/Guppy_4.2.2/01_09_20_R941_GM24631_3_Guppy_4.2.2_prom.fastq.gz
	}
	HG006 {
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG006/nanopore/Guppy_4.2.2/01_09_20_R941_GM24694_1_Guppy_4.2.2_prom.fastq.gz
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG006/nanopore/Guppy_4.2.2/01_09_20_R941_GM24694_2_Guppy_4.2.2_prom.fastq.gz
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG006/nanopore/Guppy_4.2.2/01_09_20_R941_GM24694_3_Guppy_4.2.2_prom.fastq.gz
	}
	HG007 {
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG007/nanopore/Guppy_4.2.2/01_09_20_R941_GM24695_1_Guppy_4.2.2_prom.fastq.gz
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG007/nanopore/Guppy_4.2.2/01_09_20_R941_GM24695_2_Guppy_4.2.2_prom.fastq.gz
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG007/nanopore/Guppy_4.2.2/01_09_20_R941_GM24695_3_Guppy_4.2.2_prom.fastq.gz
	}

	foreach sample {
		HG00733 HG01109 HG01243 HG02055 HG02080 HG02723 HG03098 HG03492
	} {
		lappend todo $sample [list \
			https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/${sample}/nanopore/Guppy_4.2.2/${sample}_1_Guppy_4.2.2_prom.fastq.gz \
			https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/${sample}/nanopore/Guppy_4.2.2/${sample}_2_Guppy_4.2.2_prom.fastq.gz \
			https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/${sample}/nanopore/Guppy_4.2.2/${sample}_3_Guppy_4.2.2_prom.fastq.gz \
		]
	}
	lappend todo HG01442 {
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG01442/nanopore/08_25_21_R941_HG01442_1.fast5.tar
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG01442/nanopore/08_25_21_R941_HG01442_2.fast5.tar
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG01442/nanopore/08_25_21_R941_HG01442_3.fast5.tar
	} HG02109 {
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG02109/nanopore/08_25_21_R941_HG02109_1.fast5.tar
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG02109/nanopore/08_25_21_R941_HG02109_2.fast5.tar
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG02109/nanopore/08_25_21_R941_HG02109_3.fast5.tar
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG02109/nanopore/08_25_21_R941_HG02109_4.fast5.tar
	} HG02145 {
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG02145/nanopore/08_31_21_R941_HG02145_1.fast5.tar
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG02145/nanopore/08_31_21_R941_HG02145_2.fast5.tar
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG02145/nanopore/08_31_21_R941_HG02145_3.fast5.tar
		https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG02145/nanopore/08_31_21_R941_HG02145_4.fast5.tar
	}
}

giab_getdata_job \
	-todo $todo \
	-parts 200 \
	-align minimap2 -refseq $::refseqdir/hg38 \
	nanopore-human-pangenomics $publicdir/nanopore-human-pangenomics

# extract regions
mkdir $smalltestdir/ori/nanopore-human-pangenomics_regions
foreach sample {HG002 HG003 HG004} {
	job smallfastq-$sample -deps {
		$publicdir/nanopore-human-pangenomics/$sample/map-sminimap2-${sample}_hg38.bam
	} -targets {
		$smalltestdir/ori/nanopore-human-pangenomics_regions/$sample/map-sminimap2-regions_${sample}_hg38.bam
	} -vars {
		publicdir smalltestdir sample wgsregions
	} -code {
		set dir $smalltestdir/ori/nanopore-human-pangenomics_regions/$sample
		mkdir $dir
		set regionsbam $dir/map-sminimap2-regions_${sample}_hg38.bam
		set regionsfastq $dir/regions_${sample}_hg38.fastq.gz
		exec samtools view -b -1 -h \
			$publicdir/nanopore-human-pangenomics/$sample/map-sminimap2-${sample}_hg38.bam \
			{*}$wgsregions \
			> $regionsbam
		cg bam2fastq -threads 8 $regionsbam $regionsfastq
		mkdir $dir/fastqsplit
		cg fastq_split -parts 50 $regionsfastq $dir/fastqsplit/regions_${sample}_hg38.fastq.gz
		set o [open $dir/regions.tsv w]
		puts $o chromosome\tbegin\tend
		foreach region $wgsregions {
			puts $o [join [split $region :-] \t]
		}
		close $o
		catch {
			set src [lindex [glob $publicdir/giab/truth/truth_hg38_v4.2.1/AshkenazimTrio/$sample*/*benchmark.tsv.zst] 0]
			cg regselect $src $dir/regions.tsv | cg zst > $dir/regions_[file tail $src]
			set src [lindex [glob $publicdir/giab/truth/truth_hg38_v4.2.1/AshkenazimTrio/$sample*/*bed.tsv.zst] 0]
			cg regcommon $src $dir/regions.tsv | cg zst > $dir/regions_[file tail $src]
		}
	}
}

# tandem-genotypes paper repeat test data
# ---------------------------------------
# Tandem-genotypes: robust detection of tandem repeat expansions from long DNA reads
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6425644/#Sec15
# https://github.com/mcfrith/tandem-genotypes
# https://ddbj.nig.ac.jp/resource/sra-submission/DRA007012
# human CNBP gene

mkdir $publicdir/tandem-genotypes/plasmids
cd $publicdir/tandem-genotypes/plasmids

set runs {DRR140497 DRR140498 DRR140499 DRR140500 DRR140501 DRR140502 DRR140503 DRR140504 DRR140505 DRR140506 DRR140507}
set experiments {DRX133211 DRX133212 DRX133213 DRX133214 DRX133215 DRX133216 DRX133217 DRX133218 DRX133219 DRX133220 DRX133221}
set biosamples {SAMD00128722 SAMD00128723 SAMD00128724 SAMD00128725 SAMD00128726 SAMD00128727 SAMD00128728 SAMD00128729 SAMD00128730 SAMD00128731 SAMD00128732}
# https://ddbj.nig.ac.jp/resource/biosample/SAMD00128722
set samplenames {
SAMD00128722 CAA-15
SAMD00128723 CAA-109
SAMD00128724 GGGGCC-52
SAMD00128725 GGGGCC-21
SAMD00128726 CAG-6
SAMD00128727 CAG-18
SAMD00128728 CAG-30
SAMD00128729 CAG-30-DraIII
SAMD00128730 CAG-70
SAMD00128731 CAG-130
SAMD00128732 CCTG-45
}
foreach run $runs experiment $experiments {
	wgetfile https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA007/DRA007012/$experiment/$run.fastq.bz2 $run.fastq.bz2
	cg gzip $run.fastq.bz2
}

#SCA10-subjectA	SRR2080459
#SCA10-subjectB	SRR2081063
#SCA10-subjectC-	SRR2082412
#SCA10-subjectC-	SRR2082428

# https://www.biorxiv.org/content/biorxiv/early/2018/07/24/356931.full.pdf
#
#Human whole genome nanopore (rel3) and PacBio (SRR3197748) sequence data
#from the same individual (NA12878) were downloaded from
#(https://github.com/nanopore-wgs-consortium/NA12878) and from the SRA
#database, respectively
#
#NA19240) using PromethION was downloaded from https://www.ebi.ac.uk/ena/data/view/PRJEB26791
#
#Genomic DNA from a human patient with BAFME phenotype
#
#Plasmids containing various numbers of CAG, GGGGCC and CAA used for this
#study were generated as described elsewhere19-21 and are available upon
#request (Supplemental Table1). Sequence data of the plasmids were
#deposited in (DRA007012, Table2).
#
#SCA10 sequences from three patients with spinocerebellar ataxia 10 (MIM:
#603516)9 were downloaded from SRA (subjectA: SRR2080459, subjectB:
#SRR2081063, subjectC-1: SRR2082412, subjectC-2: SRR2082428

# tandem-genotypes followup paper -> 21 ont genomes
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7791882/
# DDBJ accession DRA009852

# RNA and cDNA
# ============

# https://github.com/nanopore-wgs-consortium/NA12878/blob/master/RNA.md
# https://s3.amazonaws.com/nanopore-human-wgs/rna/fastq/NA12878-DirectRNA_All_Guppy_4.2.2.fastq.gz

# fast5 are available, but not downloading yet
mkdir $publicdir/giab/nanopore-wgs-consortium-rna
cd $publicdir/giab/nanopore-wgs-consortium-rna
foreach {sample files} {
	HG001_NA12878_directRNA {
		NA12878-DirectRNA_All_Guppy_4.2.2.fastq.gz
		NA12878-DirectRNA_All_Guppy_4.2.2_sequencing_summary.txt.gz
	}
	HG001_NA12878_cDNA {
		NA12878-cDNA_All_Guppy_4.2.2.fastq.gz
		NA12878-cDNA_All_Guppy_4.2.2_sequencing_summary.txt.gz
	}
	HG001_NA12878_ivtRNA {
		NA12878-IVT_RNA_All_Guppy_4.2.2.fastq.gz
		NA12878-IVT_RNA_All_Guppy_4.2.2_sequencing_summary.txt.gz
	}
} {
	putsvars sample
	mkdir $publicdir/giab/nanopore-wgs-consortium-rna/$sample
	cd $publicdir/giab/nanopore-wgs-consortium-rna/$sample
	job download-nanopore-wgs-consortium-$sample -deps {
	} -targets $files -vars {
		sample files
	} -code {
		foreach file $files dir {fastq summaries} {
			wgetfile https://s3.amazonaws.com/nanopore-human-wgs/rna/$dir/$file \
				$file
		}
	}
	set fastqfile [lindex $files 0]
	mkdir fastqsplit
	fastq_split_job -parts 200 $fastqfile fastqsplit/[file tail $fastqfile]
	map_job -method minimap2 -preset splice -paired 0 -threads 8 \
		map-sminimap2_splice-${sample}_hg38.bam \
		$::refseqdir/hg38 \
		$sample {*}[jobglob fastqsplit/*.fastq.gz]
	bam_index_job map-sminimap2_splice-${sample}_hg38.bam
}

# extract regions
cd $publicdir/giab/nanopore-wgs-consortium-rna
set dir $smalltestdir/ori/nanopore-wgs-consortium-rna
mkdir $smalltestdir/ori/nanopore-wgs-consortium-rna
foreach sample {HG001_NA12878_directRNA HG001_NA12878_cDNA HG001_NA12878_ivtRNA} {
	mkdir $dir/$sample
	job rna_smallfastq-$sample -deps {
		$sample/map-sminimap2_splice-${sample}_hg38.bam
	} -targets {
		$dir/$sample/map-sminimap2_splice-regions_${sample}_hg38.bam
		$dir/$sample/regions_${sample}_hg38.fastq.gz
	} -vars {
		publicdir dir smalltestdir sample rnaregions
	} -code {
		set regionsbam $dir/$sample/map-sminimap2_splice-regions_${sample}_hg38.bam
		set regionsfastq $dir/$sample/regions_${sample}_hg38.fastq.gz
		exec samtools view -b -1 -h \
			$publicdir/giab/nanopore-wgs-consortium-rna/$sample/map-sminimap2_splice-${sample}_hg38.bam \
			{*}$rnaregions \
			> $regionsbam
		exec samtools index $regionsbam
		cg bam2fastq -threads 8 $regionsbam $regionsfastq
		mkdir $dir/$sample/splitfastq
		cg fastq_split -parts 8 $regionsfastq $dir/$sample/splitfastq/regions_${sample}_hg38.fastq.gz
	}
}

# Download methylation test data
# ==============================

job methylation_testdata -cores 6 -vars {
	smalltestdir
} -code {
	cd $::smalltestdir/ori
	wget https://f5c.page.link/f5c_na12878_test
	mv f5c_na12878_test f5c_na12878_test.tar.gz
	tar xvzf f5c_na12878_test.tar.gz
	mv chr22_meth_example ont_f5c_chr22_meth_example
	cd ont_f5c_chr22_meth_example
	mkdir fast5
	mv fast5_files single_fast5_files
	single_to_multi_fast5 -i single_fast5_files -s fast5 -n 4000
	single_to_multi_fast5 -i single_fast5_files -s fast5 -n 4000
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
	#
	# run experiment to make haplotyped bam
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
	exec cg process_project -d 6 -split 1 \
		-paired 0 -clip 0 \
		-maxfastqdistr 250 \
		-aligner {minimap2} \
		-removeduplicates 0 \
		-realign 0 \
		-distrreg chr \
		-svcallers {} \
		-varcallers longshot \
		-methcallers {} \
		-hap_bam 1 \
		-threads 6 \
		-reports {} \
		tmp/meth $::refseqdir/hg19 >& tmp/meth.log
	# after analysis put haplotytyped bam back in ori
	cp -alf tmp/meth/samples/methtest/map-hlongshot-sminimap2-methtest.bam ori/ont_f5c_chr22_meth_example/
	cp -alf tmp/meth/samples/methtest/map-hlongshot-sminimap2-methtest.bam.bai ori/ont_f5c_chr22_meth_example/
}

job_wait
