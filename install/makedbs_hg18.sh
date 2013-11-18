#!/bin/sh
# the next line restarts using tclsh \
exec cg source "$0" "$@"

set build hg18
set argv [job_init {*}$argv]
if {[llength $argv]} {
	set dest [lindex $argv 0]
} else {
	set dest /complgen/refseq/
}

# download hg18
# =============
#
file mkdir ${dest}/${build}
cd ${dest}/${build}
file mkdir extra
file mkdir ${dest}/tmp/hg18

job_logdir log_jobs

# download genome
job genome_${build} -vars build -targets {genome_${build}.ifas extra/reg_${build}_fullgenome.tsv} -code {
	cg downloadgenome ${build} genome_${build}.ifas
	file rename -force reg_genome_${build}.tsv extra/reg_${build}_fullgenome.tsv
}

job genome_${build}_cindex -deps {genome_${build}.ifas} -targets {genome_${build}.ssa} -code {
	cg make_genomecindex $dep
}

job reg_${build}_sequencedgenome -vars {dest build} -deps {genome_${build}.ifas} -targets {extra/reg_${build}_sequencedgenome.tsv} -code {
	cg calcsequencedgenome $dest ${build}
}

# region databases (ucsc)
# you can explicitely download info on the databases using:
# cg downloaddbinfo ${dest} ${build} simpleRepeat microsat rmsk genomicSuperDups chainSelf

# collapse regions
foreach db {
	cytoBand evofold gwasCatalog microsat oreganno rmsk simpleRepeat targetScanS tfbsConsSites 
	tRNAs wgRna vistaEnhancers gad
	phastConsElements28way phastConsElements28wayPlacMammal phastConsElements44way
} {
	job reg_${build}_$db -targets {reg_${build}_${db}.tsv} -vars {dest build db} -code {
		cg downloaddb $dest ${build} $db
		cg collapseoverlap ucsc_${build}_${db}.tsv
		file delete ucsc_${build}_${db}.tsv
	}
}

# join regions
foreach db {
	chainSelf dgvMerged genomicSuperDups
} {
	job reg_${build}_$db -targets {reg_${build}_${db}.tsv} -vars {dest build db} -code {
		cg downloaddb $dest ${build} $db
		cg regjoin ucsc_${build}_${db}.tsv > reg_${build}_${db}.tsv.temp
		file delete ucsc_${build}_${db}.tsv
		file rename -force reg_${build}_${db}.tsv.temp reg_${build}_${db}.tsv
	}
}

# direct
foreach db {
	firstEF rnaGene
} {
	job reg_${build}_$db -targets {reg_${build}_${db}.tsv} -vars {dest build db} -code {
		cg downloaddb $dest ${build} $db
		file rename -force ucsc_${build}_${db}.tsv reg_${build}_${db}.tsv
	}
}

job reg_${build}_gwasCatalog -vars {build dest} -deps {ucsc_${build}_gwasCatalog.tsv} -targets {reg_${build}_gwasCatalog.tsv} -code {
	cg downloaddb $dest ${build} gwasCatalog
	cg select \
		-f {chrom start end trait pValue pubMedID name bin author pubDate journal title initSample replSample region genes riskAllele riskAlFreq pValueDesc orOrBeta ci95 platform cnv} \
		-nh {chrom start end name score pubMedID dbsnp bin author pubDate journal title initSample replSample region genes riskAllele riskAlFreq pValueDesc orOrBeta ci95 platform cnv} \
		ucsc_${build}_gwasCatalog.tsv ucsc_${build}_gwasCatalog.tsv.temp
	cg collapseoverlap ucsc_${build}_gwasCatalog.tsv.temp
	file rename -force reg_${build}_gwasCatalog.tsv.temp $target
	file delete ucsc_${build}_gwasCatalog.tsv
	file delete ucsc_${build}_gwasCatalog.tsv.temp
}


# other databases
foreach db {kgXref refLink} {
	job other_${build}_$db -vars {dest build db} -targets {other_${build}_${db}.tsv} -code {
		cg downloaddb $dest ${build} $db
		file rename -force ucsc_${build}_${db}.tsv other_${build}_${db}.tsv
	}
}

foreach db {
	rmsk simpleRepeat
} {
	job maketabix_${build}_$db -deps {reg_${build}_${db}.tsv} -targets {reg_${build}_${db}.tsv.gz.tbi reg_${build}_${db}.tsv.gz} -vars {build db} -code {
		cg maketabix reg_${build}_${db}.tsv
		exec gunzip -c reg_${build}_${db}.tsv.gz > reg_${build}_${db}.tsv
	}
}

# 1000 genomes
job 1000g -targets {$dest/hg18/var_hg18_1000gCHBxJPT.tsv $dest/hg18/var_hg18_1000gCEU.tsv $dest/hg18/var_hg18_1000gYRI.tsv} -vars {dest} -code {
	cd ${dest}/hg18
	cg downloaddb ${dest} hg18 1000g
}

# 1000glow (hg19)
file mkdir ${dest}/hg19
job 1000glow -targets {${dest}/hg19/var_hg19_1000glow.tsv} -vars {dest} -code {
	cd $dest/hg19
	cg downloaddb ${dest} hg19 1000glow
}

job 1000glow_liftover -deps {${dest}/hg19/var_hg19_1000glow.tsv} -targets {var_${build}_1000glow.tsv} -vars {dest build} -code {
	cg liftover ${dest}/hg19/var_hg19_1000glow.tsv ${dest}/liftover/hg19ToHg18.over.chain var_${build}_1000glow.tsv
}

# dbsnp
job dbsnp130 -targets {var_${build}_snp130.tsv} -vars {dest build} -code {
	cg downloaddb ${dest} ${build} snp130
}

# dbsnp138 (hg19)
job dbsnp138 -targets {${dest}/hg19/var_hg19_snp138.tsv} -vars {dest} -code {
	cd $dest/hg19
	cg downloaddb ${dest} hg19 snp138
}

# dbsnp138Common (hg19)
job dbsnp138Common -targets {${dest}/hg19/var_hg19_snp138Common.tsv} -vars {dest build} -code {
	cd $dest/hg19
	cg downloaddb ${dest} hg19 snp138Common
}

# liftover
job dbsnp138lift -deps {${dest}/hg19/var_hg19_snp138.tsv} -targets {var_${build}_snp138lift.tsv} -vars {dest build} -code {
	cg liftover ${dest}/hg19/var_hg19_snp138.tsv ${dest}/liftover/hg19ToHg18.over.chain ${dest}/${build}/var_${build}_snp138lift.tsv
	file delete ${dest}/${build}/var_${build}_snp138lift.tsv.unmapped
}

job dbsnp138Commonlift -deps {${dest}/hg19/var_hg19_snp138Common.tsv} -targets {var_${build}_snp138Commonlift.tsv} -vars {dest build} -code {
	cg liftover ${dest}/hg19/var_hg19_snp138Common.tsv ${dest}/liftover/hg19ToHg18.over.chain ${dest}/${build}/var_${build}_snp138Commonlift.tsv
	file delete ${dest}/${build}/var_${build}_snp138Commonlift.tsv.unmapped
}

foreach db {
	snp130 snp138lift snp138Commonlift
} {
	job maketabix_${build}_$db -deps {reg_${build}_${db}.tsv} -targets {reg_${build}_${db}.tsv.gz.tbi reg_${build}_${db}.tsv.gz} -vars {dest build db} -code {
		cg maketabix reg_${build}_${db}.tsv
		# exec gunzip -c var_${build}_${db}.tsv.gz > var_${build}_${db}.tsv
		cg select -f 'chrom start end type ref alt name freq' ${dest}/${build}/var_${build}_${db}.tsv.gz ${dest}/${build}/var_${build}_${db}.tsv
	}
}

# genes
foreach db {
	refGene ensGene knownGene genscan augustusAbinitio acembly wgEncodeGencodeManualV3 wgEncodeGencodeAutoV3
} {
	job gene_${build}_$db -targets {gene_${build}_${db}.tsv gene_${build}_${db}.tsv.gz.tbi gene_${build}_${db}.tsv.gz} -vars {dest build db} -code {
	        cg downloaddb ${dest} ${build} $db
	        cg select -s - ucsc_${build}_${db}.tsv gene_${build}_${db}.tsv.temp
		file rename -force gene_${build}_${db}.tsv.temp gene_${build}_${db}.tsv
	        file rename -force reg_${build}_${db}.info gene_${build}_${db}.info
	        cg maketabix gene_${build}_${db}.tsv
	        exec gunzip -c gene_${build}_${db}.tsv.gz > gene_${build}_${db}.tsv.temp
		file rename -force gene_${build}_${db}.tsv.temp gene_${build}_${db}.tsv
		file delete ucsc_${build}_${db}.tsv
	}
}

# homopolymer
job reg_${build}_homopolymer -deps {genome_${build}.ifas} -targets {reg_${build}_homopolymer.tsv reg_${build}_homopolymer.tsv.gz reg_${build}_homopolymer.tsv.gz.tbi reg_${build}_homopolymer.tsv.opt} -vars {dest build db} -code {
	cg extracthomopolymers genome_${build}.ifas > reg_${build}_homopolymer.tsv.temp
	file rename -force reg_${build}_homopolymer.tsv.temp reg_${build}_homopolymer.tsv
        cg maketabix reg_${build}_homopolymer.tsv
        gunzip reg_${build}_homopolymer.tsv.gz reg_${build}_homopolymer.tsv
	file_write reg_${build}_homopolymer.tsv.opt "fields\t{base size}\n"
}

# mirbase (hg19)
job reg_hg19_mirbase -targets {$dest/hg19/reg_hg19_mirbase.tsv $dest/hg19/reg_hg19_mirbase.info} -vars {dest build db} -code {
	cd $dest/hg19
	exec -ignorestderr wget -c --tries=45 --directory-prefix=${dest}/tmp/hg19 ftp://mirbase.org/pub/mirbase/20/genomes/hsa.gff2
	cg gff2sft ${dest}/tmp/hg19/hsa.gff2 ${dest}/tmp/hg19/reg_hg19_mirbase.tsv.temp
	cg select -s - ${dest}/tmp/hg19/reg_hg19_mirbase.tsv.temp ${dest}/tmp/hg19/reg_hg19_mirbase.tsv.temp2
	file rename -force ${dest}/tmp/hg19/reg_hg19_mirbase.tsv.temp2 reg_hg19_mirbase.tsv
	exec -ignorestderr wget -c ftp://mirbase.org/pub/mirbase/20/README
	file rename -force README reg_hg19_mirbase.info
}

# mirbase hg18 liftover
job reg_hg18_mirbase_liftover -deps {${dest}/hg19/reg_hg19_mirbase.tsv ${dest}/hg19/reg_hg19_mirbase.info} -targets {reg_${build}_mirbase.tsv reg_${build}_mirbase.info} -vars {dest build db} -code {
	cg liftover ${dest}/hg19/reg_hg19_mirbase.tsv ${dest}/liftover/hg19ToHg18.over.chain ${dest}/tmp/${build}/reg_${build}_mirbase.tsv.temp2
	file rename -force ${dest}/tmp/${build}/reg_${build}_mirbase.tsv.temp2 reg_${build}_mirbase.tsv
	file copy ${dest}/hg19/reg_hg19_mirbase.info reg_${build}_mirbase.info
}

# encode
# enc_transcription
job enc_transcription -targets {reg_${build}_wgEncodelogRnaSeq.tsv} -vars {dest build} -code {
	cd ${dest}/tmp/${build}
	set tables {wgEncodeCaltechRnaSeqRawSignalRep1Gm12878CellLongpolyaBb12x75 wgEncodeCaltechRnaSeqRawSignalRep1H1hescCellPapBb2R2x75 wgEncodeCaltechRnaSeqRawSignalRep2Hepg2CellPapBb2R2x75 wgEncodeCaltechRnaSeqRawSignalRep1HuvecCellPapBb2R2x75 wgEncodeCaltechRnaSeqRawSignalRep1K562CellLongpolyaBb12x75 wgEncodeCaltechRnaSeqRawSignalRep1NhekCellPapBb2R2x75}
	cg downloaddb ${dest}/tmp ${build} {*}$tables
	set todo {}
	foreach table $tables {
		cg ucscwiggle2reg -n 0.01 -p 1 -f {log10($value)} ucsc_${build}_$table.tsv
		lappend todo reg_ucsc_${build}_$table.tsv
	}
	cg collapseoverlap -o ${dest}/${build}/reg_${build}_wgEncodelogRnaSeq.tsv {*}$todo
}

# make enc_transcription info file
job enc_transcription_info -targets {reg_${build}_wgEncodeCaltechRnaSeq.info} -vars {dest build} -code {
	cd ${dest}/tmp/${build}
	cg downloaddbinfo ${dest}/tmp ${build} wgEncodeCaltechRnaSeq
	set f [open ${dest}/tmp/${build}/reg_${build}_wgEncodeCaltechRnaSeq.info]
	set o [open ${dest}/${build}/reg_${build}_wgEncodeCaltechRnaSeq.info.temp w]
	set line [gets $f]
	puts $o $line
	puts $o "\n\n== Agregation info ==\n\nThis file combines the data from 6 lines"
	puts $o "score contains the logarithm of the highest score, with a precision of 1 digit after the decimal point."
	puts $o "num contains the number of lines for which the original value is >= 0.01 (thus log >= -2.0)\n"
	fcopy $f $o
	close $f
	close $o
	file rename -force ${dest}/${build}/reg_${build}_wgEncodeCaltechRnaSeq.info.temp ${dest}/${build}/reg_${build}_wgEncodeCaltechRnaSeq.info
}

foreach {jobname resultname infosrc tables} {
	enc_H3K4Me1 wgEncodeH3k4me1 wgEncodeRegMarkEnhH3k4me1 {wgEncodeBroadChipSeqSignalGm12878H3k4me1 wgEncodeBroadChipSeqSignalH1hescH3k4me1 wgEncodeBroadChipSeqSignalHmecH3k4me1 wgEncodeBroadChipSeqSignalHsmmH3k4me1 wgEncodeBroadChipSeqSignalHuvecH3k4me1 wgEncodeBroadChipSeqSignalK562H3k4me1 wgEncodeBroadChipSeqSignalNhekH3k4me1 wgEncodeBroadChipSeqSignalNhlfH3k4me1}
	enc_H3K27Ac wgEncodeH3k27ac wgEncodeRegMarkEnhH3k27ac {wgEncodeBroadChipSeqSignalGm12878H3k27ac wgEncodeBroadChipSeqSignalHepg2H3k27ac wgEncodeBroadChipSeqSignalHmecH3k27ac wgEncodeBroadChipSeqSignalHsmmH3k27ac wgEncodeBroadChipSeqSignalHuvecH3k27ac wgEncodeBroadChipSeqSignalK562H3k27ac wgEncodeBroadChipSeqSignalNhekH3k27ac wgEncodeBroadChipSeqSignalNhlfH3k27ac}
	enc_H3K4Me3	wgEncodeH3k4me3 wgEncodeRegMarkPromoter   {wgEncodeBroadChipSeqSignalGm12878H3k4me3 wgEncodeBroadChipSeqSignalH1hescH3k4me3 wgEncodeBroadChipSeqSignalHepg2H3k4me3 wgEncodeBroadChipSeqSignalHmecH3k4me3 wgEncodeBroadChipSeqSignalHsmmH3k4me3 wgEncodeBroadChipSeqSignalHuvecH3k4me3 wgEncodeBroadChipSeqSignalK562H3k4me3 wgEncodeBroadChipSeqSignalNhekH3k4me3 wgEncodeBroadChipSeqSignalNhlfH3k4me3}
	
} {
	# make database
	job $jobname -targets {${dest}/${build}/reg_${build}_$resultname.tsv} -vars {dest build tables} -code {
		cd ${dest}/tmp/${build}
		cg downloaddb ${dest}/tmp ${build} {*}$tables
		set todo {}
		foreach table $tables {
			cg ucscwiggle2reg -n 10 -p 0 -f {5*(($value+4)/5)} ucsc_${build}_$table.tsv
			lappend todo reg_ucsc_${build}_$table.tsv
		}
		cg collapseoverlap -o $target {*}$todo
	}
	# make info file
	job ${jobname}_info -targets {${dest}/${build}/reg_${build}_$resultname.info} -vars {dest build infosrc tables} -code {
		cd ${dest}/tmp/${build}
		cg downloaddbinfo ${dest}/tmp ${build} $infosrc
		set f [open ${dest}/tmp/${build}/reg_${build}_$infosrc.info]
		set o [open $target.temp w]
		set line [gets $f]
		puts $o $line
		puts $o "\n\n== Agregation info ==\n"
		puts $o "This file combines the data from [llength $tables] cell lines"
		puts $o "score contains the highest score rounded up the next 5 fold."
		puts $o "num contains the number of cell lines for which the score is >= 10\n"
		fcopy $f $o
		close $f
		close $o
		file rename -force $target.temp $target
	}
}

# DNase Clusters and Txn Factor ChIP
job enc_RegDnaseClustered -targets {reg_${build}_wgEncodeRegDnaseClustered.tsv reg_${build}_wgEncodeRegDnaseClustered.info} -vars {dest build} -code {
	cd ${dest}/tmp/${build}
	cg downloaddb ${dest}/tmp ${build} wgEncodeRegDnaseClustered
	cg collapseoverlap ucsc_${build}_wgEncodeRegDnaseClustered.tsv
	file rename -force reg_${build}_wgEncodeRegDnaseClustered.info ${dest}/${build}
	file rename -force reg_${build}_wgEncodeRegDnaseClustered.tsv ${dest}/${build}
}

job enc_RegTfbsClustered -targets {reg_${build}_wgEncodeRegTfbsClustered.tsv reg_${build}_wgEncodeRegTfbsClustered.info} -vars {dest build} -code {
	cd ${dest}/tmp/${build}
	cg downloaddb ${dest}/tmp ${build} wgEncodeRegTfbsClustered
	cg collapseoverlap ucsc_${build}_wgEncodeRegTfbsClustered.tsv
	file rename -force reg_${build}_wgEncodeRegTfbsClustered.info ${dest}/${build}
	file rename -force reg_${build}_wgEncodeRegTfbsClustered.tsv ${dest}/${build}
}

# link local data in dir
catch {
	cplinked {*}[glob ../hg18-local/*] .
}

# extra

# dbNSFPzip (hg19)
file mkdir ${dest}/tmp/hg19
file mkdir ${dest}/hg19
job pre_var_hg19_dbnsfp -targets {${dest}/tmp/hg19/pre_var_hg19_dbnsfp} -skip {extra/var_hg19_dbnsfp.tsv extra/var_hg19_dbnsfp.tsv.opt} -vars {dest} -code {
	cd ${dest}/tmp/hg19
	if {![file exists dbNSFP2.0b4.zip]} {
		exec -ignorestderr wget -c http://dbnsfp.houstonbioinformatics.org/dbNSFPzip/dbNSFPv2.1.zip
	}
	exec unzip -o dbNSFPv2.1.zip >@ stdout 2>@ stderr
	file_write pre_var_dbnsfp.tsv [join {
		chromosome pos ref alt aaref aaalt hg18_pos genename Uniprot_acc Uniprot_id Uniprot_aapos
		Interpro_domain cds_strand refcodon SLR_test_statistic codonpos fold-degenerate Ancestral_allele
		Ensembl_geneid Ensembl_transcriptid aapos 
		SIFT_score SIFT_score_converted SIFT_pred
		Polyphen2_HDIV_score Polyphen2_HDIV_pred Polyphen2_HVAR_score Polyphen2_HVAR_pred 
		LRT_score LRT_score_converted LRT_pred
		MutationTaster_score MutationTaster_score_converted MutationTaster_pred 
		FATHMM_score FATHMM_score_converted FATHMM_pred
		GERP_NR GERP_RS PhyloP_score 29way_pi 29way_logOdds LRT_Omega UniSNP_ids 
		1000Gp1_AC 1000Gp1_AF 1000Gp1_AFR_AC 1000Gp1_AFR_AF 1000Gp1_EUR_AC 1000Gp1_EUR_AF
		1000Gp1_AMR_AC 1000Gp1_AMR_AF 1000Gp1_ASN_AC 1000Gp1_ASN_AF 
		ESP_AA_AF ESP_EA_AF
	} \t]\n
	exec cat {*}[lsort -dict [glob dbNSFP2.1_variant.chr*]] | grep -v ^# >> pre_var_dbnsfp.tsv
}

job var_hg19_dbnsfp -deps {${dest}/tmp/hg19/pre_var_hg19_dbnsfp} -targets {extra/var_hg19_dbnsfp.tsv extra/var_hg19_dbnsfp.tsv.opt} -vars {dest} -code {
	cd ${dest}/tmp/hg19
	cg select -s {chromosome pos} -f {
		chromosome {begin=$pos - 1} {end=$pos} {type="snp"} ref alt 
		aaref aaalt hg18_pos genename Uniprot_acc Uniprot_id Uniprot_aapos 
		Interpro_domain cds_strand refcodon SLR_test_statistic codonpos fold-degenerate Ancestral_allele 
		Ensembl_geneid Ensembl_transcriptid aapos 
		SIFT_score SIFT_pred 
		Polyphen2_HDIV_score Polyphen2_HDIV_pred Polyphen2_HVAR_score Polyphen2_HVAR_pred 
		LRT_score LRT_pred 
		MutationTaster_score MutationTaster_pred 
		FATHMM_score FATHMM_pred
		GERP_NR GERP_RS PhyloP_score 29way_pi 29way_logOdds LRT_Omega UniSNP_ids 
		1000Gp1_AC 1000Gp1_AF 1000Gp1_AFR_AC 1000Gp1_AFR_AF 1000Gp1_EUR_AC 1000Gp1_EUR_AF 
		1000Gp1_AMR_AC 1000Gp1_AMR_AF 1000Gp1_ASN_AC 1000Gp1_ASN_AF 
		ESP_AA_AF ESP_EA_AF
	} pre_var_dbnsfp.tsv pre2_var_hg19_dbnsfp.tsv
	cg groupby {chromosome begin end type} pre2_var_hg19_dbnsfp.tsv pre3_var_hg19_dbnsfp.tsv
	cg select -f {
		chromosome begin end type {ref=lindex($ref,0)} alt 
		SIFT_score  SIFT_pred
		Polyphen2_HDIV_score Polyphen2_HDIV_pred Polyphen2_HVAR_score Polyphen2_HVAR_pred 
		LRT_score LRT_pred 
		MutationTaster_score MutationTaster_pred 
		FATHMM_score FATHMM_pred
		GERP_NR GERP_RS PhyloP_score 29way_pi 29way_logOdds LRT_Omega 
		ESP_AA_AF ESP_EA_AF
	} pre3_var_hg19_dbnsfp.tsv var_hg19_dbnsfp.tsv.temp
	# move dbNSFPzip files to target
	file copy -force var_hg19_dbnsfp.tsv.temp ${dest}/hg19/extra/var_hg19_dbnsfp.tsv
	file_write ${dest}/hg19/extra/var_hg19_dbnsfp.tsv.opt "fields\t{SIFT_score Polyphen2_HDIV_score Polyphen2_HDIV_pred Polyphen2_HVAR_score Polyphen2_HVAR_pred LRT_score LRT_pred MutationTaster_score MutationTaster_pred FATHMM_score GERP_NR GERP_RS PhyloP_score 29way_pi 29way_logOdds LRT_Omega ESP_AA_AF ESP_EA_AF}"
	file copy dbNSFP2.1.readme.txt ${dest}/hg19/extra/var_hg19_dbnsfp.info
}

# move dbNSFP hg19 files to target hg18
job dbNSFP_hg18_liftover -deps {$dest/hg19/extra/var_hg19_dbnsfp.tsv} -targets {extra/var_${build}_dbnsfp.tsv extra/var_${build}_dbnsfp.tsv.opt extra/var_${build}_dbnsfp.info} -vars {dest build db} -code {
	cg liftover $dest/hg19/extra/var_hg19_dbnsfp.tsv ${dest}/liftover/hg19ToHg18.over.chain ${dest}/tmp/${build}/var_${build}_dbnsfp.tsv.temp
	file rename -force ${dest}/tmp/${build}/var_${build}_dbnsfp.tsv.temp ${dest}/${build}/extra/var_${build}_dbnsfp.tsv
	file copy -force $dest/hg19/extra/var_hg19_dbnsfp.tsv.opt ${dest}/${build}/extra/var_${build}_dbnsfp.tsv.opt
	file copy -force $dest/hg19/extra/var_hg19_dbnsfp.info ${dest}/${build}/extra/var_${build}_dbnsfp.info
}

job_wait

# todo
# cd /complgen/backup2013-10/hg18/
# cp reg_hg18_conserved.tsv reg_hg18_consnoncoding.tsv reg_hg18_refgene.tsv /complgen/hg18/extra
