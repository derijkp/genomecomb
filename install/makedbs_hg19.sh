#!/bin/sh
# the next line restarts using tclsh \
exec cg source "$0" "$@"

set build hg19
set argv [job_init {*}$argv]
if {[llength $argv]} {
	set dest [lindex $argv 0]
} else {
	set dest /complgen/refseq/
}


# download hg19
# =============
#
file mkdir ${dest}/${build}
cd ${dest}/${build}
file mkdir extra

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
	cytoBand evofold microsat oreganno rmsk simpleRepeat targetScanS tfbsConsSites 
	tRNAs wgRna vistaEnhancers gad
	phastConsElements46way phastConsElements46wayPlacental phastConsElements46wayPrimates
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

job reg_${build}_gwasCatalog -vars {build dest} -targets {reg_${build}_gwasCatalog.tsv} -code {
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
# 1000g (hg18)
job 1000gh18 -targets {$dest/hg18/var_hg18_1000gCHBxJPT.tsv $dest/hg18/var_hg18_1000gCEU.tsv $dest/hg18/var_hg18_1000gYRI.tsv} -vars {dest} -code {
	cd ${dest}/hg18
	cg downloaddb ${dest} hg18 1000g
}

job 1000gliftover -deps {$dest/hg18/var_hg18_1000gCHBxJPT.tsv $dest/hg18/var_hg18_1000gCEU.tsv $dest/hg18/var_hg18_1000gYRI.tsv} -targets {var_${build}_1000gCHBxJPT.tsv var_${build}_1000gCEU.tsv var_${build}_1000gYRI.tsv} -vars {dest build} -code {
	cg liftover ${dest}/hg18/var_hg18_1000gCHBxJPT.tsv ${dest}/liftover/hg18ToHg19.over.chain var_${build}_1000gCHBxJPT.tsv
	cg liftover ${dest}/hg18/var_hg18_1000gCEU.tsv ${dest}/liftover/hg18ToHg19.over.chain var_${build}_1000gCEU.tsv
	cg liftover ${dest}/hg18/var_hg18_1000gYRI.tsv ${dest}/liftover/hg18ToHg19.over.chain var_${build}_1000gYRI.tsv
	file delete var_hg19_1000gCEU.tsv.unmapped var_hg19_1000gCHBxJPT.tsv.unmapped var_hg19_1000gYRI.tsv.unmapped
}

job 1000glow -targets {${dest}/hg19/var_hg19_1000glow.tsv} -vars {dest} -code {
	cd $dest/hg19
	cg downloaddb ${dest} hg19 1000glow
}

# dbsnp
job dbsnp138 -targets {${dest}/hg19/var_hg19_snp138.tsv} -vars {dest} -code {
	cd $dest/hg19
	cg downloaddb ${dest} hg19 snp138
}

job dbsnp138Common -targets {${dest}/hg19/var_hg19_snp138Common.tsv} -vars {dest build} -code {
	cd $dest/hg19
	cg downloaddb ${dest} hg19 snp138Common
}

foreach db {
	snp138 snp138Common
} {
	job maketabix_${build}_$db -deps {var_${build}_${db}.tsv} -targets {var_${build}_${db}.tsv.gz.tbi var_${build}_${db}.tsv.gz} -vars {dest build db} -code {
		cg maketabix var_${build}_${db}.tsv
		# exec gunzip -c var_${build}_${db}.tsv.gz > var_${build}_${db}.tsv
		cg select -f {chrom start end type ref alt name freq} ${dest}/${build}/var_${build}_${db}.tsv.gz ${dest}/${build}/var_${build}_${db}.tsv
	}
}

# genes
foreach db {
	refGene ensGene knownGene genscan acembly
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


# gencode
job gene_${build}_gencode -targets {gene_${build}_gencode.tsv gene_${build}_gencode.tsv.gz gene_${build}_gencode.tsv.gz.tbi gene_${build}_gencode.info} -vars {dest build db} -code {
	exec -ignorestderr wget -c ftp://ftp.sanger.ac.uk/pub/gencode/release_18/gencode.v18.annotation.gtf.gz
	cg gtf2sft gencode.v18.annotation.gtf.gz gene_${build}_gencode.tsv.temp
	cg select -s - gene_${build}_gencode.tsv.temp gene_${build}_gencode.tsv.temp2
	file delete gene_${build}_gencode.tsv.temp
	file rename -force gene_${build}_gencode.tsv.temp2 gene_${build}_gencode.tsv
	file delete gene_${build}_gencode.tsv.gz
	cg maketabix gene_${build}_gencode.tsv
	exec gunzip -c gene_${build}_gencode.tsv.gz > gene_${build}_gencode.tsv
	exec -ignorestderr wget -c ftp://ftp.sanger.ac.uk/pub/gencode/_README.TXT
	file rename -force _README.TXT gene_${build}_gencode.info
	file delete gencode.v18.annotation.gtf.gz
}

# homopolymer
job reg_${build}_homopolymer -deps {genome_${build}.ifas} -targets {reg_${build}_homopolymer.tsv reg_${build}_homopolymer.tsv.gz reg_${build}_homopolymer.tsv.gz.tbi reg_${build}_homopolymer.tsv.opt} -vars {dest build db} -code {
	cg extracthomopolymers genome_${build}.ifas > reg_${build}_homopolymer.tsv.temp
	file rename -force reg_${build}_homopolymer.tsv.temp reg_${build}_homopolymer.tsv
        cg maketabix reg_${build}_homopolymer.tsv
        gunzip reg_${build}_homopolymer.tsv.gz reg_${build}_homopolymer.tsv
	file_write reg_${build}_homopolymer.tsv.opt "fields\t{base size}\n"
}

# mirbase
job reg_hg19_mirbase -targets {$dest/hg19/reg_hg19_mirbase.tsv $dest/hg19/reg_hg19_mirbase.info} -vars {dest build db} -code {
	cd $dest/hg19
	exec -ignorestderr wget -c --tries=45 --directory-prefix=${dest}/tmp/hg19 ftp://mirbase.org/pub/mirbase/20/genomes/hsa.gff2
	cg gff2sft ${dest}/tmp/hg19/hsa.gff2 ${dest}/tmp/hg19/reg_hg19_mirbase.tsv.temp
	cg select -s - ${dest}/tmp/hg19/reg_hg19_mirbase.tsv.temp ${dest}/tmp/hg19/reg_hg19_mirbase.tsv.temp2
	file rename -force ${dest}/tmp/hg19/reg_hg19_mirbase.tsv.temp2 reg_hg19_mirbase.tsv
	exec -ignorestderr wget -c ftp://mirbase.org/pub/mirbase/20/README
	file rename -force README reg_hg19_mirbase.info
}

# GERP
job GERP -targets {extra/reg_${build}_GERP.tsv extra/reg_${build}_GERP.info} -vars {dest build tables} -code {
	cd ${dest}/tmp/${build}
	set table allHg19RS_BW
	cg downloaddbinfo ${dest}/tmp ${build} $table
	file rename -force reg_${build}_$table.info ${dest}/extra/${build}/reg_${build}_GERP.info
	cg downloaddb ${dest}/tmp ${build} $table
	cg ucscwb2reg -p 1 -f {} ucsc_${build}_$table.tsv
	cg select -s - reg_ucsc_${build}_$table.tsv ${dest}/${build}/extra/reg_${build}_GERP.tsv.temp
	mv ${dest}/${build}/extra/reg_${build}_GERP.tsv.temp ${dest}/${build}/extra/reg_${build}_GERP.tsv
#	file rename -force reg_ucsc_${build}_$table.tsv ${dest}/${build}/extra/reg_${build}_GERP.tsv

}

# encode
foreach {jobname resultname infosrc tables} {
	enc_transcription wgEncodeCaltechRnaSeq wgEncodeCaltechRnaSeq {wgEncodeRegTxnCaltechRnaSeqGm12878R2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqH1hescR2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqHelas3R2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqHepg2R2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqHsmmR2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqHuvecR2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqK562R2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqNhekR2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqNhlfR2x75Il200SigPooled}
	enc_H3K4Me1 wgEncodeH3k4me1 wgEncodeRegMarkEnhH3k4me1 {wgEncodeBroadHistoneGm12878H3k4me1StdSig wgEncodeBroadHistoneH1hescH3k4me1StdSig wgEncodeBroadHistoneHsmmH3k4me1StdSig wgEncodeBroadHistoneHuvecH3k4me1StdSig wgEncodeBroadHistoneK562H3k4me1StdSig wgEncodeBroadHistoneNhekH3k4me1StdSig wgEncodeBroadHistoneNhlfH3k4me1StdSig}
	enc_H3K4Me3	wgEncodeH3k4me3 wgEncodeRegMarkPromoter   {wgEncodeBroadHistoneGm12878H3k4me3StdSig wgEncodeBroadHistoneH1hescH3k4me3StdSig wgEncodeBroadHistoneHsmmH3k4me3StdSig wgEncodeBroadHistoneHuvecH3k4me3StdSig wgEncodeBroadHistoneK562H3k4me3StdSig wgEncodeBroadHistoneNhekH3k4me3StdSig wgEncodeBroadHistoneNhlfH3k4me3StdSig}
	enc_H3K27Ac wgEncodeH3k27ac wgEncodeRegMarkEnhH3k27ac {wgEncodeBroadHistoneGm12878H3k27acStdSig wgEncodeBroadHistoneH1hescH3k27acStdSig wgEncodeBroadHistoneHsmmH3k27acStdSig wgEncodeBroadHistoneHuvecH3k27acStdSig wgEncodeBroadHistoneK562H3k27acStdSig wgEncodeBroadHistoneNhekH3k27acStdSig wgEncodeBroadHistoneNhlfH3k27acStdSig}
	
} {
	# make database
	job $jobname -targets {${dest}/${build}/reg_${build}_$resultname.tsv} -vars {dest build tables} -code {
		cd ${dest}/tmp/${build}
		cg downloaddb ${dest}/tmp ${build} {*}$tables
		set todo {}
		foreach table $tables {
			cg ucscwb2reg -n 10 -p 0 -f {5*(($value+4)/5)} ucsc_${build}_$table.tsv
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
job enc_RegDnaseClustered -targets {reg_${build}_wgEncodeRegDnaseClusteredV2.tsv} -vars {dest build} -code {
	cd ${dest}/tmp/${build}
	cg downloaddb ${dest}/tmp ${build} wgEncodeRegDnaseClusteredV2
	cg collapseoverlap ucsc_${build}_wgEncodeRegDnaseClusteredV2.tsv
	file rename -force reg_${build}_wgEncodeRegDnaseClusteredV2.tsv ${dest}/${build}
	cg downloaddbinfo ${dest}/tmp ${build} wgEncodeRegDnaseClusteredV2
	file rename -force reg_${build}_wgEncodeRegDnaseClusteredV2.info ${dest}/${build}/reg_${build}_wgEncodeRegDnaseClusteredV2.info
}

job enc_RegTfbsClustered -targets {reg_${build}_wgEncodeRegTfbsClusteredV3.tsv} -vars {dest build} -code {
	cd ${dest}/tmp/${build}
	cg downloaddb ${dest}/tmp ${build} wgEncodeRegTfbsClusteredV3
	cg select -s - -f {chrom	start	end	name	score} ucsc_${build}_wgEncodeRegTfbsClusteredV3.tsv pucsc_${build}_wgEncodeRegTfbsClusteredV3.tsv
	cg collapseoverlap pucsc_${build}_wgEncodeRegTfbsClusteredV3.tsv
	file rename -force reg_pucsc_${build}_wgEncodeRegTfbsClusteredV3.tsv ${dest}/${build}/reg_${build}_wgEncodeRegTfbsClusteredV3.tsv
	cg downloaddbinfo ${dest}/tmp ${build} wgEncodeRegTfbsClusteredV3
	file rename -force reg_${build}_wgEncodeRegTfbsClusteredV3.info ${dest}/${build}/reg_${build}_wgEncodeRegTfbsClusteredV3.info
}

# link local data in dir
catch {
	foreach file [glob ../hg19-local/*] {
		file delete extra/[file tail $file]
		cplinked $file [file tail $file]
	}
}

# extra dir
# targetseq exome
job reg_exome_targetseq -targets {extra/reg_hg19_exome_targetseq.tsv} -code {
	cd extra
	exec -ignorestderr wget -c http://tools.invitrogen.com/content/sfs/manuals/TargetSeq_exome_named_targets_hg19.bed
	cg bed2sft TargetSeq_exome_named_targets_hg19.bed ureg_hg19_exome_targetseq.tsv
	cg select -s {chromosome begin end} ureg_hg19_exome_targetseq.tsv sreg_hg19_exome_targetseq.tsv
	cg collapseoverlap -o reg_hg19_exome_targetseq.tsv sreg_hg19_exome_targetseq.tsv
	file delete TargetSeq_exome_named_targets_hg19.bed sreg_hg19_exome_targetseq.tsv ureg_hg19_exome_targetseq.tsv
}
# cg select -f '* {id=NR "-" $name}' reg_hg19_exome_targetseq.tsv | less

# dbNSFPzip
file mkdir ${dest}/tmp/hg19
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

job extragenome -deps {genome_${build}.ifas genome_${build}.ifas.index genome_${build}.ssa} -vars build \
-targets {extra/genome_${build}.ifas extra/genome_${build}.ifas.index extra/genome_${build}.ssa} -code {
	exec ln -s genome_${build}.ifas extra/genome_${build}.ifas
	exec ln -s genome_${build}.ifas.index extra/genome_${build}.ifas.index 
	exec ln -s genome_${build}.ssa extra/genome_${build}.ssa
}

# genome in extra
catch {
	foreach file [glob genome_*] {
		file delete extra/[file tail $file]
		mklink $file extra/[file tail $file]
	}
}

job_wait

# todo
# cd /complgen/refseq/backup/hg19/extra
# cp reg_hg19_conserved.tsv reg_hg19_consnoncoding.tsv reg_hg19_refgene.tsv /complgen/refseq/hg19/extra
# cp reg_hg19_exome_SureSelectV4.tsv reg_hg19_exome_targetseq.tsv /complgen/refseq/hg19/extra/
# rs genome_hg19* /complgen/refseq/hg19/extra
