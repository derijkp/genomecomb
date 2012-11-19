#!/bin/sh
# the next line restarts using tclsh \
exec cgsh "$0" "$@"

set build hg19
set dest /complgen/refseq/

file mkdir ${dest}/${build}/extra
cd ${dest}/${build}/extra

# targetseq exome
wget http://tools.invitrogen.com/content/sfs/manuals/TargetSeq_exome_named_targets_hg19.bed
cg bed2sft TargetSeq_exome_named_targets_hg19.bed ureg_exome_targetseq.tsv
cg select -s 'chromosome begin end' ureg_exome_targetseq.tsv sreg_exome_targetseq.tsv
cg collapseoverlap -o reg_exome_targetseq.tsv sreg_exome_targetseq.tsv
rm TargetSeq_exome_named_targets_hg19.bed sreg_exome_targetseq.tsv ureg_exome_targetseq.tsv
# cg select -f '* {id=NR "-" $name}' reg_exome_targetseq.tsv | less

# dbNSFPzip
cd ${dest}/tmp/hg19
if {![file exists ${dest}/tmp/hg19/pre_var_hg19_dbnsfp.tsv]} {
	if {![file exists dbNSFP2.0b4.zip]} {
		exec wget -c http://dbnsfp.houstonbioinformatics.org/dbNSFPzip/dbNSFP2.0b4.zip
	}
	exec unzip dbNSFP2.0b4.zip >@ stdout 2>@ stderr
	file_write pre_var_dbnsfp.tsv "chromosome\tpos\tref\talt\taaref\taaalt\thg18_pos\tgenename\tUniprot_acc\tUniprot_id\tUniprot_aapos\tInterpro_domain\tcds_strand\trefcodon\tSLR_test_statistic\tcodonpos\tfold-degenerate\tAncestral_allele\tEnsembl_geneid\tEnsembl_transcriptid\taapos\tSIFT_score\tPolyphen2_HDIV_score\tPolyphen2_HDIV_pred\tPolyphen2_HVAR_score\tPolyphen2_HVAR_pred\tLRT_score\tLRT_pred\tMutationTaster_score\tMutationTaster_pred\tGERP_NR\tGERP_RS\tPhyloP_score\t29way_pi\t29way_logOdds\tLRT_Omega\tUniSNP_ids\t1000Gp1_AC\t1000Gp1_AF\t1000Gp1_AFR_AC\t1000Gp1_AFR_AF\t1000Gp1_EUR_AC\t1000Gp1_EUR_AF\t1000Gp1_AMR_AC\t1000Gp1_AMR_AF\t1000Gp1_ASN_AC\t1000Gp1_ASN_AF\tESP_AA_AF\tESP_EA_AF\n"
	exec cat {*}[glob dbNSFP2.0b4_variant.chr*] | grep -v ^# >> pre_var_dbnsfp.tsv
}
cg select -s {chromosome pos} -f {chromosome {begin=$pos - 1} {end=$pos} {type="snp"} ref alt aaref aaalt hg18_pos genename Uniprot_acc Uniprot_id Uniprot_aapos Interpro_domain cds_strand refcodon SLR_test_statistic codonpos fold-degenerate Ancestral_allele Ensembl_geneid Ensembl_transcriptid aapos SIFT_score Polyphen2_HDIV_score Polyphen2_HDIV_pred Polyphen2_HVAR_score Polyphen2_HVAR_pred LRT_score LRT_pred MutationTaster_score MutationTaster_pred GERP_NR GERP_RS PhyloP_score 29way_pi 29way_logOdds LRT_Omega UniSNP_ids 1000Gp1_AC 1000Gp1_AF 1000Gp1_AFR_AC 1000Gp1_AFR_AF 1000Gp1_EUR_AC 1000Gp1_EUR_AF 1000Gp1_AMR_AC 1000Gp1_AMR_AF 1000Gp1_ASN_AC 1000Gp1_ASN_AF ESP_AA_AF ESP_EA_AF} pre_var_dbnsfp.tsv pre2_var_hg19_dbnsfp.tsv
cg groupby {chromosome begin end type} pre2_var_hg19_dbnsfp.tsv pre3_var_hg19_dbnsfp.tsv
cg select -f {chromosome begin end type {ref=lindex($ref,0)} alt SIFT_score Polyphen2_HDIV_score Polyphen2_HDIV_pred Polyphen2_HVAR_score Polyphen2_HVAR_pred LRT_score LRT_pred MutationTaster_score MutationTaster_pred GERP_NR GERP_RS PhyloP_score 29way_pi 29way_logOdds LRT_Omega ESP_AA_AF ESP_EA_AF} pre3_var_hg19_dbnsfp.tsv var_hg19_dbnsfp.tsv.temp
# move dbNSFPzip files to target
file copy -force var_hg19_dbnsfp.tsv.temp ${dest}/${build}/extra/var_${build}_dbnsfp.tsv
file_write ${dest}/${build}/extra/var_${build}_dbnsfp.tsv.opt "fields\t{SIFT_score Polyphen2_HDIV_score Polyphen2_HDIV_pred Polyphen2_HVAR_score Polyphen2_HVAR_pred LRT_score LRT_pred MutationTaster_score MutationTaster_pred GERP_NR GERP_RS PhyloP_score 29way_pi 29way_logOdds LRT_Omega ESP_AA_AF ESP_EA_AF}"
file copy dbNSFP2.0b4.readme.txt ${dest}/${build}/extra/var_${build}_dbnsfp.info
