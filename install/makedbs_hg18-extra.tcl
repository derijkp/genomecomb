#!/bin/sh
# the next line restarts using tclsh \
exec cgsh "$0" "$@"

set build hg18
set dest /complgen/refseq/

file mkdir ${dest}/${build}/extra
cd ${dest}/${build}/extra

# dbNSFPzip
if {![file exists ${dest}/tmp/hg19/pre_var_hg19_dbnsfp.tsv]} {
	mkdir ${dest}/tmp/hg19
	cd ${dest}/tmp/hg19
	exec wget http://dbnsfp.houstonbioinformatics.org/dbNSFPzip/dbNSFP2.0b3.zip
	exec wget http://dbnsfp.houstonbioinformatics.org/dbNSFPzip/dbNSFP2.0b3.readme.txt
	exec unzip dbNSFP2.0b3.zip
	file_write pre_var_dbnsfp.tsv "chromosome\tpos\tref\talt\taaref\taaalt\thg18_pos\tgenename\tUniprot_acc\tUniprot_id\tUniprot_aapos\tInterpro_domain\tcds_strand\trefcodon\tSLR_test_statistic\tcodonpos\tfold-degenerate\tAncestral_allele\tEnsembl_geneid\tEnsembl_transcriptid\taapos\tSIFT_score\tPolyphen2_HDIV_score\tPolyphen2_HDIV_pred\tPolyphen2_HVAR_score\tPolyphen2_HVAR_pred\tLRT_score\tLRT_pred\tMutationTaster_score\tMutationTaster_pred\tGERP_NR\tGERP_RS\tPhyloP_score\t29way_pi\t29way_logOdds\tLRT_Omega\tUniSNP_ids\t1000Gp1_AC\t1000Gp1_AF\t1000Gp1_AFR_AC\t1000Gp1_AFR_AF\t1000Gp1_EUR_AC\t1000Gp1_EUR_AF\t1000Gp1_AMR_AC\t1000Gp1_AMR_AF\t1000Gp1_ASN_AC\t1000Gp1_ASN_AF\tESP5400_AA_AF\tESP5400_EA_AF"
	exec cat dbNSFP2.0b3_variant.chr* | grep -v ^# >> pre_var_dbnsfp.tsv
	cd ${dest}/${build}/extra
}
cg select -s {chromosome hg18_pos} -f {chromosome {begin=$hg18_pos - 1} {end=$hg18_pos} {type="snp"} ref alt SIFT_score Polyphen2_HDIV_score Polyphen2_HDIV_pred Polyphen2_HVAR_score Polyphen2_HVAR_pred LRT_score LRT_pred MutationTaster_score MutationTaster_pred GERP_NR GERP_RS PhyloP_score 29way_pi 29way_logOdds LRT_Omega ESP5400_AA_AF ESP5400_EA_AF} pre_var_dbnsfp.tsv pre2_var_${build}_dbnsfp.tsv
cg groupby {chromosome begin end type} pre2_var_${build}_dbnsfp.tsv pre3_var_${build}_dbnsfp.tsv
cg select -f {chromosome begin end type {ref=lindex($ref,0)} alt SIFT_score Polyphen2_HDIV_score Polyphen2_HDIV_pred Polyphen2_HVAR_score Polyphen2_HVAR_pred LRT_score LRT_pred MutationTaster_score MutationTaster_pred GERP_NR GERP_RS PhyloP_score 29way_pi 29way_logOdds LRT_Omega ESP5400_AA_AF ESP5400_EA_AF} pre3_var_${build}_dbnsfp.tsv var_${build}_dbnsfp.tsv.temp
file rename var_${build}_dbnsfp.tsv.temp ${dest}/${build}/extra/var_${build}_dbnsfp.tsv
file_write ${dest}/${build}/extra/var_${build}_dbnsfp.tsv.opt "fields\t{SIFT_score Polyphen2_HDIV_score Polyphen2_HDIV_pred Polyphen2_HVAR_score Polyphen2_HVAR_pred LRT_score LRT_pred MutationTaster_score MutationTaster_pred GERP_NR GERP_RS PhyloP_score 29way_pi 29way_logOdds LRT_Omega ESP5400_AA_AF ESP5400_EA_AF}"

cd ${dest}/${build}/extra
