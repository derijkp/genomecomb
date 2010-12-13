# region databases (ucsc)
cd /complgen/refseq/hg18
cg downloaddb /complgen/refseq/hg18 hg18 simpleRepeat microsat rmsk genomicSuperDups chainSelf
cg downloaddb /complgen/refseq/hg18 hg18 omimGene phastConsElements44way oreganno rnaGene tRNAs tfbsConsSites targetScanS evofold
cg downloaddb /complgen/refseq/hg18 hg18 cytoBand dgv gwasCatalog kgXref phastConsElements28wayPlacMammal phastConsElements28way wgRna refLink vistaEnhancers gad
mv ucsc_hg18_gwasCatalog.tsv ucsc_hg18_gwasCatalog.tsv.temp
cg select \
	-f 'chrom start end trait pValue pubMedID name bin author pubDate journal title initSample replSample region genes riskAllele riskAlFreq pValueDesc orOrBeta ci95 platform cnv' \
	-nh 'chrom start end name score pubMedID dbsnp bin author pubDate journal title initSample replSample region genes riskAllele riskAlFreq pValueDesc orOrBeta ci95 platform cnv' \
	ucsc_hg18_gwasCatalog.tsv.temp ucsc_hg18_gwasCatalog.tsv

# collapse regions
cg collapseoverlap ucsc_hg18_cytoBand.tsv ucsc_hg18_evofold.tsv ucsc_hg18_gwasCatalog.tsv ucsc_hg18_microsat.tsv ucsc_hg18_omimGene.tsv ucsc_hg18_oreganno.tsv ucsc_hg18_phastConsElements28wayPlacMammal.tsv ucsc_hg18_phastConsElements28way.tsv ucsc_hg18_phastConsElements44way.tsv ucsc_hg18_rmsk.tsv ucsc_hg18_rnaGene.tsv ucsc_hg18_simpleRepeat.tsv ucsc_hg18_targetScanS.tsv ucsc_hg18_tfbsConsSites.tsv ucsc_hg18_tRNAs.tsv ucsc_hg18_wgRna.tsv ucsc_hg18_vistaEnhancers.tsv ucsc_hg18_gad.tsv
cg regjoin ucsc_hg18_chainSelf.tsv > reg_hg18_chainSelf.tsv
cg regjoin ucsc_hg18_dgv.tsv > reg_hg18_dgv.tsv
cg regjoin ucsc_hg18_genomicSuperDups.tsv > reg_hg18_genomicSuperDups.tsv
# shorten some names
mv ucsc_hg18_phastConsElements28wayPlacMammal.tsv ucsc_hg18_phastCons28P.tsv
mv  reg_hg18_phastConsElements28wayPlacMammal.tsv  reg_hg18_phastCons28P.tsv
mv ucsc_hg18_phastConsElements28way.tsv ucsc_hg18_phastCons28.tsv
mv  reg_hg18_phastConsElements28way.tsv  reg_hg18_phastCons28.tsv
mv ucsc_hg18_phastConsElements44way.tsv ucsc_hg18_phastCons44.tsv
mv  reg_hg18_phastConsElements44way.tsv  reg_hg18_phastCons44.tsv

# 1000 genomes
cg downloaddb /complgen/refseq/hg18 hg18 1000g
# dbdnp
cg downloaddb /complgen/refseq/hg18 hg18 snp130

# genes
cg downloaddb /complgen/refseq/hg18 hg18 refGene ensGene knownGene genscan
mv ucsc_hg18_refGene.tsv gene_hg18_refGene.tsv
mv ucsc_hg18_ensGene.tsv gene_hg18_ensGene.tsv
mv ucsc_hg18_knownGene.tsv gene_hg18_knownGene.tsv




## annovar
## -------
#cd /complgen/refseq/annovar
#perl annotate_variation.pl -downdb -buildver hg18 gene /complgen/refseq/hg18
#perl annotate_variation.pl -downdb -buildver hg18 knownGene /complgen/refseq/hg18
#perl annotate_variation.pl -downdb -buildver hg18 ensGene /complgen/refseq/hg18
#perl annotate_variation.pl -downdb -buildver hg18 avsift /complgen/refseq/hg18
#
##!/bin/sh
#cd /complgen/refseq
#cg select -q '$chrom ~ /^chr[0-9XYM][0-9]*$/' -s 'chrom chromStart chromEnd' <  ucsc_ori/_data_db_ucsc-exapted_repeats.tsv  > regdb-exapted_repeats.tsv
#cg regjoin oregdb-exapted_repeats.tsv > regdb-exapted_repeats.tsv
#cg select -q '$chrom ~ /^chr[0-9XYM][0-9]*$/' -s 'chrom chromStart chromEnd' < ucsc_ori/_data_db_ucsc-microsatelite.tsv > regdb-microsatelite.tsv
#cg regjoin oregdb-microsatelite.tsv > regdb-microsatelite.tsv
#cg select -q '$genoName ~ /^chr[0-9XYM][0-9]*$/' -s 'genoName genoStart genoEnd' < ucsc_ori/_data_db_ucsc-repeatmasker.tsv > regdb-repeatmasker.tsv
#cg regjoin oregdb-repeatmasker.tsv > regdb-repeatmasker.tsv
#cg select -q '$chrom ~ /^chr[0-9XYM][0-9]*$/' -s 'chrom chromStart chromEnd' < ucsc_ori/_data_db_ucsc-rnagenes.tsv > regdb-rnagenes.tsv
#cg regjoin oregdb-rnagenes.tsv > regdb-rnagenes.tsv
#cg select -q '$chrom ~ /^chr[0-9XYM][0-9]*$/' -s 'chrom chromStart chromEnd' < ucsc_ori/_data_db_ucsc-segdups.tsv > regdb-segdups.tsv
#cg regjoin oregdb-segdups.tsv > regdb-segdups.tsv
#cg select -q '$chrom ~ /^chr[0-9XYM][0-9]*$/' -s 'chrom chromStart chromEnd' < ucsc_ori/_data_db_ucsc-simple_repeats.tsv > regdb-simple_repeats.tsv
#cg regjoin oregdb-simple_repeats.tsv > regdb-simple_repeats.tsv
#cg select -q '$tName ~ /^chr[0-9XYM][0-9]*$/' -s 'tName tStart tEnd' < ucsc_ori/_data_db_ucsc-selfchain.tsv > temp-selfchain.tsv
#cg regjoin temp-selfchain.tsv > regdb-selfchain.tsv
#cg covered regdb-selfchain.tsv > regdb-selfchain.covered
#
