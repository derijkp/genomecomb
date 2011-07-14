# region databases (ucsc)
cd /complgen/refseq/hg18

# download genome
cg downloadgenome hg18 genome_hg18.ifas
cg make_genomecindex genome_hg18.ifas
mkdir extra
mv reg_genome_hg18.tsv extra/reg_hg18_fullgenome.tsv


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
cg maketabix /complgen/refseq/hg18/var_hg18_snp130.tsv
gunzip -c /complgen/refseq/hg18/var_hg18_snp130.tsv.gz > /complgen/refseq/hg18/var_hg18_snp130.tsv

# genes
cg downloaddb /complgen/refseq/hg18 hg18 refGene ensGene knownGene genscan
mv ucsc_hg18_refGene.tsv gene_hg18_refGene.tsv
cg maketabix gene_hg18_refGene.tsv
mv ucsc_hg18_ensGene.tsv gene_hg18_ensGene.tsv
cg maketabix gene_hg18_knownGene.tsv
mv ucsc_hg18_knownGene.tsv gene_hg18_knownGene.tsv
cg maketabix gene_hg18_ensGene.tsv
echo -e "genecol\tproteinID\ntranscriptcol\tname" > gene_hg18_knownGene.tsv.opt
# todo: also get kgXref, and translate
cg downloaddb /complgen/refseq/hg18 hg18 genscan
mv ucsc_hg18_genscan.tsv gene_hg18_genscan.tsv
cg maketabix gene_hg18_genscan.tsv


# homopolymer
cd /complgen/refseq/hg18
cg extracthomopolymers genome_hg18.ifas > reg_hg18_homopolymer.tsv.temp
mv reg_hg18_homopolymer.tsv.temp reg_hg18_homopolymer.tsv
echo -e "fields\t{base size}" > reg_hg18_homopolymer.tsv.opt
