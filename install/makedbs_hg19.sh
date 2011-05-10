export dest=/complgen/refseq/hg19
# download hg19
# =============
#
mkdir $dest
cd $dest

# download genome
cg downloadgenome hg19 genome_hg19.ifas
cg make_genomecindex genome_hg19.ifas
mkdir extra
mv reg_genome_hg19.tsv extra/reg_hg19_fullgenome.tsv

# region databases (ucsc)
cg downloaddb $dest hg19 simpleRepeat microsat rmsk genomicSuperDups chainSelf
cg downloaddb $dest hg19 omimGene oreganno tRNAs targetScanS evofold
# cg downloaddb phastConsElements44way rnaGene tfbsConsSites kgXref phastConsElements28wayPlacMammal refLink phastConsElements28way
cg downloaddb $dest hg19 cytoBand dgv gwasCatalog wgRna vistaEnhancers gad
mv ucsc_hg19_gwasCatalog.tsv ucsc_hg19_gwasCatalog.tsv.temp
cg select \
	-f 'chrom start end trait pValue pubMedID name bin author pubDate journal title initSample replSample region genes riskAllele riskAlFreq pValueDesc orOrBeta ci95 platform cnv' \
	-nh 'chrom start end name score pubMedID dbsnp bin author pubDate journal title initSample replSample region genes riskAllele riskAlFreq pValueDesc orOrBeta ci95 platform cnv' \
	ucsc_hg19_gwasCatalog.tsv.temp ucsc_hg19_gwasCatalog.tsv

# collapse regions
cg collapseoverlap ucsc_hg19_cytoBand.tsv ucsc_hg19_evofold.tsv ucsc_hg19_gwasCatalog.tsv ucsc_hg19_microsat.tsv ucsc_hg19_omimGene.tsv ucsc_hg19_oreganno.tsv ucsc_hg19_rmsk.tsv ucsc_hg19_simpleRepeat.tsv ucsc_hg19_targetScanS.tsv ucsc_hg19_tRNAs.tsv ucsc_hg19_wgRna.tsv ucsc_hg19_vistaEnhancers.tsv ucsc_hg19_gad.tsv
# cg collapseoverlap ucsc_hg19_tfbsConsSites.tsv ucsc_hg19_phastConsElements28wayPlacMammal.tsv ucsc_hg19_phastConsElements44way.tsv ucsc_hg19_rnaGene.tsv ucsc_hg19_phastConsElements28way.tsv
cg regjoin ucsc_hg19_chainSelf.tsv > reg_hg19_chainSelf.tsv
cg regjoin ucsc_hg19_dgv.tsv > reg_hg19_dgv.tsv
cg regjoin ucsc_hg19_genomicSuperDups.tsv > reg_hg19_genomicSuperDups.tsv
# shorten some names
#mv ucsc_hg19_phastConsElements28wayPlacMammal.tsv ucsc_hg19_phastCons28P.tsv
#mv  reg_hg19_phastConsElements28wayPlacMammal.tsv  reg_hg19_phastCons28P.tsv
#mv ucsc_hg19_phastConsElements28way.tsv ucsc_hg19_phastCons28.tsv
#mv  reg_hg19_phastConsElements28way.tsv  reg_hg19_phastCons28.tsv
#mv ucsc_hg19_phastConsElements44way.tsv ucsc_hg19_phastCons44.tsv
#mv  reg_hg19_phastConsElements44way.tsv  reg_hg19_phastCons44.tsv

# 1000 genomes
#cg downloaddb $dest hg19 1000g
# dbsnp
cg downloaddb $dest hg19 snp132

# genes
cg downloaddb $dest hg19 refGene ensGene knownGene genscan
mv ucsc_hg19_refGene.tsv gene_hg19_refGene.tsv
mv ucsc_hg19_ensGene.tsv gene_hg19_ensGene.tsv
mv ucsc_hg19_knownGene.tsv gene_hg19_knownGene.tsv

# homopolymer
cd /complgen/refseq/hg19
cg extracthomopolymers genome_hg19.ifas > reg_hg19_homopolymer.tsv.temp
mv reg_hg19_homopolymer.tsv.temp reg_hg19_homopolymer.tsv
