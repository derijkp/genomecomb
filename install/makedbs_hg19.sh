export build=hg19
export dest=/complgen/refseq/

# download hg19
# =============
#
mkdir -p ${dest}/${build}
cd ${dest}/${build}

# download genome
cg downloadgenome ${build} genome_${build}.ifas
cg make_genomecindex genome_${build}.ifas
mkdir extra
mv reg_genome_${build}.tsv extra/reg_${build}_fullgenome.tsv
cg calcsequencedgenome ${dest} ${build}

# region databases (ucsc)
cg downloaddb ${dest} ${build} simpleRepeat microsat rmsk genomicSuperDups chainSelf
cg downloaddb ${dest} ${build} omimGene oreganno tRNAs targetScanS evofold refLink
cg downloaddb ${dest} ${build} cytoBand dgv gwasCatalog kgXref wgRna vistaEnhancers gad tfbsConsSites
# cg downloaddb ${dest} ${build} phastConsElements28way phastConsElements28wayPlacMammal phastConsElements44way rnaGene
mv ucsc_${build}_gwasCatalog.tsv ucsc_${build}_gwasCatalog.tsv.temp
cg select \
	-f 'chrom start end trait pValue pubMedID name bin author pubDate journal title initSample replSample region genes riskAllele riskAlFreq pValueDesc orOrBeta ci95 platform cnv' \
	-nh 'chrom start end name score pubMedID dbsnp bin author pubDate journal title initSample replSample region genes riskAllele riskAlFreq pValueDesc orOrBeta ci95 platform cnv' \
	ucsc_${build}_gwasCatalog.tsv.temp ucsc_${build}_gwasCatalog.tsv

# collapse regions
cg collapseoverlap ucsc_${build}_cytoBand.tsv ucsc_${build}_evofold.tsv ucsc_${build}_gwasCatalog.tsv ucsc_${build}_microsat.tsv ucsc_${build}_omimGene.tsv ucsc_${build}_oreganno.tsv ucsc_${build}_rmsk.tsv ucsc_${build}_simpleRepeat.tsv ucsc_${build}_targetScanS.tsv ucsc_${build}_tfbsConsSites.tsv ucsc_${build}_tRNAs.tsv ucsc_${build}_wgRna.tsv ucsc_${build}_vistaEnhancers.tsv ucsc_${build}_gad.tsv
# cg collapseoverlap ucsc_${build}_phastConsElements28wayPlacMammal.tsv ucsc_${build}_phastConsElements28way.tsv ucsc_${build}_phastConsElements44way.tsv ucsc_${build}_rnaGene.tsv
cg regjoin ucsc_${build}_chainSelf.tsv > reg_${build}_chainSelf.tsv
cg regjoin ucsc_${build}_dgv.tsv > reg_${build}_dgv.tsv
cg regjoin ucsc_${build}_genomicSuperDups.tsv > reg_${build}_genomicSuperDups.tsv
# shorten some names
#mv ucsc_${build}_phastConsElements28wayPlacMammal.tsv ucsc_${build}_phastCons28P.tsv
#mv  reg_${build}_phastConsElements28wayPlacMammal.tsv  reg_${build}_phastCons28P.tsv
#mv ucsc_${build}_phastConsElements28way.tsv ucsc_${build}_phastCons28.tsv
#mv  reg_${build}_phastConsElements28way.tsv  reg_${build}_phastCons28.tsv
#mv ucsc_${build}_phastConsElements44way.tsv ucsc_${build}_phastCons44.tsv
#mv  reg_${build}_phastConsElements44way.tsv  reg_${build}_phastCons44.tsv

# 1000 genomes
cg downloaddb ${dest}/hg18 hg18 1000g
cg liftover ${dest}/hg18/var_hg18_1000gCHBxJPT.tsv ${dest}/liftover/hg18ToHg19.over.chain ${dest}/${build}/var_${build}_1000gCHBxJPT.tsv
cg liftover ${dest}/hg18/var_hg18_1000gCEU.tsv ${dest}/liftover/hg18ToHg19.over.chain ${dest}/${build}/var_${build}_1000gCEU.tsv
cg liftover ${dest}/hg18/var_hg18_1000gYRI.tsv ${dest}/liftover/hg18ToHg19.over.chain ${dest}/${build}/var_${build}_1000gYRI.tsv
rm var_hg19_1000gCEU.tsv.unmapped var_hg19_1000gCHBxJPT.tsv.unmapped var_hg19_1000gYRI.tsv.unmapped

cg downloaddb ${dest} ${build} 1000glow

# dbsnp
cg downloaddb ${dest} ${build} snp132
cg bgzip ${dest}/${build}/var_${build}_snp132.tsv
cg maketabix ${dest}/${build}/var_${build}_snp132.tsv.gz
gunzip -c ${dest}/${build}/var_${build}_snp132.tsv.gz > ${dest}/${build}/var_${build}_snp132.tsv

# genes
cg downloaddb ${dest} ${build} refGene ensGene knownGene genscan
mv ucsc_${build}_refGene.tsv gene_${build}_refGene.tsv
cg maketabix gene_${build}_refGene.tsv
mv ucsc_${build}_ensGene.tsv gene_${build}_ensGene.tsv
cg maketabix gene_${build}_ensGene.tsv
echo -e "genecol\tproteinID\ntranscriptcol\tname" > gene_${build}_knownGene.tsv.opt
mv ucsc_${build}_knownGene.tsv gene_${build}_knownGene.tsv
cg maketabix gene_${build}_knownGene.tsv
# todo: also get kgXref, and translate
cg downloaddb ${dest} ${build} genscan
mv ucsc_${build}_genscan.tsv gene_${build}_genscan.tsv
cg maketabix gene_${build}_genscan.tsv

# homopolymer
cd ${dest}/${build}
cg extracthomopolymers genome_${build}.ifas > reg_${build}_homopolymer.tsv.temp
mv reg_${build}_homopolymer.tsv.temp reg_${build}_homopolymer.tsv
echo -e "fields\t{base size}" > reg_${build}_homopolymer.tsv.opt
