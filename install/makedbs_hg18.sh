export dest=/complgen/refseq
export build=hg18

# download hg18
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
cg downloaddb ${dest} ${build} omimGene oreganno tRNAs targetScanS evofold
cg downloaddb ${dest} ${build} cytoBand dgv gwasCatalog wgRna vistaEnhancers gad tfbsConsSites
cg downloaddb ${dest} ${build} phastConsElements28way phastConsElements28wayPlacMammal phastConsElements44way rnaGene 
# you can download info on the databases using:
# cg downloaddbinfo ${dest} ${build} simpleRepeat microsat rmsk genomicSuperDups chainSelf
mv ucsc_${build}_gwasCatalog.tsv ucsc_${build}_gwasCatalog.tsv.temp

cg select \
	-f 'chrom start end trait pValue pubMedID name bin author pubDate journal title initSample replSample region genes riskAllele riskAlFreq pValueDesc orOrBeta ci95 platform cnv' \
	-nh 'chrom start end name score pubMedID dbsnp bin author pubDate journal title initSample replSample region genes riskAllele riskAlFreq pValueDesc orOrBeta ci95 platform cnv' \
	ucsc_${build}_gwasCatalog.tsv.temp ucsc_${build}_gwasCatalog.tsv

rm reg_${build}_omimGene.info
cg downloaddbinfo ${dest} ${build} omimGene2
mv reg_${build}_omimGene2.info reg_${build}_omimGene.info

# other databases
cg downloaddb ${dest} ${build} refLink kgXref
mv ucsc_${build}_kgXref.tsv other_${build}_kgXref.tsv
mv ucsc_${build}_refLink.tsv other_${build}_refLink.tsv

# collapse regions
cg collapseoverlap ucsc_${build}_cytoBand.tsv ucsc_${build}_evofold.tsv ucsc_${build}_gwasCatalog.tsv ucsc_${build}_microsat.tsv ucsc_${build}_omimGene.tsv ucsc_${build}_oreganno.tsv ucsc_${build}_rmsk.tsv ucsc_${build}_simpleRepeat.tsv ucsc_${build}_targetScanS.tsv ucsc_${build}_tfbsConsSites.tsv ucsc_${build}_tRNAs.tsv ucsc_${build}_wgRna.tsv ucsc_${build}_vistaEnhancers.tsv ucsc_${build}_gad.tsv
cg collapseoverlap ucsc_${build}_phastConsElements28wayPlacMammal.tsv ucsc_${build}_phastConsElements28way.tsv ucsc_${build}_phastConsElements44way.tsv ucsc_${build}_rnaGene.tsv
cg regjoin ucsc_${build}_chainSelf.tsv > reg_${build}_chainSelf.tsv
cg regjoin ucsc_${build}_dgv.tsv > reg_${build}_dgv.tsv
cg regjoin ucsc_${build}_genomicSuperDups.tsv > reg_${build}_genomicSuperDups.tsv
# shorten some names
mv ucsc_${build}_phastConsElements28wayPlacMammal.tsv ucsc_${build}_phastCons28P.tsv
mv  reg_${build}_phastConsElements28wayPlacMammal.tsv  reg_${build}_phastCons28P.tsv
mv ucsc_${build}_phastConsElements28way.tsv ucsc_${build}_phastCons28.tsv
mv  reg_${build}_phastConsElements28way.tsv  reg_${build}_phastCons28.tsv
mv ucsc_${build}_phastConsElements44way.tsv ucsc_${build}_phastCons44.tsv
mv  reg_${build}_phastConsElements44way.tsv  reg_${build}_phastCons44.tsv
mv  reg_${build}_phastConsElements28wayPlacMammal.info  reg_${build}_phastCons28P.info
mv  reg_${build}_phastConsElements28way.info  reg_${build}_phastCons28.info
mv  reg_${build}_phastConsElements44way.info  reg_${build}_phastCons44.info

cg downloaddb ${dest} ${build} firstEF 
mv ucsc_${build}_firstEF.tsv reg_${build}_firstEF.tsv

# 1000 genomes
cg downloaddb ${dest} ${build} 1000g
cg downloaddb $dest hg19 1000glow
cg liftover ${dest}/hg19/var_hg19_1000glow.tsv ${dest}/liftover/hg19ToHg18.over.chain ${dest}/${build}/var_${build}_1000glow.tsv

# dbsnp
cg downloaddb ${dest} ${build} snp130
cg bgzip ${dest}/${build}/var_${build}_snp130.tsv
cg maketabix ${dest}/${build}/var_${build}_snp130.tsv
gunzip -c ${dest}/${build}/var_${build}_snp130.tsv.gz > ${dest}/${build}/var_${build}_snp130.tsv

# dbdnp132 with liftover
if [ ! -f ${dest}/hg19/var_hg19_snp132.tsv ]; then
	cg downloaddb ${dest}/hg19 hg19 snp132;
fi
cg liftover ${dest}/hg19/var_hg19_snp132.tsv ${dest}/liftover/hg19ToHg18.over.chain ${dest}/${build}/var_${build}_snp132lift.tsv
rm ${dest}/${build}/var_${build}_snp132lift.tsv.unmapped
cg maketabix ${dest}/${build}/var_${build}_snp132lift.tsv
gunzip -c ${dest}/${build}/var_${build}_snp132lift.tsv.gz > ${dest}/${build}/var_${build}_snp132lift.tsv

# genes
cg downloaddb ${dest} ${build} refGene ensGene knownGene genscan
mv ucsc_${build}_refGene.tsv gene_${build}_refGene.tsv
mv reg_${build}_refGene.info gene_${build}_refGene.info
cg maketabix gene_${build}_refGene.tsv
mv ucsc_${build}_ensGene.tsv gene_${build}_ensGene.tsv
mv reg_${build}_ensGene.info gene_${build}_ensGene.info
cg maketabix gene_${build}_ensGene.tsv
echo -e "genecol\tproteinID\ntranscriptcol\tname" > gene_${build}_knownGene.tsv.opt
mv ucsc_${build}_knownGene.tsv gene_${build}_knownGene.tsv
mv reg_${build}_knownGene.info gene_${build}_knownGene.info
cg maketabix gene_${build}_knownGene.tsv
# todo: also get kgXref, and translate
cg downloaddb ${dest} ${build} genscan
mv ucsc_${build}_genscan.tsv gene_${build}_genscan.tsv
mv reg_${build}_genscan.info gene_${build}_genscan.info
cg maketabix gene_${build}_genscan.tsv

# cg downloaddbinfo ${dest} ${build} rnaGene augustusAbinitio acembly firstEF
cg downloaddb ${dest} ${build} rnaGene augustusAbinitio acembly
mv ucsc_${build}_rnaGene.tsv reg_${build}_rnaGene.tsv
mv ucsc_${build}_augustusAbinitio.tsv gene_${build}_augustusAbinitio.tsv
rm reg_${build}_augustusAbinitio.info
mv ucsc_${build}_acembly.tsv gene_${build}_acembly.tsv
mv reg_${build}_acembly.info gene_${build}_acembly.info

# homopolymer
cd ${dest}/${build}
cg extracthomopolymers genome_${build}.ifas > reg_${build}_homopolymer.tsv.temp
mv reg_${build}_homopolymer.tsv.temp reg_${build}_homopolymer.tsv
echo -e "fields\t{base size}" > reg_${build}_homopolymer.tsv.opt

# encode
cd ${dest}/tmp/${build}
# enc_transcription
cg downloaddb ${dest}/tmp ${build} wgEncodeCaltechRnaSeqRawSignalRep1Gm12878CellLongpolyaBb12x75 wgEncodeCaltechRnaSeqRawSignalRep1H1hescCellPapBb2R2x75 wgEncodeCaltechRnaSeqRawSignalRep2Hepg2CellPapBb2R2x75 wgEncodeCaltechRnaSeqRawSignalRep1HuvecCellPapBb2R2x75 wgEncodeCaltechRnaSeqRawSignalRep1K562CellLongpolyaBb12x75 wgEncodeCaltechRnaSeqRawSignalRep1NhekCellPapBb2R2x75
cg ucscwiggle2reg -n 0.01 -p 1 -f 'log10($value)' ucsc_${build}_wgEncodeCaltechRnaSeqRawSignalRep1Gm12878CellLongpolyaBb12x75.tsv
cg ucscwiggle2reg -n 0.01 -p 1 -f 'log10($value)' ucsc_${build}_wgEncodeCaltechRnaSeqRawSignalRep1Gm12878CellLongpolyaBb12x75.tsv
cg ucscwiggle2reg -n 0.01 -p 1 -f 'log10($value)' ucsc_${build}_wgEncodeCaltechRnaSeqRawSignalRep1H1hescCellPapBb2R2x75.tsv
cg ucscwiggle2reg -n 0.01 -p 1 -f 'log10($value)' ucsc_${build}_wgEncodeCaltechRnaSeqRawSignalRep2Hepg2CellPapBb2R2x75.tsv
cg ucscwiggle2reg -n 0.01 -p 1 -f 'log10($value)' ucsc_${build}_wgEncodeCaltechRnaSeqRawSignalRep1HuvecCellPapBb2R2x75.tsv
cg ucscwiggle2reg -n 0.01 -p 1 -f 'log10($value)' ucsc_${build}_wgEncodeCaltechRnaSeqRawSignalRep1K562CellLongpolyaBb12x75.tsv
cg ucscwiggle2reg -n 0.01 -p 1 -f 'log10($value)' ucsc_${build}_wgEncodeCaltechRnaSeqRawSignalRep1NhekCellPapBb2R2x75.tsv
cg collapseoverlap -o ${dest}/${build}/reg_${build}_wgEncodelogRnaSeq.tsv reg_ucsc_${build}_wgEncodeCaltechRnaSeqRawSignalRep1Gm12878CellLongpolyaBb12x75.tsv reg_ucsc_${build}_wgEncodeCaltechRnaSeqRawSignalRep1H1hescCellPapBb2R2x75.tsv reg_ucsc_${build}_wgEncodeCaltechRnaSeqRawSignalRep2Hepg2CellPapBb2R2x75.tsv reg_ucsc_${build}_wgEncodeCaltechRnaSeqRawSignalRep1HuvecCellPapBb2R2x75.tsv reg_ucsc_${build}_wgEncodeCaltechRnaSeqRawSignalRep1K562CellLongpolyaBb12x75.tsv reg_ucsc_${build}_wgEncodeCaltechRnaSeqRawSignalRep1NhekCellPapBb2R2x75.tsv

# enc_H3K4Me1
cg downloaddb ${dest} ${build} wgEncodeBroadChipSeqSignalGm12878H3k4me1 wgEncodeBroadChipSeqSignalH1hescH3k4me1 wgEncodeBroadChipSeqSignalHmecH3k4me1 wgEncodeBroadChipSeqSignalHsmmH3k4me1 wgEncodeBroadChipSeqSignalHuvecH3k4me1 wgEncodeBroadChipSeqSignalK562H3k4me1 wgEncodeBroadChipSeqSignalNhekH3k4me1 wgEncodeBroadChipSeqSignalNhlfH3k4me1
cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalGm12878H3k4me1.tsv
cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalH1hescH3k4me1.tsv
cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalHmecH3k4me1.tsv
cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalHsmmH3k4me1.tsv
cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalHuvecH3k4me1.tsv
cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalK562H3k4me1.tsv
cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalNhekH3k4me1.tsv
cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalNhlfH3k4me1.tsv
cg collapseoverlap -o ${dest}/${build}/reg_${build}_wgEncodeH3k4me1.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalGm12878H3k4me1.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalH1hescH3k4me1.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalHmecH3k4me1.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalHsmmH3k4me1.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalHuvecH3k4me1.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalK562H3k4me1.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalNhekH3k4me1.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalNhlfH3k4me1.tsv

# enc_H3K27Ac
cg downloaddb ${dest} ${build} wgEncodeBroadChipSeqSignalGm12878H3k27ac wgEncodeBroadChipSeqSignalHepg2H3k27ac wgEncodeBroadChipSeqSignalHmecH3k27ac wgEncodeBroadChipSeqSignalHsmmH3k27ac wgEncodeBroadChipSeqSignalHuvecH3k27ac wgEncodeBroadChipSeqSignalK562H3k27ac wgEncodeBroadChipSeqSignalNhekH3k27ac wgEncodeBroadChipSeqSignalNhlfH3k27ac
cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalGm12878H3k27ac.tsv
cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalHepg2H3k27ac.tsv
cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalHmecH3k27ac.tsv
cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalHsmmH3k27ac.tsv
cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalHuvecH3k27ac.tsv
cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalK562H3k27ac.tsv
cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalNhekH3k27ac.tsv
cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalNhlfH3k27ac.tsv
# collapse
cg collapseoverlap -o ${dest}/${build}/reg_${build}_wgEncodeH3k27ac.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalGm12878H3k27ac.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalHepg2H3k27ac.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalHmecH3k27ac.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalHsmmH3k27ac.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalHuvecH3k27ac.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalK562H3k27ac.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalNhekH3k27ac.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalNhlfH3k27ac.tsv

# enc_H3K4Me3
cg downloaddb ${dest} ${build} wgEncodeBroadChipSeqSignalGm12878H3k4me3 wgEncodeBroadChipSeqSignalHepg2H3k4me3 wgEncodeBroadChipSeqSignalHmecH3k4me3 wgEncodeBroadChipSeqSignalHsmmH3k4me3 wgEncodeBroadChipSeqSignalHuvecH3k4me3 wgEncodeBroadChipSeqSignalK562H3k4me3 wgEncodeBroadChipSeqSignalNhekH3k4me3 wgEncodeBroadChipSeqSignalNhlfH3k4me3
cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalGm12878H3k4me3.tsv
cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalHepg2H3k4me3.tsv
cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalHmecH3k4me3.tsv
cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalHsmmH3k4me3.tsv
cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalHuvecH3k4me3.tsv
cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalK562H3k4me3.tsv
cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalNhekH3k4me3.tsv
cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalNhlfH3k4me3.tsv
cg collapseoverlap -o ${dest}/${build}/reg_${build}_wgEncodeH3k4me3.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalGm12878H3k4me3.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalHepg2H3k4me3.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalHmecH3k4me3.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalHsmmH3k4me3.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalHuvecH3k4me3.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalK562H3k4me3.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalNhekH3k4me3.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalNhlfH3k4me3.tsv

# DNase Clusters and Txn Factor ChIP
cg downloaddb ${dest} ${build} wgEncodeRegDnaseClustered wgEncodeRegTfbsClustered
cg collapseoverlap ucsc_${build}_wgEncodeRegDnaseClustered.tsv ucsc_${build}_wgEncodeRegTfbsClustered.tsv
mv reg_${build}_wgEncodeRegDnaseClustered.tsv reg_${build}_wgEncodeRegTfbsClustered.tsv ${dest}/${build}
