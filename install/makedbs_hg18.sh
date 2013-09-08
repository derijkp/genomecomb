export dest=/complgen/refseq
export build=hg18

# download hg18
# =============
#
mkdir -p ${dest}/${build}
cd ${dest}/${build}

# download genome
if [ -f "${dest}/${build}/genome_${build}.ifas" ]; then
	echo "Skipping ${dest}/${build}/genome_${build}.ifas, file exists"
else
	cg downloadgenome ${build} genome_${build}.ifas
	cg make_genomecindex genome_${build}.ifas
fi
mkdir extra
mv reg_genome_${build}.tsv extra/reg_${build}_fullgenome.tsv
if [ -f "${dest}/${build}/extra/reg_${build}_sequencedgenome.tsv" ]; then
	echo "Skipping ${dest}/${build}/extra/reg_${build}_sequencedgenome.tsv, file exists"
else
	cg calcsequencedgenome ${dest} ${build}
fi

# region databases (ucsc)
cg downloaddb ${dest} ${build} simpleRepeat microsat rmsk genomicSuperDups chainSelf
cg downloaddb ${dest} ${build} oreganno tRNAs targetScanS evofold
cg downloaddb ${dest} ${build} cytoBand dgv gwasCatalog wgRna vistaEnhancers gad tfbsConsSites
cg downloaddb ${dest} ${build} phastConsElements28way phastConsElements28wayPlacMammal phastConsElements44way rnaGene 
# you can explicitely download info on the databases using:
# cg downloaddbinfo ${dest} ${build} simpleRepeat microsat rmsk genomicSuperDups chainSelf
mv ucsc_${build}_gwasCatalog.tsv ucsc_${build}_gwasCatalog.tsv.temp
cg maketabix ${dest}/${build}/reg_${build}_rmsk.tsv
gunzip -c ${dest}/${build}/reg_${build}_rmsk.tsv.gz > ${dest}/${build}/reg_${build}_rmsk.tsv
cg maketabix ${dest}/${build}/reg_${build}_simpleRepeat.tsv
gunzip -c ${dest}/${build}/reg_${build}_simpleRepeat.tsv.gz > ${dest}/${build}/reg_${build}_simpleRepeat.tsv

cg select \
	-f 'chrom start end trait pValue pubMedID name bin author pubDate journal title initSample replSample region genes riskAllele riskAlFreq pValueDesc orOrBeta ci95 platform cnv' \
	-nh 'chrom start end name score pubMedID dbsnp bin author pubDate journal title initSample replSample region genes riskAllele riskAlFreq pValueDesc orOrBeta ci95 platform cnv' \
	ucsc_${build}_gwasCatalog.tsv.temp ucsc_${build}_gwasCatalog.tsv

# other databases
cg downloaddb ${dest} ${build} refLink kgXref
mv ucsc_${build}_kgXref.tsv other_${build}_kgXref.tsv
mv ucsc_${build}_refLink.tsv other_${build}_refLink.tsv

# collapse regions
cg collapseoverlap ucsc_${build}_cytoBand.tsv ucsc_${build}_evofold.tsv ucsc_${build}_gwasCatalog.tsv ucsc_${build}_microsat.tsv ucsc_${build}_oreganno.tsv ucsc_${build}_rmsk.tsv ucsc_${build}_simpleRepeat.tsv ucsc_${build}_targetScanS.tsv ucsc_${build}_tfbsConsSites.tsv ucsc_${build}_tRNAs.tsv ucsc_${build}_wgRna.tsv ucsc_${build}_vistaEnhancers.tsv ucsc_${build}_gad.tsv
cg collapseoverlap ucsc_${build}_phastConsElements28wayPlacMammal.tsv ucsc_${build}_phastConsElements28way.tsv ucsc_${build}_phastConsElements44way.tsv ucsc_${build}_rnaGene.tsv
cg regjoin ucsc_${build}_chainSelf.tsv > reg_${build}_chainSelf.tsv
cg regjoin ucsc_${build}_dgv.tsv > reg_${build}_dgv.tsv
cg regjoin ucsc_${build}_genomicSuperDups.tsv > reg_${build}_genomicSuperDups.tsv

cg downloaddb ${dest} ${build} firstEF 
mv ucsc_${build}_firstEF.tsv reg_${build}_firstEF.tsv

# 1000 genomes
cg downloaddb ${dest} ${build} 1000g
mkdir ${dest}/hg19
cg downloaddb $dest hg19 1000glow
if [ -f "${dest}/${build}/var_${build}_1000glow.tsv" ]; then
	echo "Skipping ${dest}/${build}/var_${build}_1000glow.tsv, file exists"
else
	cg liftover ${dest}/hg19/var_hg19_1000glow.tsv ${dest}/liftover/hg19ToHg18.over.chain ${dest}/${build}/var_${build}_1000glow.tsv
fi

# dbsnp
cg downloaddb ${dest} ${build} snp130
if [ -f "${dest}/${build}/var_${build}_snp130.tsv.gz.tbi" ]; then
	echo "Skipping ${dest}/${build}/var_${build}_snp130.tsv.tbi, file exists"
else
	cg maketabix ${dest}/${build}/var_${build}_snp130.tsv
	# gunzip -c ${dest}/${build}/var_${build}_snp130.tsv.gz > ${dest}/${build}/var_${build}_snp130.tsv
	cg select -f 'chrom start end type ref alt name freq' ${dest}/${build}/var_${build}_snp130.tsv.gz ${dest}/${build}/var_${build}_snp130.tsv
fi

# dbsnp135 with liftover
if [ -f "${dest}/${build}/var_${build}_snp135lift.tsv" ]; then
	echo "Skipping ${dest}/${build}/var_${build}_snp135lift.tsv, file exists"
else
	if [ ! -f ${dest}/hg19/var_hg19_snp135.tsv ]; then
		cg downloaddb ${dest} hg19 snp135;
	fi
	cg liftover ${dest}/hg19/var_hg19_snp135.tsv ${dest}/liftover/hg19ToHg18.over.chain ${dest}/${build}/var_${build}_snp135lift.tsv
	rm ${dest}/${build}/var_${build}_snp135lift.tsv.unmapped
fi
if [ -f "${dest}/${build}/var_${build}_snp135lift.tsv.gz.tbi" ]; then
	echo "Skipping ${dest}/${build}/var_${build}_snp135lift.tsv.gz.tbi, file exists"
else
	cg maketabix ${dest}/${build}/var_${build}_snp135lift.tsv
	gunzip -c ${dest}/${build}/var_${build}_snp135lift.tsv.gz > ${dest}/${build}/var_${build}_snp135lift.tsv
fi

# genes
for database in refGene ensGene knownGene genscan augustusAbinitio acembly wgEncodeGencodeManualV3 wgEncodeGencodeAutoV3; do
    if [ -f "${dest}/${build}/gene_${build}_${database}.tsv" ]; then
        echo "Skipping ${dest}/${build}/gene_${build}_${database}.tsv, file exists"
    else
        cg downloaddb ${dest} ${build} $database
        cg select -s - ucsc_${build}_${database}.tsv gene_${build}_${database}.tsv
        mv reg_${build}_${database}.info gene_${build}_${database}.info
        cg maketabix gene_${build}_${database}.tsv
        gunzip -c gene_${build}_${database}.tsv.gz > gene_${build}_${database}.tsv
    fi
done

cg downloaddb ${dest} ${build} rnaGene
mv ucsc_${build}_rnaGene.tsv reg_${build}_rnaGene.tsv

# homopolymer
cd ${dest}/${build}
if [ -f "reg_${build}_homopolymer.tsv" ]; then
	echo "Skipping reg_${build}_homopolymer.tsv, file exists"
else
	cg extracthomopolymers genome_${build}.ifas > reg_${build}_homopolymer.tsv.temp
	mv reg_${build}_homopolymer.tsv.temp reg_${build}_homopolymer.tsv
fi
echo -e "fields\t{base size}" > reg_${build}_homopolymer.tsv.opt

# encode
cd ${dest}/tmp/${build}
# enc_transcription
if [ -f "${dest}/${build}/reg_${build}_wgEncodelogRnaSeq.tsv" ]; then
	echo "Skipping ${dest}/${build}/reg_${build}_wgEncodelogRnaSeq.tsv, file exists"
else
	cg downloaddb ${dest}/tmp ${build} wgEncodeCaltechRnaSeqRawSignalRep1Gm12878CellLongpolyaBb12x75 wgEncodeCaltechRnaSeqRawSignalRep1H1hescCellPapBb2R2x75 wgEncodeCaltechRnaSeqRawSignalRep2Hepg2CellPapBb2R2x75 wgEncodeCaltechRnaSeqRawSignalRep1HuvecCellPapBb2R2x75 wgEncodeCaltechRnaSeqRawSignalRep1K562CellLongpolyaBb12x75 wgEncodeCaltechRnaSeqRawSignalRep1NhekCellPapBb2R2x75
	cg ucscwiggle2reg -n 0.01 -p 1 -f 'log10($value)' ucsc_${build}_wgEncodeCaltechRnaSeqRawSignalRep1Gm12878CellLongpolyaBb12x75.tsv
	cg ucscwiggle2reg -n 0.01 -p 1 -f 'log10($value)' ucsc_${build}_wgEncodeCaltechRnaSeqRawSignalRep1Gm12878CellLongpolyaBb12x75.tsv
	cg ucscwiggle2reg -n 0.01 -p 1 -f 'log10($value)' ucsc_${build}_wgEncodeCaltechRnaSeqRawSignalRep1H1hescCellPapBb2R2x75.tsv
	cg ucscwiggle2reg -n 0.01 -p 1 -f 'log10($value)' ucsc_${build}_wgEncodeCaltechRnaSeqRawSignalRep2Hepg2CellPapBb2R2x75.tsv
	cg ucscwiggle2reg -n 0.01 -p 1 -f 'log10($value)' ucsc_${build}_wgEncodeCaltechRnaSeqRawSignalRep1HuvecCellPapBb2R2x75.tsv
	cg ucscwiggle2reg -n 0.01 -p 1 -f 'log10($value)' ucsc_${build}_wgEncodeCaltechRnaSeqRawSignalRep1K562CellLongpolyaBb12x75.tsv
	cg ucscwiggle2reg -n 0.01 -p 1 -f 'log10($value)' ucsc_${build}_wgEncodeCaltechRnaSeqRawSignalRep1NhekCellPapBb2R2x75.tsv
	cg collapseoverlap -o ${dest}/${build}/reg_${build}_wgEncodelogRnaSeq.tsv reg_ucsc_${build}_wgEncodeCaltechRnaSeqRawSignalRep1Gm12878CellLongpolyaBb12x75.tsv reg_ucsc_${build}_wgEncodeCaltechRnaSeqRawSignalRep1H1hescCellPapBb2R2x75.tsv reg_ucsc_${build}_wgEncodeCaltechRnaSeqRawSignalRep2Hepg2CellPapBb2R2x75.tsv reg_ucsc_${build}_wgEncodeCaltechRnaSeqRawSignalRep1HuvecCellPapBb2R2x75.tsv reg_ucsc_${build}_wgEncodeCaltechRnaSeqRawSignalRep1K562CellLongpolyaBb12x75.tsv reg_ucsc_${build}_wgEncodeCaltechRnaSeqRawSignalRep1NhekCellPapBb2R2x75.tsv
fi
# make info file
cg downloaddbinfo ${dest}/tmp ${build} wgEncodeCaltechRnaSeq
head -1 ${dest}/tmp/${build}/reg_hg18_wgEncodeCaltechRnaSeq.info > ${dest}/${build}/reg_hg18_wgEncodeCaltechRnaSeq.info
echo $'\n\n== Agregation info ==\n\nThis file combines the data from 7 lines' >> ${dest}/${build}/reg_hg18_wgEncodeCaltechRnaSeq.info
echo $'score contains the logarithm of the highest score, with a precision of 1 digit after the decimal point.' >> ${dest}/${build}/reg_hg18_wgEncodeCaltechRnaSeq.info
echo $'num contains the number of lines for which the original value is >= 0.01 (thus log >= -2.0)\n' >> ${dest}/${build}/reg_hg18_wgEncodeCaltechRnaSeq.info
head -n -2 ${dest}/tmp/${build}/reg_hg18_wgEncodeCaltechRnaSeq.info >> ${dest}/${build}/reg_hg18_wgEncodeCaltechRnaSeq.info

# enc_H3K4Me1
if [ -f "${dest}/${build}/reg_${build}_wgEncodeH3k4me1.tsv" ]; then
	echo "Skipping ${dest}/${build}/reg_${build}_wgEncodeH3k4me1.tsv, file exists"
else
	cg downloaddb ${dest}/tmp ${build} wgEncodeBroadChipSeqSignalGm12878H3k4me1 wgEncodeBroadChipSeqSignalH1hescH3k4me1 wgEncodeBroadChipSeqSignalHmecH3k4me1 wgEncodeBroadChipSeqSignalHsmmH3k4me1 wgEncodeBroadChipSeqSignalHuvecH3k4me1 wgEncodeBroadChipSeqSignalK562H3k4me1 wgEncodeBroadChipSeqSignalNhekH3k4me1 wgEncodeBroadChipSeqSignalNhlfH3k4me1
	cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalGm12878H3k4me1.tsv
	cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalH1hescH3k4me1.tsv
	cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalHmecH3k4me1.tsv
	cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalHsmmH3k4me1.tsv
	cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalHuvecH3k4me1.tsv
	cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalK562H3k4me1.tsv
	cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalNhekH3k4me1.tsv
	cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalNhlfH3k4me1.tsv
	cg collapseoverlap -o ${dest}/${build}/reg_${build}_wgEncodeH3k4me1.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalGm12878H3k4me1.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalH1hescH3k4me1.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalHmecH3k4me1.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalHsmmH3k4me1.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalHuvecH3k4me1.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalK562H3k4me1.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalNhekH3k4me1.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalNhlfH3k4me1.tsv
fi
# make info file
cg downloaddbinfo ${dest}/tmp ${build} wgEncodeRegMarkEnhH3k4me1
head -1 ${dest}/tmp/${build}/reg_hg18_wgEncodeRegMarkEnhH3k4me1.info > ${dest}/${build}/reg_hg18_wgEncodeH3k4me1.info
echo $'\n\n== Agregation info ==\n' >> ${dest}/${build}/reg_hg18_wgEncodeH3k4me1.info
echo $'This file combines the data from 8 cell lines' >> ${dest}/${build}/reg_hg18_wgEncodeH3k4me1.info
echo $'score contains the highest score rounded up the next 5 fold.' >> ${dest}/${build}/reg_hg18_wgEncodeH3k4me1.info
echo $'num contains the number of cell lines for which the score is >= 10\n' >> ${dest}/${build}/reg_hg18_wgEncodeH3k4me1.info
head -n -2 ${dest}/tmp/${build}/reg_hg18_wgEncodeRegMarkEnhH3k4me1.info >> ${dest}/${build}/reg_hg18_wgEncodeH3k4me1.info

# enc_H3K27Ac
if [ -f "${dest}/${build}/reg_${build}_wgEncodeH3k27ac.tsv" ]; then
	echo "Skipping ${dest}/${build}/reg_${build}_wgEncodeH3k27ac.tsv, file exists"
else
	cg downloaddb ${dest}/tmp ${build} wgEncodeBroadChipSeqSignalGm12878H3k27ac wgEncodeBroadChipSeqSignalHepg2H3k27ac wgEncodeBroadChipSeqSignalHmecH3k27ac wgEncodeBroadChipSeqSignalHsmmH3k27ac wgEncodeBroadChipSeqSignalHuvecH3k27ac wgEncodeBroadChipSeqSignalK562H3k27ac wgEncodeBroadChipSeqSignalNhekH3k27ac wgEncodeBroadChipSeqSignalNhlfH3k27ac
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
fi
# make info file
cg downloaddbinfo ${dest}/tmp ${build} wgEncodeRegMarkEnhH3k27ac
head -1 ${dest}/tmp/${build}/reg_hg18_wgEncodeRegMarkEnhH3k27ac.info > ${dest}/${build}/reg_hg18_wgEncodeH3k27ac.info
echo $'\n\n== Agregation info ==\n' >> ${dest}/${build}/reg_hg18_wgEncodeH3k27ac.info
echo $'This file combines the data from 8 cell lines' >> ${dest}/${build}/reg_hg18_wgEncodeH3k27ac.info
echo $'score contains the highest score rounded up the next 5 fold.' >> ${dest}/${build}/reg_hg18_wgEncodeH3k27ac.info
echo $'num contains the number of cell lines for which the score is >= 10\n' >> ${dest}/${build}/reg_hg18_wgEncodeH3k27ac.info
head -n -2 ${dest}/tmp/${build}/reg_hg18_wgEncodeRegMarkEnhH3k27ac.info >> ${dest}/${build}/reg_hg18_wgEncodeH3k27ac.info

# enc_H3K4Me3
if [ -f "${dest}/${build}/reg_${build}_wgEncodeH3k4me3.tsv" ]; then
	echo "Skipping ${dest}/${build}/reg_${build}_wgEncodeH3k4me3.tsv, file exists"
else
	cg downloaddb ${dest}/tmp ${build} wgEncodeBroadChipSeqSignalGm12878H3k4me3 wgEncodeBroadChipSeqSignalHepg2H3k4me3 wgEncodeBroadChipSeqSignalHmecH3k4me3 wgEncodeBroadChipSeqSignalHsmmH3k4me3 wgEncodeBroadChipSeqSignalHuvecH3k4me3 wgEncodeBroadChipSeqSignalK562H3k4me3 wgEncodeBroadChipSeqSignalNhekH3k4me3 wgEncodeBroadChipSeqSignalNhlfH3k4me3
	cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalGm12878H3k4me3.tsv
	cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalHepg2H3k4me3.tsv
	cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalHmecH3k4me3.tsv
	cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalHsmmH3k4me3.tsv
	cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalHuvecH3k4me3.tsv
	cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalK562H3k4me3.tsv
	cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalNhekH3k4me3.tsv
	cg ucscwiggle2reg -n 10 -p 0 -f '5*(($value+4)/5)' ucsc_${build}_wgEncodeBroadChipSeqSignalNhlfH3k4me3.tsv
	cg collapseoverlap -o ${dest}/${build}/reg_${build}_wgEncodeH3k4me3.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalGm12878H3k4me3.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalHepg2H3k4me3.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalHmecH3k4me3.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalHsmmH3k4me3.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalHuvecH3k4me3.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalK562H3k4me3.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalNhekH3k4me3.tsv reg_ucsc_${build}_wgEncodeBroadChipSeqSignalNhlfH3k4me3.tsv
fi
# make info file
cg downloaddbinfo ${dest}/tmp ${build} wgEncodeRegMarkPromoter
head -1 ${dest}/tmp/${build}/reg_hg18_wgEncodeRegMarkPromoter.info > ${dest}/${build}/reg_hg18_wgEncodeH3k4me3.info
echo $'\n\n== Agregation info ==\n' >> ${dest}/${build}/reg_hg18_wgEncodeH3k4me3.info
echo $'This file combines the data from 9 cell lines' >> ${dest}/${build}/reg_hg18_wgEncodeH3k4me3.info
echo $'score contains the highest score rounded up the next 5 fold.' >> ${dest}/${build}/reg_hg18_wgEncodeH3k4me3.info
echo $'num contains the number of cell lines for which the score is >= 10\n' >> ${dest}/${build}/reg_hg18_wgEncodeH3k4me3.info
head -n -2 ${dest}/tmp/${build}/reg_hg18_wgEncodeRegMarkPromoter.info >> ${dest}/${build}/reg_hg18_wgEncodeH3k4me3.info

# DNase Clusters and Txn Factor ChIP
if [ -f "${dest}/${build}/reg_${build}_wgEncodeRegDnaseClustered.tsv" ]; then
	echo "Skipping ${dest}/${build}/reg_${build}_wgEncodeRegDnaseClustered.tsv, file exists"
else
	cg downloaddb ${dest}/tmp ${build} wgEncodeRegDnaseClustered
	cg collapseoverlap ucsc_${build}_wgEncodeRegDnaseClustered.tsv
	mv reg_${build}_wgEncodeRegDnaseClustered.tsv ${dest}/${build}
fi

if [ -f "${dest}/${build}/reg_${build}_wgEncodeRegTfbsClustered.tsv" ]; then
	echo "Skipping ${dest}/${build}/reg_${build}_wgEncodeRegTfbsClustered.tsv, file exists"
else
	cg downloaddb ${dest}/tmp ${build} wgEncodeRegTfbsClustered
	cg collapseoverlap ucsc_${build}_wgEncodeRegTfbsClustered.tsv
	mv reg_${build}_wgEncodeRegTfbsClustered.tsv ${dest}/${build}
fi

# mirbase
wget --tries=45 --directory-prefix=${dest}/tmp/hg19 ftp://mirbase.org/pub/mirbase/19/genomes/hsa.gff2
cg gff2sft ${dest}/tmp/hg19/hsa.gff2 ${dest}/tmp/hg19/reg_hg19_mirbase.tsv.temp
cg select -s - ${dest}/tmp/hg19/reg_hg19_mirbase.tsv.temp ${dest}/tmp/hg19/reg_hg19_mirbase.tsv.temp2
rm ${dest}/tmp/${build}/reg_${build}_mirbase.tsv.temp2 || true
cg liftover ${dest}/tmp/hg19/reg_hg19_mirbase.tsv.temp2 ${dest}/liftover/hg19ToHg18.over.chain ${dest}/tmp/${build}/reg_${build}_mirbase.tsv.temp2
mv ${dest}/tmp/${build}/reg_${build}_mirbase.tsv.temp2 reg_${build}_mirbase.tsv
wget ftp://mirbase.org/pub/mirbase/19/README
mv README reg_${build}_mirbase.info
