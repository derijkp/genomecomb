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
if [ -f "${dest}/${build}/extra/reg_${build}_sequencedgenome.tsv" ]; then
	echo "Skipping ${dest}/${build}/extra/reg_${build}_sequencedgenome.tsv, file exists"
else
	cg calcsequencedgenome ${dest} ${build}
fi

# region databases (ucsc)
cg downloaddb ${dest} ${build} simpleRepeat microsat rmsk genomicSuperDups chainSelf
cg downloaddb ${dest} ${build} oreganno tRNAs targetScanS evofold
cg downloaddb ${dest} ${build} cytoBand dgv gwasCatalog wgRna vistaEnhancers gad tfbsConsSites
cg downloaddb ${dest} ${build} phastConsElements46way phastConsElements46wayPlacental phastConsElements46wayPrimates
# you can explicitely download info on the databases using:
# cg downloaddbinfo ${dest} ${build} simpleRepeat microsat rmsk genomicSuperDups chainSelf
mv ucsc_${build}_gwasCatalog.tsv ucsc_${build}_gwasCatalog.tsv.temp
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
cg collapseoverlap ucsc_${build}_phastConsElements46way.tsv ucsc_${build}_phastConsElements46wayPlacental.tsv ucsc_${build}_phastConsElements46wayPrimates.tsv
cg regjoin ucsc_${build}_chainSelf.tsv > reg_${build}_chainSelf.tsv
cg regjoin ucsc_${build}_dgv.tsv > reg_${build}_dgv.tsv
cg regjoin ucsc_${build}_genomicSuperDups.tsv > reg_${build}_genomicSuperDups.tsv

# 1000 genomes
cg downloaddb ${dest} hg18 1000g
cg liftover ${dest}/hg18/var_hg18_1000gCHBxJPT.tsv ${dest}/liftover/hg18ToHg19.over.chain ${dest}/${build}/var_${build}_1000gCHBxJPT.tsv
cg liftover ${dest}/hg18/var_hg18_1000gCEU.tsv ${dest}/liftover/hg18ToHg19.over.chain ${dest}/${build}/var_${build}_1000gCEU.tsv
cg liftover ${dest}/hg18/var_hg18_1000gYRI.tsv ${dest}/liftover/hg18ToHg19.over.chain ${dest}/${build}/var_${build}_1000gYRI.tsv
rm var_hg19_1000gCEU.tsv.unmapped var_hg19_1000gCHBxJPT.tsv.unmapped var_hg19_1000gYRI.tsv.unmapped

cg downloaddb ${dest} ${build} 1000glow

# dbsnp
if [ -f "${dest}/${build}/var_${build}_snp135.tsv" ]; then
	echo "Skipping ${dest}/${build}/var_${build}_snp135.tsv, file exists"
else
	cg downloaddb ${dest} ${build} snp135
fi
if [ -f "${dest}/${build}/var_${build}_snp135.tsv.tbi" ]; then
	echo "Skipping ${dest}/${build}/var_${build}_snp135.tsv.tbi, file exists"
else
	cg bgzip ${dest}/${build}/var_${build}_snp135.tsv
	cg maketabix ${dest}/${build}/var_${build}_snp135.tsv.gz
	# gunzip -c ${dest}/${build}/var_${build}_snp135.tsv.gz > ${dest}/${build}/var_${build}_snp135.tsv
	cg select -f 'chrom start end type ref alt name freq' ${dest}/${build}/var_${build}_snp135.tsv.gz ${dest}/${build}/var_${build}_snp135.tsv
fi

# genes
for database in refGene ensGene knownGene genscan acembly ; do
    if [ -f "${dest}/${build}/gene_${build}_${database}.tsv" ]; then
        echo "Skipping ${dest}/${build}/gene_${build}_${database}.tsv, file exists"
    else
        cg downloaddb ${dest} ${build} $database
        cg select -s 'chrom start end' -q '$chrom !~ /_/' ucsc_${build}_${database}.tsv gene_${build}_${database}.tsv
        mv reg_${build}_${database}.info gene_${build}_${database}.info
        cg maketabix gene_${build}_${database}.tsv
        gunzip -c gene_${build}_${database}.tsv.gz > gene_${build}_${database}.tsv
    fi
done

# homopolymer
cd ${dest}/${build}
if [ -f "reg_${build}_homopolymer.tsv" ]; then
	echo "Skipping reg_${build}_homopolymer.tsv, file exists"
else
	cg extracthomopolymers genome_${build}.ifas > reg_${build}_homopolymer.tsv.temp
	mv reg_${build}_homopolymer.tsv.temp reg_${build}_homopolymer.tsv
        cg maketabix reg_${build}_homopolymer.tsv
        gunzip -c reg_${build}_homopolymer.tsv.gz > reg_${build}_homopolymer.tsv
fi
echo -e "fields\t{base size}" > reg_${build}_homopolymer.tsv.opt
