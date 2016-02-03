#!/bin/sh
# the next line restarts using tclsh \
exec cg source "$0" "$@"

set build hg19
if {![info exists argv]} {set argv {}}
set argv [job_init {*}$argv]
if {[llength $argv]} {
	set dest [lindex $argv 0]
} else {
	set dest /complgen/refseq
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
	cg downloadgenome ${build} genome_${build}.ifas {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 M X Y}
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
	file delete -force ucsc_${build}_gwasCatalog.tsv.temp
	cg select \
		-f {chrom start end trait pValue pubMedID name bin author pubDate journal title initSample replSample region genes riskAllele riskAlFreq pValueDesc orOrBeta ci95 platform cnv} \
		-nh {chrom start end name score pubMedID dbsnp bin author pubDate journal title initSample replSample region genes riskAllele riskAlFreq pValueDesc orOrBeta ci95 platform cnv} \
		ucsc_${build}_gwasCatalog.tsv ucsc_${build}_gwasCatalog.tsv.temp
	cg collapseoverlap ucsc_${build}_gwasCatalog.tsv.temp
	file rename -force reg_${build}_gwasCatalog.tsv.temp $target
	file delete ucsc_${build}_gwasCatalog.tsv
	file delete ucsc_${build}_gwasCatalog.tsv.temp
}

foreach db {
	rmsk simpleRepeat
} {
	job maketabix_${build}_$db -deps {reg_${build}_${db}.tsv} -targets {reg_${build}_${db}.tsv.gz.tbi reg_${build}_${db}.tsv.gz} -vars {build db} -code {
		cg maketabix reg_${build}_${db}.tsv
		exec gunzip -c reg_${build}_${db}.tsv.gz > reg_${build}_${db}.tsv
	}
}

## 1000 genomes
## 1000g (hg18)
#job 1000gh18 -targets {$dest/hg18/var_hg18_1000gCHBxJPT.tsv $dest/hg18/var_hg18_1000gCEU.tsv $dest/hg18/var_hg18_1000gYRI.tsv} -vars {dest} -code {
#	cd ${dest}/hg18
#	cg downloaddb ${dest} hg18 1000g
#}
#
#job 1000gliftover -deps {$dest/hg18/var_hg18_1000gCHBxJPT.tsv $dest/hg18/var_hg18_1000gCEU.tsv $dest/hg18/var_hg18_1000gYRI.tsv} -targets {var_${build}_1000gCHBxJPT.tsv var_${build}_1000gCEU.tsv var_${build}_1000gYRI.tsv} -vars {dest build} -code {
#	cg liftover ${dest}/hg18/var_hg18_1000gCHBxJPT.tsv var_${build}_1000gCHBxJPT.tsv ${dest}/liftover/hg18ToHg19.over.chain
#	cg liftover ${dest}/hg18/var_hg18_1000gCEU.tsv var_${build}_1000gCEU.tsv ${dest}/liftover/hg18ToHg19.over.chain
#	cg liftover ${dest}/hg18/var_hg18_1000gYRI.tsv var_${build}_1000gYRI.tsv ${dest}/liftover/hg18ToHg19.over.chain
#	file delete var_hg19_1000gCEU.tsv.unmapped var_hg19_1000gCHBxJPT.tsv.unmapped var_hg19_1000gYRI.tsv.unmapped
#}
#
#job 1000glow -targets {${dest}/hg19/var_hg19_1000glow.tsv ${dest}/hg19/extra/var_hg19_1000glow.tsv.opt} -vars {dest} -code {
#	cd $dest/hg19
#	cg downloaddb ${dest} hg19 1000glow
#	cplinked $target $dest/$build/extra/var_${build}_1000glow.tsv
#	file_write $dest/$build/extra/var_${build}_1000glow.tsv.opt "fields\t{AMR_AF ASN_AF AFR_AF EUR_AF}\n"
#}

job 1000g3 -targets {${dest}/hg19/var_hg19_1000g3.tsv ${dest}/hg19/extra/var_hg19_1000g3.tsv.opt} -vars {dest build} -code {
	cd $dest/hg19
	cg downloaddb ${dest} hg19 1000g3
	cplinked $target $dest/$build/extra/var_${build}_1000g3.tsv
	file_write $dest/$build/extra/var_${build}_1000g3.tsv.opt "fields\t{EUR_AF AMR_AF EAS_AF SAS_AF AFR_AF}\n"
}

# dbsnp
job dbsnp144 -targets {${dest}/$build/var_${build}_snp144.tsv} -vars {dest build} -code {
	cd $dest/$build
	cg downloaddb ${dest} $build snp144
}

job dbsnp144Common -targets {${dest}/${build}/var_${build}_snp144Common.tsv} -vars {dest build} -code {
	cd $dest/${build}
	cg downloaddb ${dest} ${build} snp144Common
}

foreach db {
	snp144 snp144Common
} {
	job maketabix_${build}_$db -deps {var_${build}_${db}.tsv} -targets {var_${build}_${db}.tsv.gz.tbi var_${build}_${db}.tsv.gz} -vars {dest build db} -code {
		cg maketabix var_${build}_${db}.tsv
		# exec gunzip -c var_${build}_${db}.tsv.gz > var_${build}_${db}.tsv
		cg select -f {chrom start end type ref alt name freq} ${dest}/${build}/var_${build}_${db}.tsv.gz ${dest}/${build}/var_${build}_${db}.tsv
	}
}

job clinvar -targets {${dest}/${build}/var_${build}_clinvar.tsv} -vars {dest build} -code {
	cd $dest/tmp/${build}
	wgetfile ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
	cg vcf2tsv clinvar.vcf.gz clinvar.tsv
	set f [open clinvar.tsv]
	set header [tsv_open $f comment]
	close $f
	if {![regexp reference=GRCh37 $comment]} {
		error "clinvar.vcf.gz is from a different reference genome version"
	}
	cg collapsealleles clinvar.tsv > $target.temp
	file_write [gzroot $target].opt.temp "fields\t{CLNACC CLNDBN}\nheaderfields\t{clinvar_acc clinvar_disease}\n"
	file rename -force [gzroot $target].opt.temp [gzroot $target].opt
	file rename -force $target.temp $target
}

# genes
foreach db {
	refGene ensGene knownGene wgEncodeGencodeCompV19 genscan acembly
} {
	job gene_${build}_$db -targets {gene_${build}_${db}.tsv gene_${build}_${db}.tsv.gz.tbi gene_${build}_${db}.tsv.gz} -vars {dest build db} -code {
		if {$db eq "wgEncodeGencodeCompV19"} {set dbname gencode} else {set dbname $db}
		if {$db in "genscan acembly"} {set geneidcol name} else {set geneidcol name2}
	        cg downloaddb ${dest} ${build} $db
		if {$db in "ensGene knownGene"} {
			unset -nocomplain a
			if {$db eq "ensGene"} {
				cg downloaddb ${dest} ${build} ensemblToGeneName
				array set a [split [string trim [file_read ucsc_hg19_ensemblToGeneName.tsv]] "\n\t"]
			} else {
				cg downloaddb ${dest} ${build} kgXref
				array set a [split [string trim [cg select -f {kgID geneSymbol} ucsc_hg19_kgXref.tsv]] "\n\t"]
			}
			catch {close $f} ; catch {close $o}
			set f [open ucsc_${build}_${db}.tsv]
			set o [open gene_${build}_${dbname}.tsv.temp2 w]
			set header [split [gets $f] \t]
			set namepos [lsearch $header name]
			lappend header geneid
			puts $o [join $header \t]
			while {[gets $f line] != -1} {
				set line [split $line \t]
				set name [lindex $line $namepos]
				set geneid [get a($name) $name]
				lappend line $geneid
				puts $o [join $line \t]
			}
			close $o
			close $f
		        cg select -s - -f {chrom start end strand geneid *} gene_${build}_${dbname}.tsv.temp2 gene_${build}_${dbname}.tsv.temp
			file delete gene_${build}_${dbname}.tsv.temp2 ucsc_hg19_ensemblToGeneName.tsv
		} else {
		        cg select -s - -f [list chrom start end strand "geneid=\$$geneidcol" *] ucsc_${build}_${db}.tsv gene_${build}_${dbname}.tsv.temp
		}
		file rename -force gene_${build}_${dbname}.tsv.temp gene_${build}_${dbname}.tsv
	        file rename -force reg_${build}_${db}.info gene_${build}_${dbname}.info
	        cg maketabix gene_${build}_${dbname}.tsv
	        exec gunzip -c gene_${build}_${dbname}.tsv.gz > gene_${build}_${dbname}.tsv.temp
		file rename -force gene_${build}_${dbname}.tsv.temp gene_${build}_${dbname}.tsv
		file delete ucsc_${build}_${dbname}.tsv
		cg index gene_${build}_${dbname}.tsv
	}
}

job reg_${build}_genes \
-deps {gene_${build}_refGene.tsv gene_${build}_ensGene.tsv gene_${build}_knownGene.tsv gene_${build}_gencode.tsv} \
-targets {extra/reg_${build}_genes.tsv} \
-code {
	exec cg cat -fields {chrom start end geneid} {*}$deps | cg select -s {chrom start end geneid} -f {chrom {start=$start - 2000} {end=$end + 2000} geneid} | cg regcollapse > $target.temp
	file rename -force $target.temp $target
}

job reg_${build}_phenotype -deps {extra/reg_${build}_genes.tsv} \
-targets {reg_${build}_phenotype.tsv geneannot_${build}_phenotype.tsv} -vars {dest build} -code {
	# get target2, uses biomart (ensembl geneset) and clinvar for gene-phenotype correlations
	cg downloaddb ${dest} ${build} phenotype
	# to reg file
	cg geneannot2reg $dep $target2 $target.temp
	file rename -force $target.temp $target
}

job reg_${build}_go -deps {extra/reg_${build}_genes.tsv} \
-targets {extra/reg_${build}_go.tsv extra/geneannot_${build}_go.tsv} -code {
	cg downloadmart $target2 hsapiens_gene_ensembl gene_ensembl_config {hgnc_symbol name_1006 namespace_1003}
	cg geneannot2reg $dep $target2 $target.temp
	file rename -force $target.temp $target
	file delete temp_go.tsv.temp
}

foreach {name id disease pattern} {
	ad	C0002395	{Alzheimer's Disease} Alzheimer
	ftd	C0338451 {Frontotemporal Dementia} frontotemp
	als	C0002736	{Amyotrophic Lateral Sclerosis} {amyotrophic lateral}
	pd	C0030567	{Parkinson Disease} parkinson
	ep	C0014544	{Epilepsy} epilep
	cmt	C0007959	{Charcot-Marie-Tooth Disease} charcot
	pn	C0031117	{Peripheral Neuropathy} {peripher.*neuropat|motor.*neuropat|neuropat.*motor|peripher.*nerv}
} {
	job reg_${build}_biograph_$name -vars {id pattern name disease} \
	-deps {extra/reg_${build}_genes.tsv geneannot_${build}_phenotype.tsv} \
	-targets {extra/reg_${build}_biograph_$name.tsv extra/geneannot_${build}_biograph_$name.tsv extra/reg_${build}_biograph_$name.info} -code {
		set genefile [tempfile]
		cg select -q "\$phenotype_description regexp \"(?i)$pattern\"" $dep2 $genefile
		cg downloadbiograph $target2 $id $genefile
		cg geneannot2reg $dep $target2 $target.temp
		file rename -force $target.temp $target
		set temp "gene ranking\n============\ngene ranking for $disease\n\nThe following scores will be given to variants in genes that"
		append temp "0: are directly linked to the disease in the ensembl gene database or in the clinvar database\n"
		append temp "1: have a known relation to the disease in biograph (can be association, etc.)\n"
		append temp ">1: position in genelist sorted according to biograph prioritization (there are 18182 genes)"
		file_write $target3 $temp
	}
}

# homopolymer
job reg_${build}_homopolymer -deps {genome_${build}.ifas} -targets {reg_${build}_homopolymer.tsv reg_${build}_homopolymer.tsv.gz reg_${build}_homopolymer.tsv.gz.tbi reg_${build}_homopolymer.tsv.opt} -vars {dest build db} -code {
	cg extracthomopolymers genome_${build}.ifas > reg_${build}_homopolymer.tsv.temp
	file rename -force reg_${build}_homopolymer.tsv.temp reg_${build}_homopolymer.tsv
        cg maketabix reg_${build}_homopolymer.tsv
        gunzip reg_${build}_homopolymer.tsv.gz reg_${build}_homopolymer.tsv
	file_write reg_${build}_homopolymer.tsv.opt "fields\t{base size}\n"
}

# mirdmg
job reg_${build}_mirdmg -targets {$dest/${build}/mir_${build}_mirdmg.tsv $dest/${build}/mir_${build}_mirdmg.info} -vars {dest build db} -code {
	wgetfile http://mirna.bioinf.be/public/mir_hg19_mirdmg.tsv mir_hg19_mirdmg.tsv.temp
	file rename mir_hg19_mirdmg.tsv.temp mir_hg19_mirdmg.tsv
}

# mirbase
#job reg_${build}_mirbase -targets {$dest/${build}/reg_${build}_mirbase.tsv $dest/${build}/reg_${build}_mirbase.tsv.opt $dest/${build}/reg_${build}_mirbase.info} -vars {dest build db} -code {
#	set organism hsa
#	cd $dest/${build}
#	file_write $dest/${build}/reg_${build}_mirbase.tsv.opt "fields\t{ID}\n"
#	exec -ignorestderr wget -c --tries=45 --directory-prefix=${dest}/tmp/${build} ftp://mirbase.org/pub/mirbase/20/genomes/$organism.gff2
#	cg gff2sft ${dest}/tmp/${build}/$organism.gff2 ${dest}/tmp/${build}/reg_${build}_mirbase.tsv.temp
#	cg select -s - ${dest}/tmp/${build}/reg_${build}_mirbase.tsv.temp ${dest}/tmp/${build}/reg_${build}_mirbase.tsv.temp2
#	file rename -force ${dest}/tmp/${build}/reg_${build}_mirbase.tsv.temp2 reg_${build}_mirbase.tsv
#	exec -ignorestderr wget -c ftp://mirbase.org/pub/mirbase/20/README
#	file rename -force README reg_${build}_mirbase.info
#}

# exome variant server
job reg_hg19_evs -targets {$dest/hg19/extra/var_hg19_evs.tsv $dest/hg19/extra/var_hg19_evs.tsv.opt $dest/hg19/extra/var_hg19_evs.tsv.info} -vars {dest build db} -code {
	cg downloaddb $dest/tmp hg19 evs http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.protein-hgvs-update.snps_indels.vcf.tar.gz
}

# exac
job reg_hg19_exac -targets {$dest/hg19/extra/var_hg19_exac.tsv $dest/hg19/extra/var_hg19_exac.tsv.opt $dest/hg19/extra/var_hg19_exac.tsv.info} -vars {dest build db} -code {
	cg downloaddb $dest/tmp hg19 exac ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/ExAC.r0.3.sites.vep.vcf.gz
}

# CADD
job reg_hg19_cadd -targets {$dest/hg19/extra/var_hg19_cadd.tsv $dest/hg19/extra/var_hg19_cadd.tsv.opt $dest/hg19/extra/var_hg19_cadd.info} -vars {dest build db} -code {
	wgetfile http://krishna.gs.washington.edu/download/CADD/v1.2/whole_genome_SNVs.tsv.gz ${dest}/tmp/hg19/whole_genome_SNVs.tsv.gz
	wgetfile http://cadd.gs.washington.edu/home ${dest}/tmp/hg19/var_hg19_cadd.tsv.info
	file rename -force ${dest}/tmp/hg19/var_hg19_cadd.tsv.info ${dest}/$build/extra/var_hg19_cadd.tsv.info
	file_write $dest/$build/extra/var_${build}_cadd.tsv.opt "fields\t{score pscore}\n"
	exec cg select -hc 1 -f {chromosome=$Chrom {begin=$Pos - 1} end=$Pos type="snp" ref=$Ref alt=$Alt score=$RawScore pscore=$PHRED} ${dest}/tmp/hg19/whole_genome_SNVs.tsv.gz | cg collapsealleles | lz4c - > $dest/tmp/hg19/var_hg19_cadd.tsv.lz4.temp
	file rename -force $dest/tmp/hg19/var_hg19_cadd.tsv.lz4.temp $dest/hg19/extra/var_hg19_cadd.tsv.lz4
}

# GERP
job GERP -targets {extra/reg_${build}_GERP.tsv extra/reg_${build}_GERP.info} -vars {dest build tables} -code {
	cd ${dest}/tmp/${build}
	set table allHg19RS_BW
	cg downloaddbinfo ${dest}/tmp ${build} $table
	file rename -force reg_${build}_$table.info ${dest}/${build}/extra/reg_${build}_GERP.info
	cg downloaddb ${dest}/tmp ${build} $table
	cg ucscwb2reg -p 1 -f {} ucsc_${build}_$table.tsv
	cg select -s - reg_ucsc_${build}_$table.tsv ${dest}/${build}/extra/reg_${build}_GERP.tsv.temp
	file rename -force ${dest}/${build}/extra/reg_${build}_GERP.tsv.temp ${dest}/${build}/extra/reg_${build}_GERP.tsv
	cg lz4 ${dest}/${build}/extra/reg_${build}_GERP.tsv
#	file rename -force reg_ucsc_${build}_$table.tsv ${dest}/${build}/extra/reg_${build}_GERP.tsv

}

# encode
foreach {jobname resultname infosrc tables} {
	enc_transcription wgEncodeCaltechRnaSeq wgEncodeCaltechRnaSeq {wgEncodeRegTxnCaltechRnaSeqGm12878R2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqH1hescR2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqHelas3R2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqHepg2R2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqHsmmR2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqHuvecR2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqK562R2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqNhekR2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqNhlfR2x75Il200SigPooled}
	enc_H3K27Ac wgEncodeH3k27ac wgEncodeBroadHistoneViewSignal {wgEncodeBroadHistoneGm12878H3k27acStdSig wgEncodeBroadHistoneH1hescH3k27acStdSig wgEncodeBroadHistoneHsmmH3k27acStdSig wgEncodeBroadHistoneHuvecH3k27acStdSig wgEncodeBroadHistoneK562H3k27acStdSig wgEncodeBroadHistoneNhekH3k27acStdSig wgEncodeBroadHistoneNhlfH3k27acStdSig}
	enc_H3K4Me1 wgEncodeH3k4me1 wgEncodeRegMarkH3k4me1 {wgEncodeBroadHistoneGm12878H3k4me1StdSig wgEncodeBroadHistoneH1hescH3k4me1StdSig wgEncodeBroadHistoneHsmmH3k4me1StdSig wgEncodeBroadHistoneHuvecH3k4me1StdSig wgEncodeBroadHistoneK562H3k4me1StdSig wgEncodeBroadHistoneNhekH3k4me1StdSig wgEncodeBroadHistoneNhlfH3k4me1StdSig}
	enc_H3K4Me3	wgEncodeH3k4me3 wgEncodeRegMarkH3k4me3   {wgEncodeBroadHistoneGm12878H3k4me3StdSig wgEncodeBroadHistoneH1hescH3k4me3StdSig wgEncodeBroadHistoneHsmmH3k4me3StdSig wgEncodeBroadHistoneHuvecH3k4me3StdSig wgEncodeBroadHistoneK562H3k4me3StdSig wgEncodeBroadHistoneNhekH3k4me3StdSig wgEncodeBroadHistoneNhlfH3k4me3StdSig}
} {
	cd $dest/$build
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
		file mkdir ${dest}/tmp/$build
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
foreach file [glob -nocomplain ../hg19-local/*] {
	catch {
		file delete extra/[file tail $file]
		cplinked $file [file tail $file]
	}
}

# extra dir
# targetseq exome
job reg_exome_targetseq -targets {extra/reg_hg19_exome_targetseq.tsv} -code {
	cd extra
	wgetfile http://tools.invitrogen.com/content/sfs/manuals/TargetSeq_exome_named_targets_hg19.bed
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
	set url ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFPv3.1a.zip
	set file [file tail $url]
	if {![file exists $file]} {
		wgetfile $url
	}
	exec unzip -o $file >@ stdout 2>@ stderr
	set files [glob dbNSFP*_variant.chr*]
	set f [open [lindex $files 0]]
	set header [split [string range [gets $f] 1 end] \t]
	close $f
	set header [list_change $header {
		hg19_chr chromosome hg19_pos(1-based) end
		chr hg38_chr pos(1-coor) hg38_pos
		hg18_pos(1-coor) hg18_pos {SLR_test_statistic } SLR_test_statistic
		GERP++_NR GERP_NR GERP++_RS GERP_RS 
		ESP6500_AA_AF ESP_AA_AF ESP6500_EA_AF ESP_EA_AF
	}]
	file_write pre_var_dbnsfp.tsv.temp [join $header \t]\n
	exec cat {*}[lsort -dict $files] | grep -v ^# >> pre_var_dbnsfp.tsv.temp
	file rename -force pre_var_dbnsfp.tsv.temp pre_var_dbnsfp.tsv
}

job var_hg19_dbnsfp -deps {${dest}/tmp/hg19/pre_var_hg19_dbnsfp} -targets {extra/var_hg19_dbnsfp.tsv extra/var_hg19_dbnsfp.tsv.opt} -vars {dest} -code {
	cd ${dest}/tmp/hg19
	cg select -s {chromosome end} -f {
		chromosome {begin=$end - 1} end {type="snp"} ref alt
		MetaSVM_score MetaSVM_pred
		MetaSVM_score MetaSVM_pred MetaLR_score MetaLR_pred
		SIFT_score SIFT_pred
		Polyphen2_HDIV_score Polyphen2_HDIV_pred Polyphen2_HVAR_score Polyphen2_HVAR_pred 
		LRT_score LRT_pred 
		MutationTaster_score MutationTaster_pred 
		FATHMM_score FATHMM_pred
		GERP_NR GERP_RS phyloP7way_vertebrate SiPhy_29way_pi SiPhy_29way_logOdds LRT_Omega 
	} pre_var_dbnsfp.tsv pre2_var_hg19_dbnsfp.tsv.temp
	file rename -force pre2_var_hg19_dbnsfp.tsv.temp pre2_var_hg19_dbnsfp.tsv
	cg groupby {chromosome begin end type} pre2_var_hg19_dbnsfp.tsv pre3_var_hg19_dbnsfp.tsv.temp
	file rename -force pre3_var_hg19_dbnsfp.tsv.temp pre3_var_hg19_dbnsfp.tsv
	cg select -f {
		chromosome begin end type {ref=lindex($ref,0)} alt 
		MetaSVM_score MetaSVM_pred MetaLR_score MetaLR_pred
		SIFT_score SIFT_pred
		Polyphen2_HDIV_score Polyphen2_HDIV_pred Polyphen2_HVAR_score Polyphen2_HVAR_pred 
		LRT_score LRT_pred 
		MutationTaster_score MutationTaster_pred 
		FATHMM_score FATHMM_pred
		GERP_NR GERP_RS MetaLR_score MetaLR_pred SiPhy_29way_pi SiPhy_29way_logOdds LRT_Omega 
	} pre3_var_hg19_dbnsfp.tsv var_hg19_dbnsfp.tsv.temp
	# move dbNSFPzip files to target
	file rename -force var_hg19_dbnsfp.tsv.temp ${dest}/hg19/extra/var_hg19_dbnsfp.tsv
	file_write ${dest}/hg19/extra/var_hg19_dbnsfp.tsv.opt "fields\t{SIFT_score Polyphen2_HDIV_score Polyphen2_HDIV_pred Polyphen2_HVAR_score Polyphen2_HVAR_pred MetaSVM_score MetaSVM_pred MetaLR_score MetaLR_pred LRT_score LRT_pred MutationTaster_score MutationTaster_pred FATHMM_score GERP_NR GERP_RS MetaLR_score MetaLR_pred SiPhy_29way_pi SiPhy_29way_logOdds LRT_Omega ESP_AA_AF ESP_EA_AF}"
	file copy -force [glob dbNSFP*.readme.txt] ${dest}/hg19/extra/var_hg19_dbnsfp.info
}

job extragenome -deps {genome_${build}.ifas genome_${build}.ifas.index genome_${build}.ssa} -vars build \
-targets {extra/genome_${build}.ifas extra/genome_${build}.ifas.index extra/genome_${build}.ssa} -code {
	mklink genome_${build}.ifas extra/genome_${build}.ifas
	mklink genome_${build}.ifas.index extra/genome_${build}.ifas.index 
	mklink genome_${build}.ssa extra/genome_${build}.ssa
}
# genome in extra
foreach file [glob genome_*] {
	catch {
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
