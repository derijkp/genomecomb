#!/bin/sh
# the next line restarts using tclsh \
exec cg source "$0" "$@"

set build hg18
if {![info exists argv]} {set argv {}}
set argv [job_init {*}$argv]
if {[llength $argv]} {
	set dest [lindex $argv 0]
} else {
	set dest /complgen/refseq/
}

# download hg18
# =============
#
file mkdir ${dest}/${build}
cd ${dest}/${build}
file mkdir extra
file mkdir ${dest}/tmp/hg18

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
	cytoBand evofold gwasCatalog microsat oreganno rmsk simpleRepeat targetScanS tfbsConsSites 
	tRNAs wgRna vistaEnhancers gad
	phastConsElements28way phastConsElements28wayPlacMammal phastConsElements44way
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

# direct
foreach db {
	firstEF rnaGene
} {
	job reg_${build}_$db -targets {reg_${build}_${db}.tsv} -vars {dest build db} -code {
		cg downloaddb $dest ${build} $db
		file rename -force ucsc_${build}_${db}.tsv reg_${build}_${db}.tsv
	}
}

job reg_${build}_gwasCatalog -vars {build dest} -deps {ucsc_${build}_gwasCatalog.tsv} -targets {reg_${build}_gwasCatalog.tsv} -code {
	cg downloaddb $dest ${build} gwasCatalog
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
#job 1000g -targets {$dest/hg18/var_hg18_1000gCHBxJPT.tsv $dest/hg18/var_hg18_1000gCEU.tsv $dest/hg18/var_hg18_1000gYRI.tsv} -vars {dest} -code {
#	cd ${dest}/hg18
#	cg downloaddb ${dest} hg18 1000g
#}
#
## 1000glow (hg19)
#file mkdir ${dest}/hg19
#job 1000glow -targets {${dest}/hg19/var_hg19_1000glow.tsv} -vars {dest} -code {
#	cd $dest/hg19
#	cg downloaddb ${dest} hg19 1000glow
#}
#
#job 1000glow_liftover -deps {${dest}/hg19/var_hg19_1000glow.tsv} -targets {var_${build}_1000glow.tsv} -vars {dest build} -code {
#	cg liftover ${dest}/hg19/var_hg19_1000glow.tsv var_${build}_1000glow.tsv ${dest}/liftover/hg19ToHg18.over.chain
#}

job 1000g3 -targets {${dest}/hg19/var_hg19_1000g3.tsv ${dest}/hg19/extra/var_hg19_1000g3.tsv.opt} -vars {dest build} -code {
	cd $dest/hg19
	cg downloaddb ${dest} hg19 1000g3
	cplinked $target $dest/$build/extra/var_${build}_1000g3.tsv
	file_write $dest/$build/extra/var_${build}_1000g3.tsv.opt "fields\t{EUR_AF AMR_AF EAS_AF SAS_AF AFR_AF}\n"
}

job 1000g3_liftover -deps {${dest}/hg19/var_hg19_1000g3.tsv} -targets {var_${build}_1000g3.tsv} -vars {dest build} -code {
	cg liftover ${dest}/hg19/var_hg19_1000g3.tsv var_${build}_1000g3.tsv ${dest}/liftover/hg19ToHg18.over.chain
}

# dbsnp
job dbsnp130 -targets {var_${build}_snp130.tsv} -vars {dest build} -code {
	cg downloaddb ${dest} ${build} snp130
}

# dbsnp144 (hg19)
job dbsnp144 -targets {${dest}/hg19/var_hg19_snp144.tsv} -vars {dest} -code {
	cd $dest/hg19
	cg downloaddb ${dest} hg19 snp144
}

# dbsnp144Common (hg19)
job dbsnp144Common -targets {${dest}/hg19/var_hg19_snp144Common.tsv} -vars {dest build} -code {
	cd $dest/hg19
	cg downloaddb ${dest} hg19 snp144Common
}

# liftover
job dbsnp144lift -deps {${dest}/hg19/var_hg19_snp144.tsv} -targets {var_${build}_snp144lift.tsv} -vars {dest build} -code {
	cg liftover -split 0 ${dest}/hg19/var_hg19_snp144.tsv ${dest}/${build}/var_${build}_snp144lift.tsv ${dest}/liftover/hg19ToHg18.over.chain
	file delete ${dest}/${build}/var_${build}_snp144lift.tsv.unmapped
}

job dbsnp144Commonlift -deps {${dest}/hg19/var_hg19_snp144Common.tsv} -targets {var_${build}_snp144Commonlift.tsv} -vars {dest build} -code {
	cg liftover ${dest}/hg19/var_hg19_snp144Common.tsv ${dest}/${build}/var_${build}_snp144Commonlift.tsv ${dest}/liftover/hg19ToHg18.over.chain
	file delete ${dest}/${build}/var_${build}_snp144Commonlift.tsv.unmapped
}

foreach db {
	snp130 snp144lift snp144Commonlift
} {
	job maketabix_${build}_$db -deps {reg_${build}_${db}.tsv} -targets {reg_${build}_${db}.tsv.gz.tbi reg_${build}_${db}.tsv.gz} -vars {dest build db} -code {
		cg maketabix reg_${build}_${db}.tsv
		# exec gunzip -c var_${build}_${db}.tsv.gz > var_${build}_${db}.tsv
		cg select -f 'chrom start end type ref alt name freq' ${dest}/${build}/var_${build}_${db}.tsv.gz ${dest}/${build}/var_${build}_${db}.tsv
	}
}

# genes
foreach db {
	refGene ensGene knownGene wgEncodeGencodeManualV3 wgEncodeGencodeAutoV3 genscan augustusAbinitio acembly
} {
	job gene_${build}_$db -targets {gene_${build}_${db}.tsv gene_${build}_${db}.tsv.gz.tbi gene_${build}_${db}.tsv.gz} -vars {dest build db} -code {
		if {$db eq "wgEncodeGencodeCompV19"} {set dbname gencode} else {set dbname $db}
		if {$db in "genscan acembly augustusAbinitio"} {set geneidcol name} else {set geneidcol name2}
	        cg downloaddb ${dest} ${build} $db
		if {$db in "ensGene knownGene"} {
			unset -nocomplain a
			if {$db eq "ensGene"} {
				cg downloaddb ${dest} ${build} knownToEnsembl
				cg downloaddb ${dest} ${build} kgXref
				unset -nocomplain kga
				array set kga [split [string trim [cg select -f {kgID geneSymbol} ucsc_${build}_kgXref.tsv]] "\n\t"]
				list_foreach {kg ens} [split [string trim [file_read ucsc_hg18_knownToEnsembl.tsv]] \n] {
					set a($ens) [get kga($kg) $ens]
				}
				unset -nocomplain kga
			} else {
				cg downloaddb ${dest} ${build} kgXref
				array set a [split [string trim [cg select -f {kgID geneSymbol} ucsc_${build}_kgXref.tsv]] "\n\t"]
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
			file delete gene_${build}_${dbname}.tsv.temp2 ucsc_hg18_ensemblToGeneName.tsv
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

job reg_${build}_genes -targets {extra/reg_${build}_genes.tsv} \
-deps {gene_${build}_refGene.tsv gene_${build}_ensGene.tsv gene_${build}_knownGene.tsv gene_${build}_wgEncodeGencodeAutoV3.tsv gene_${build}_wgEncodeGencodeManualV3.tsv} \
-code {
	exec cg cat -fields {chrom start end geneid} {*}$deps | cg select -s {chrom start end geneid} | cg regcollapse > $target.temp
	file rename -force $target.temp $target
}

# homopolymer
job reg_${build}_homopolymer -deps {genome_${build}.ifas} -targets {reg_${build}_homopolymer.tsv reg_${build}_homopolymer.tsv.gz reg_${build}_homopolymer.tsv.gz.tbi reg_${build}_homopolymer.tsv.opt} -vars {dest build db} -code {
	cg extracthomopolymers genome_${build}.ifas > reg_${build}_homopolymer.tsv.temp
	file rename -force reg_${build}_homopolymer.tsv.temp reg_${build}_homopolymer.tsv
        cg maketabix reg_${build}_homopolymer.tsv
        gunzip reg_${build}_homopolymer.tsv.gz reg_${build}_homopolymer.tsv
	file_write reg_${build}_homopolymer.tsv.opt "fields\t{base size}\n"
}

# mirbase (hg19)
job reg_hg19_mirbase -targets {$dest/hg19/reg_hg19_mirbase.tsv $dest/hg19/reg_hg19_mirbase.tsv.opt $dest/hg19/reg_hg19_mirbase.info} -vars {dest build db} -code {
	cd $dest/hg19
	file_write $dest/hg19/reg_hg19_mirbase.tsv.opt "fields\t{ID}\n"
	exec -ignorestderr wget -c --tries=45 --directory-prefix=${dest}/tmp/hg19 ftp://mirbase.org/pub/mirbase/20/genomes/hsa.gff2
	cg gff2sft ${dest}/tmp/hg19/hsa.gff2 ${dest}/tmp/hg19/reg_hg19_mirbase.tsv.temp
	cg select -s - ${dest}/tmp/hg19/reg_hg19_mirbase.tsv.temp ${dest}/tmp/hg19/reg_hg19_mirbase.tsv.temp2
	file rename -force ${dest}/tmp/hg19/reg_hg19_mirbase.tsv.temp2 reg_hg19_mirbase.tsv
	exec -ignorestderr wget -c ftp://mirbase.org/pub/mirbase/20/README
	file rename -force README reg_hg19_mirbase.info
}

# mirbase hg18 liftover
job reg_hg18_mirbase_liftover -deps {${dest}/hg19/reg_hg19_mirbase.tsv ${dest}/hg19/reg_hg19_mirbase.info} \
-targets {reg_${build}_mirbase.tsv reg_${build}_mirbase.tsv.opt reg_${build}_mirbase.info} -vars {dest build db} -code {
	file_write reg_${build}_mirbase.tsv.opt "fields\t{ID}\n"
	cg liftover ${dest}/hg19/reg_hg19_mirbase.tsv ${dest}/tmp/${build}/reg_${build}_mirbase.tsv.temp2 ${dest}/liftover/hg19ToHg18.over.chain
	file rename -force ${dest}/tmp/${build}/reg_${build}_mirbase.tsv.temp2 reg_${build}_mirbase.tsv
	file copy ${dest}/hg19/reg_hg19_mirbase.info reg_${build}_mirbase.info
}

# encode
# enc_transcription
job enc_transcription -targets {reg_${build}_wgEncodelogRnaSeq.tsv} -vars {dest build} -code {
	cd ${dest}/tmp/${build}
	set tables {wgEncodeCaltechRnaSeqRawSignalRep1Gm12878CellLongpolyaBb12x75 wgEncodeCaltechRnaSeqRawSignalRep1H1hescCellPapBb2R2x75 wgEncodeCaltechRnaSeqRawSignalRep2Hepg2CellPapBb2R2x75 wgEncodeCaltechRnaSeqRawSignalRep1HuvecCellPapBb2R2x75 wgEncodeCaltechRnaSeqRawSignalRep1K562CellLongpolyaBb12x75 wgEncodeCaltechRnaSeqRawSignalRep1NhekCellPapBb2R2x75}
	cg downloaddb ${dest}/tmp ${build} {*}$tables
	set todo {}
	foreach table $tables {
		cg ucscwiggle2reg -n 0.01 -p 1 -f {log10($value)} ucsc_${build}_$table.tsv
		lappend todo reg_ucsc_${build}_$table.tsv
	}
	cg collapseoverlap -o ${dest}/${build}/reg_${build}_wgEncodelogRnaSeq.tsv {*}$todo
}

# make enc_transcription info file
job enc_transcription_info -targets {reg_${build}_wgEncodeCaltechRnaSeq.info} -vars {dest build} -code {
	cd ${dest}/tmp/${build}
	cg downloaddbinfo ${dest}/tmp ${build} wgEncodeCaltechRnaSeq
	set f [open ${dest}/tmp/${build}/reg_${build}_wgEncodeCaltechRnaSeq.info]
	set o [open ${dest}/${build}/reg_${build}_wgEncodeCaltechRnaSeq.info.temp w]
	set line [gets $f]
	puts $o $line
	puts $o "\n\n== Agregation info ==\n\nThis file combines the data from 6 lines"
	puts $o "score contains the logarithm of the highest score, with a precision of 1 digit after the decimal point."
	puts $o "num contains the number of lines for which the original value is >= 0.01 (thus log >= -2.0)\n"
	fcopy $f $o
	close $f
	close $o
	file rename -force ${dest}/${build}/reg_${build}_wgEncodeCaltechRnaSeq.info.temp ${dest}/${build}/reg_${build}_wgEncodeCaltechRnaSeq.info
}

foreach {jobname resultname infosrc tables} {
	enc_H3K4Me1 wgEncodeH3k4me1 wgEncodeRegMarkEnhH3k4me1 {wgEncodeBroadChipSeqSignalGm12878H3k4me1 wgEncodeBroadChipSeqSignalH1hescH3k4me1 wgEncodeBroadChipSeqSignalHmecH3k4me1 wgEncodeBroadChipSeqSignalHsmmH3k4me1 wgEncodeBroadChipSeqSignalHuvecH3k4me1 wgEncodeBroadChipSeqSignalK562H3k4me1 wgEncodeBroadChipSeqSignalNhekH3k4me1 wgEncodeBroadChipSeqSignalNhlfH3k4me1}
	enc_H3K27Ac wgEncodeH3k27ac wgEncodeRegMarkEnhH3k27ac {wgEncodeBroadChipSeqSignalGm12878H3k27ac wgEncodeBroadChipSeqSignalHepg2H3k27ac wgEncodeBroadChipSeqSignalHmecH3k27ac wgEncodeBroadChipSeqSignalHsmmH3k27ac wgEncodeBroadChipSeqSignalHuvecH3k27ac wgEncodeBroadChipSeqSignalK562H3k27ac wgEncodeBroadChipSeqSignalNhekH3k27ac wgEncodeBroadChipSeqSignalNhlfH3k27ac}
	enc_H3K4Me3	wgEncodeH3k4me3 wgEncodeRegMarkPromoter   {wgEncodeBroadChipSeqSignalGm12878H3k4me3 wgEncodeBroadChipSeqSignalH1hescH3k4me3 wgEncodeBroadChipSeqSignalHepg2H3k4me3 wgEncodeBroadChipSeqSignalHmecH3k4me3 wgEncodeBroadChipSeqSignalHsmmH3k4me3 wgEncodeBroadChipSeqSignalHuvecH3k4me3 wgEncodeBroadChipSeqSignalK562H3k4me3 wgEncodeBroadChipSeqSignalNhekH3k4me3 wgEncodeBroadChipSeqSignalNhlfH3k4me3}
	
} {
	# make database
	job $jobname -targets {${dest}/${build}/reg_${build}_$resultname.tsv} -vars {dest build tables} -code {
		cd ${dest}/tmp/${build}
		cg downloaddb ${dest}/tmp ${build} {*}$tables
		set todo {}
		foreach table $tables {
			cg ucscwiggle2reg -n 10 -p 0 -f {5*(($value+4)/5)} ucsc_${build}_$table.tsv
			lappend todo reg_ucsc_${build}_$table.tsv
		}
		cg collapseoverlap -o $target {*}$todo
	}
	# make info file
	job ${jobname}_info -targets {${dest}/${build}/reg_${build}_$resultname.info} -vars {dest build infosrc tables} -code {
		cd ${dest}/tmp/${build}
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
job enc_RegDnaseClustered -targets {reg_${build}_wgEncodeRegDnaseClustered.tsv reg_${build}_wgEncodeRegDnaseClustered.info} -vars {dest build} -code {
	cd ${dest}/tmp/${build}
	cg downloaddb ${dest}/tmp ${build} wgEncodeRegDnaseClustered
	cg collapseoverlap ucsc_${build}_wgEncodeRegDnaseClustered.tsv
	file rename -force reg_${build}_wgEncodeRegDnaseClustered.info ${dest}/${build}
	file rename -force reg_${build}_wgEncodeRegDnaseClustered.tsv ${dest}/${build}
}

job enc_RegTfbsClustered -targets {reg_${build}_wgEncodeRegTfbsClustered.tsv reg_${build}_wgEncodeRegTfbsClustered.info} -vars {dest build} -code {
	cd ${dest}/tmp/${build}
	cg downloaddb ${dest}/tmp ${build} wgEncodeRegTfbsClustered
	cg collapseoverlap ucsc_${build}_wgEncodeRegTfbsClustered.tsv
	file rename -force reg_${build}_wgEncodeRegTfbsClustered.info ${dest}/${build}
	file rename -force reg_${build}_wgEncodeRegTfbsClustered.tsv ${dest}/${build}
}

# link local data in dir
catch {
	cplinked {*}[glob ../hg18-local/*] .
}

# extra

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

# move dbNSFP hg19 files to target hg18
job dbNSFP_hg18_liftover -deps {$dest/hg19/extra/var_hg19_dbnsfp.tsv} -targets {extra/var_${build}_dbnsfp.tsv extra/var_${build}_dbnsfp.tsv.opt extra/var_${build}_dbnsfp.info} -vars {dest build db} -code {
	file delete ${dest}/tmp/${build}/var_${build}_dbnsfp.tsv.temp
	cg liftover -split 0 $dest/hg19/extra/var_hg19_dbnsfp.tsv ${dest}/tmp/${build}/var_${build}_dbnsfp.tsv.temp ${dest}/liftover/hg18ToHg19.over.chain
	file copy -force $dest/hg19/extra/var_hg19_dbnsfp.tsv.opt ${dest}/${build}/extra/var_${build}_dbnsfp.tsv.opt
	file copy -force $dest/hg19/extra/var_hg19_dbnsfp.info ${dest}/${build}/extra/var_${build}_dbnsfp.info
	file rename -force ${dest}/tmp/${build}/var_${build}_dbnsfp.tsv.temp ${dest}/${build}/extra/var_${build}_dbnsfp.tsv
}

job_wait

# todo
# cd /complgen/backup2013-10/hg18/
# cp reg_hg18_conserved.tsv reg_hg18_consnoncoding.tsv reg_hg18_refgene.tsv /complgen/hg18/extra
