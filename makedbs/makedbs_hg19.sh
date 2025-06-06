#!/bin/sh
# the next line restarts using tclsh \
exec cg source "$0" "$@"

set build hg19
set defaultdest /complgen/refseqnew

# settings
set genomeurl {}
set par {chromosome	begin	end	name
X	60001	2699520	PAR1
X	154931044	155260560	PAR2
Y	10001	2649520	PAR1
Y	59034050	59363566	PAR2
}
set organelles {chromosome
chrM
}
set mirbasegenome hsa-20:hg19
set dbsnpversion 151
set 1000g3url ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz
set 1000g3readmeurl ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/README_phase3_callset_20150220
set 1000g3build hg19
set clinvarurl ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20200615.vcf.gz
set clinvarpapuurl ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20200615_papu.vcf.gz
set kaviarurl http://s3-us-west-2.amazonaws.com/kaviar-160204-public/Kaviar-160204-Public-hg19-trim.vcf.tar
set kaviarbuild hg19
set evsurl http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.protein-hgvs-update.snps_indels.vcf.tar.gz
set evsbuild hg19
set exacurl ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/ExAC.r1.sites.vep.vcf.gz
set exacbuild hg19
set caddversion 1.6
set caddurl http://krishna.gs.washington.edu/download/CADD/v$caddversion/GRCh37/whole_genome_SNVs.tsv.gz
set caddbuild hg19
set gnomadbuild hg19
set gnomadversion 2.1.1
set gnomadbaseurl https://storage.googleapis.com/gnomad-public/release/$gnomadversion/vcf
set gnomadexbuild hg19
set gnomadexversion 2.1.1
set gnomadexurl https://storage.googleapis.com/gnomad-public/release/$gnomadexversion/vcf/exomes/gnomad.exomes.r$gnomadexversion.sites.vcf.bgz
set gnomadlofbuild hg19
set gnomadlofversion 2.1.1
set gnomadlof https://storage.googleapis.com/gnomad-public/release/$gnomadlofversion/constraint/gnomad.v$gnomadlofversion.lof_metrics.by_gene.txt.bgz
set ccrversion 2.20180420
set ccrurl https://s3.us-east-2.amazonaws.com/ccrs/ccrs/ccrs.autosomes.v${ccrversion}.bed.gz
set ccrbuild hg19
set dbnsfpurl ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFP4.0a.zip
set dbnsfpbuild hg38
set regionsdb_collapse {
	cytoBand evofold microsat oreganno rmsk simpleRepeat targetScanS tfbsConsSites 
	tRNAs wgRna vistaEnhancers gad
	phastConsElements46way phastConsElements46wayPlacental phastConsElements46wayPrimates
	phastConsElements100way wgEncodeDacMapabilityConsensusExcludable
}
set refSeqFuncElemsurl https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/105.20201022/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz
set regionsdb_join {
	chainSelf dgvMerged genomicSuperDups
}
set gencodeversion 34
set genesdb [list \
	{refGene int reg} \
	[list wgEncodeGencodeBasicV${gencodeversion} gencode extra int reg] \
	{ensGene extra int} \
	{knownGene extra int reg} \
	[list wgEncodeGencodeCompV${gencodeversion} cgencode extra] \
	{genscan extra} \
	{acembly extra} \
	{augustusGene extra} \
	{lincRNAsTranscripts lincRNA} \
]
set gtfurl https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.ensGene.gtf.gz
set gtffile genes_hg19_ensGene.gtf.gz

# keep actual command line used for log
set cmdline [clean_cmdline cg [info script] {*}$args]

# arguments
if {![info exists argv]} {set argv {}}
set argv [job_init {*}$argv]
foreach {dest webcache} $argv break
if {![info exists dest]} {set dest $defaultdest}
if {![info exists webcache]} {set webcache $dest/webcache}
if {[info exists webcache]} {set env(webcache) $webcache}
set dest [file_absolute $dest]

# prepare
putslog "Installing in $dest/$build"
file mkdir ${dest}/${build}
cd ${dest}/${build}
file mkdir extra

logverbose 2
set_job_logdir log_jobs
job_logfile ${dest}/${build}/log_makedbs_${build} ${dest}/${build} $cmdline

# download
# ========
#

# makerefdb
# ---------
# first part is done with the more generic cg makerefdb, but we call it with makerefdb_job to run its jobs under the already started job manager
makerefdb_job \
	-genomeurl $genomeurl \
	-refSeqFuncElemsurl $refSeqFuncElemsurl \
	-genesdb $genesdb \
	-pseudoautosomal $par \
	-organelles $organelles \
	-dbsnp $dbsnpversion \
	-mirbase $mirbase \
	/complgen/refseqnew/$build >@ stdout 2>@ stderr

# rest after this is human specific
# ---------------------------------

job gtf -targets {
	$gtffile
} -vars {gtfurl gtffile} -code {
	wgetfile $gtfurl $gtffile.temp
	file rename $gtffile.temp $gtffile
}

job reg_${build}_gwasCatalog -targets {
	reg_${build}_gwasCatalog.tsv
} -vars {build dest} -code {
	set target [gzroot $target].zst
	cg download_ucsc $target.ucsc ${build} gwasCatalog
	cg select \
		-f {chrom start end trait pValue pubMedID name bin author pubDate journal title initSample replSample region genes riskAllele riskAlFreq pValueDesc orOrBeta ci95 platform cnv} \
		-nh {chrom start end name score pubMedID dbsnp bin author pubDate journal title initSample replSample region genes riskAllele riskAlFreq pValueDesc orOrBeta ci95 platform cnv} \
		$target.ucsc $target.temp
	cg regcollapse $target.temp > $target.temp2
	file delete $target.temp
	file delete $target.ucsc
	file rename -force -- $target.ucsc.info [gzroot $target].info
	compress $target.temp2 $target
}

job sniffles_trf -optional 1 -targets {
	extra/sniffles_hg19.trf.bed
} -vars {} -code {
	wgetfile https://github.com/fritzsedlazeck/Sniffles/raw/master/annotations/human_hs37d5.trf.bed extra/sniffles_hg19.trf.bed.temp
	exec cg bed2tsv extra/sniffles_hg19.trf.bed.temp \
		| cg select -f {{chromosome=if($chromosome regexp "^chr",$chromosome,"chr$chromosome")} begin end} \
		| cg select -s chromosome \
		| cg tsv2bed > extra/sniffles_hg19.trf.bed.temp2
	file rename -force extra/sniffles_hg19.trf.bed.temp2 extra/sniffles_hg19.trf.bed
	file delete extra/sniffles_hg19.trf.bed.temp
}

## 1000 genomes
job 1000g3 -targets {
	var_${build}_1000g3.tsv
	extra/var_${build}_1000g3.tsv
	extra/var_${build}_1000g3.tsv.opt
} -vars {dest 1000g3url 1000g3readmeurl 1000g3build build} -code {
	set target [gzroot $target].zst
	if {$1000g3build eq $build} {
		cg download_1000g3 $target $1000g3url $1000g3readmeurl
	} else {
		cg download_1000g3 $target.$1000g3build.zst $1000g3url $1000g3readmeurl
		liftover_refdb $target.$1000g3build.zst $target $dest $1000g3build $build
	}
	cg zindex $target
	cplinked $target extra/var_${build}_1000g3.tsv.zst
	cplinked $target.zsti extra/var_${build}_1000g3.tsv.zst.zsti
	file_write extra/var_${build}_1000g3.tsv.opt "fields\t{EUR_AF AMR_AF EAS_AF SAS_AF AFR_AF}\n"
}

job clinvar -targets {
	var_${build}_clinvar.tsv
} -vars {dest build clinvarurl clinvarpapuurl} -code {
	set target [gzroot $target].zst
	cg download_clinvar --stack 1 $target $build $clinvarurl $clinvarpapuurl
	cg zindex $target
}

if {[info exists kaviarurl]} {
	job kaviar -targets {
		var_${build}_kaviar.tsv
	} -vars {dest build kaviarurl kaviarbuild} -code {
		set target [gzroot $target].zst
		if {$kaviarbuild eq $build} {
			cg download_kaviar $target $kaviarurl
		} else {
			cg download_kaviar $target.$kaviarbuild.zst $kaviarurl
			liftover_refdb $target.$kaviarbuild.zst $target $dest $kaviarbuild $build
		}
	}
}

job reg_${build}_phenotype -deps {
	extra/reg_${build}_genes.tsv
} -targets {
	extra/reg_${build}_phenotype.tsv
	extra/geneannot_${build}_phenotype.tsv
} -vars {dest build} -code {
	set target [gzroot $target].zst
	# get target2, uses biomart (ensembl geneset) and clinvar for gene-phenotype correlations
	cg download_phenotype $target2 ${build}
	# to reg file
	cg geneannot2reg $dep $target2 $target.temp
	compress $target.temp $target
}

job reg_${build}_go -deps {
	extra/reg_${build}_genes.tsv
} -targets {
	extra/reg_${build}_go.tsv
	extra/geneannot_${build}_go.tsv
} -code {
	set target [gzroot $target].zst
	cg download_mart $target2.zst hsapiens_gene_ensembl gene_ensembl_config {hgnc_symbol name_1006 namespace_1003}
	cg geneannot2reg $dep $target2.zst $target.temp
	compress $target.temp $target
	file delete temp_go.tsv.temp
}

# biograph is no longer accessible
#foreach {name id disease pattern} {
#	ad	C0002395	{Alzheimer's Disease} Alzheimer
#	ftd	C0338451 {Frontotemporal Dementia} frontotemp
#	als	C0002736	{Amyotrophic Lateral Sclerosis} {amyotrophic lateral}
#	pd	C0030567	{Parkinson Disease} parkinson
#	ep	C0014544	{Epilepsy} epilep
#	cmt	C0007959	{Charcot-Marie-Tooth Disease} charcot
#	pn	C0031117	{Peripheral Neuropathy} {peripher.*neuropat|motor.*neuropat|neuropat.*motor|peripher.*nerv}
#} {
#	job reg_${build}_biograph_$name -vars {id pattern name disease} \
#	-deps {extra/reg_${build}_genes.tsv extra/geneannot_${build}_phenotype.tsv} \
#	-targets {extra/reg_${build}_biograph_$name.tsv extra/geneannot_${build}_biograph_$name.tsv extra/reg_${build}_biograph_$name.info} -code {
#		set genefile [tempfile]
#		cg select -q "\$phenotype_description regexp \"(?i)$pattern\"" $dep2 $genefile
#		cg download_biograph $target2.zst $id $genefile
#		cg geneannot2reg $dep $target2.zst $target.temp
#		cg zst -i 1 $target.temp
#		file rename -force -- $target.temp.zst [gzroot $target].zst
#		file rename -force -- $target.temp.zst.zsti [gzroot $target].zst.zsti
#		set temp "gene ranking\n============\ngene ranking for $disease\n\nThe following scores will be given to variants in genes that"
#		append temp "0: are directly linked to the disease in the ensembl gene database or in the clinvar database\n"
#		append temp "1: have a known relation to the disease in biograph (can be association, etc.)\n"
#		append temp ">1: position in genelist sorted according to biograph prioritization (there are 18182 genes)"
#		file_write $target3 $temp
#	}
#}

# exome variant server
job var_${build}_evs -targets {
	extra/var_${build}_evs.tsv
	extra/var_${build}_evs.tsv.opt
	extra/var_${build}_evs.tsv.info
} -vars {dest build db evsurl evsbuild} -code {
	set target [gzroot $target].zst
	if {$evsbuild eq $build} {
		cg download_evs $target $evsurl
	} else {
		cg download_evs $target.$evsbuild.zst $evsurl
		liftover_refdb $target.$evsbuild.zst $target $dest $evsbuild $build
	}
}

# exac
job var_${build}_exac -targets {
	extra/var_${build}_exac.tsv
	extra/var_${build}_exac.tsv.opt
	extra/var_${build}_exac.tsv.info
} -vars {dest build db exacurl exacbuild} -code {
	set target [gzroot $target].zst
	if {$exacbuild eq $build} {
		cg download_exac --stack 1 --verbose 2 $target $exacurl
	} else {
		file mkdir $target.temp
		cg download_exac --stack 1 --verbose 2 $target.temp/$exacbuild.tsv.zst $exacurl
		liftover_refdb $target.temp/$exacbuild.tsv.zst $target $dest $exacbuild $build 0
		file delete -force $target.temp
	}
}

# GERP
job GERP -targets {
	extra/bcol_${build}_GERP.bcol
	extra/bcol_${build}_GERP.bcol.info
} -vars {dest build tables} -code {
	set table allHg19RS_BW
	set tempdir $target.temp
	file mkdir $tempdir
	set tail [file tail $target]
	cg download_ucsc $tempdir/ucsc_${build}_$table.tsv ${build} $table
	cg ucscwb2reg -p 1 -f {} $tempdir/ucsc_${build}_$table.tsv $tempdir/reg_ucsc_${build}_$table.tsv
	cg select -s - $tempdir/reg_ucsc_${build}_$table.tsv $tempdir/reg_${build}_GERP.tsv.temp
	exec cg bcol make --precision 1 --compress 9 -t f -p begin -e end -c chromosome $tempdir/$tail score < $tempdir/reg_${build}_GERP.tsv.temp
	file rename -force -- $tempdir/$tail.bin.zst $target.bin.zst
	file rename -force -- $tempdir/$tail.bin.zst.zsti $target.bin.zst.zsti
	file rename -force -- $tempdir/$tail $target
	file rename -force -- $tempdir/ucsc_${build}_$table.tsv.info [gzroot $target].info
	file delete -force $tempdir
}

# encode
foreach {jobname resultname infosrc tables} {
	enc_transcription wgEncodeCaltechRnaSeq wgEncodeCaltechRnaSeq {wgEncodeRegTxnCaltechRnaSeqGm12878R2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqH1hescR2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqHelas3R2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqHepg2R2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqHsmmR2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqHuvecR2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqK562R2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqNhekR2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqNhlfR2x75Il200SigPooled}
	enc_H3K27Ac wgEncodeH3k27ac wgEncodeBroadHistone {wgEncodeBroadHistoneGm12878H3k27acStdSig wgEncodeBroadHistoneH1hescH3k27acStdSig wgEncodeBroadHistoneHsmmH3k27acStdSig wgEncodeBroadHistoneHuvecH3k27acStdSig wgEncodeBroadHistoneK562H3k27acStdSig wgEncodeBroadHistoneNhekH3k27acStdSig wgEncodeBroadHistoneNhlfH3k27acStdSig}
	enc_H3K4Me1 wgEncodeH3k4me1 wgEncodeRegMarkH3k4me1 {wgEncodeBroadHistoneGm12878H3k4me1StdSig wgEncodeBroadHistoneH1hescH3k4me1StdSig wgEncodeBroadHistoneHsmmH3k4me1StdSig wgEncodeBroadHistoneHuvecH3k4me1StdSig wgEncodeBroadHistoneK562H3k4me1StdSig wgEncodeBroadHistoneNhekH3k4me1StdSig wgEncodeBroadHistoneNhlfH3k4me1StdSig}
	enc_H3K4Me3	wgEncodeH3k4me3 wgEncodeRegMarkH3k4me3   {wgEncodeBroadHistoneGm12878H3k4me3StdSig wgEncodeBroadHistoneH1hescH3k4me3StdSig wgEncodeBroadHistoneHsmmH3k4me3StdSig wgEncodeBroadHistoneHuvecH3k4me3StdSig wgEncodeBroadHistoneK562H3k4me3StdSig wgEncodeBroadHistoneNhekH3k4me3StdSig wgEncodeBroadHistoneNhlfH3k4me3StdSig}
} {
	# make database
	job $jobname -targets {
		bcol_${build}_$resultname.bcol
	} -vars {dest build tables resultname} -code {
		set tempdir $target.temp
		file mkdir $tempdir
		set todo {}
		foreach table $tables {
			lappend todo $tempdir/reg_ucsc_${build}_$table.tsv
			if {[file exists $tempdir/reg_ucsc_${build}_$table.tsv]} {
				putslog "$tempdir/reg_ucsc_${build}_$table.tsv already there"
				continue
			}
			putslog "Downloading $tempdir/ucsc_$table.tsv"
			cg download_ucsc $tempdir/ucsc_$table.tsv $build $table 2>@ stderr
			cg ucscwb2reg -n 10 -p 0 -f {5*(($value+4)/5)} $tempdir/ucsc_$table.tsv $tempdir/reg_ucsc_${build}_$table.tsv 2>@ stderr
		}
		if {![file exists $tempdir/reg_${build}_$resultname.tsv]} {
			cg regcollapse -o $tempdir/reg_${build}_$resultname.tsv.temp {*}$todo
			file rename -force -- $tempdir/reg_${build}_$resultname.tsv.temp $tempdir/reg_${build}_$resultname.tsv
		}
		cg bcol make --compress 9 -t iu -p begin -e end -c chromosome $tempdir/bcol_${build}_$resultname.bcol score < $tempdir/reg_${build}_$resultname.tsv
		file rename -force -- $tempdir/bcol_${build}_$resultname.bcol.bin.zst bcol_${build}_$resultname.bcol.bin.zst
		file rename -force -- $tempdir/bcol_${build}_$resultname.bcol.bin.zst.zsti bcol_${build}_$resultname.bcol.bin.zst.zsti
		file rename -force -- $tempdir/bcol_${build}_$resultname.bcol bcol_${build}_$resultname.bcol
		file delete -force $tempdir
	}
	# make info file
	job ${jobname}_info -targets {bcol_${build}_$resultname.bcol.info} -vars {dest build infosrc tables resultname} -code {
		cg download_ucscinfo $target.temp ${build} $infosrc
		set c [file_read $target.temp]
		set c [string_change $c [list {== Description ==} [deindent [subst {
			== Agregation info ==
			This file combines the data from [llength $tables] cell lines ([join $tables .])
			score contains the highest score rounded up the next 5 fold.
			num contains the number of cell lines for which the score is >= 10
			
			== Description ==
		}]]]]
		file_write $target.temp2 $c
		file rename -force -- $target.temp2 $target
		file delete $target.temp
	}
}

# DNase Clusters and Txn Factor ChIP
job enc_RegDnaseClustered -targets {
	reg_${build}_wgEncodeRegDnaseClustered.tsv
	reg_${build}_wgEncodeRegDnaseClustered.tsv.info
} -vars {dest build} -code {
	cg download_ucsc $target.ucsc $build wgEncodeRegDnaseClusteredV3
	cg download_ucscinfo $target.ucsc.info $build wgEncodeRegDnaseClustered
	cg regcollapse $target.ucsc | cg zst > $target.temp.zst
	file rename -force -- $target.ucsc.info $target.info
	file rename -force -- $target.temp.zst $target.zst
	file delete $target.ucsc
	cg zstindex $target.zst
}

job enc_RegTfbsClustered -targets {
	reg_${build}_wgEncodeRegTfbsClusteredV3.tsv
} -vars {dest build} -code {
	cg download_ucsc $target.ucsc ${build} wgEncodeRegTfbsClusteredV3
	cg select -s - -f {chrom	start	end	name	score} $target.ucsc $target.temp
	cg regcollapse $target.temp | cg zst > $target.temp2.zst
	file rename -force -- $target.temp2.zst $target.zst
	file rename -force -- $target.ucsc.info $target.info
	file delete -force $target.ucsc $target.temp
	cg zstindex $target.zst
}

# link local data in dir
#foreach file [glob -nocomplain ../${build}-local/*] {
#	catch {
#		file delete extra/[file tail $file]
#		cplinked $file [file tail $file]
#	}
#}

# extra dir
cd ${dest}/${build}

# dbNSFP
job var_${build}_dbnsfp -targets {
	extra/var_${build}_dbnsfp.tsv
	extra/var_${build}_dbnsfp.tsv.opt
} -vars {dest build dbnsfpurl dbnsfpbuild} -code {
	set target [gzroot $target].zst
	cg download_dbnsfp $target $build $dbnsfpurl $dbnsfpbuild
}

# ccr
# ---
job ccr -deps {
} -vars {
	build ccrversion ccrurl ccrbuild dest
} -targets {
	reg_${build}_ccr.tsv
} -code {
	file_write $target.info [subst [deindent {
		= CCR (constrained coding regions) =
		
		== Download info ==
		dbname	ccr
		version	$ccrversion
		citation	Havrilla, J.M., Pedersen, B.S., Layer, R.M. & Quinlan, A.R. A map of constrained coding regions in the human genome. Nature Genetics (2018). doi:10.1038/s41588-018-0294-6
		license	cite
		source	$ccrurl
		time	[timestamp]
		
		== Description ==
		
		More info on https://github.com/quinlan-lab/ccrhtml
		and in https://www.nature.com/articles/s41588-018-0294-6
		
		== Category ==
		Annotation
	}]]
	file_write $target.opt "fields\t{ccr_pct}\n"
	file mkdir $target.temp
	set tail [file tail $ccrurl]
	wgetfile $ccrurl $target.temp/$tail
	if {$build ne $ccrbuild} {
		cg select -s - -overwrite 1 -hc 1 -f {chrom start end {ccr_pct=format("%.2f",$ccr_pct)} *} $target.temp/$tail $target.$ccrbuild.temp
		liftover_refdb $target.$ccrbuild.temp $target.zst $dest $ccrbuild $build
	} else {
		cg select -s - -overwrite 1 -hc 1 -f {chrom start end {ccr_pct=format("%.2f",$ccr_pct)} *} $target.temp/$tail $target.temp.zst
		file rename -force -- $target.temp.zst $target.zst
	}
	cg zstindex $target.zst
	file delete -force $target.temp
}

# gnomad lof
# ----------
job lofgnomad -deps {
} -vars {
	build gnomadlof gnomadlofbuild gnomadlofversion dest
} -targets {
	reg_${build}_lofgnomad.tsv
} -code {
	file_write $target.info [subst [deindent {
		= lofgnomad (loss of function tolerence) =
		
		== Download info ==
		dbname	lofgnomad
		version	$gnomadlofversion
		citation	Lek M., Karczewski K., Exome Aggregation Consortium. Analysis of protein-coding genetic variation in 60,706 humans. Nature volume 536, pages 285-291 (2016)
		license	cite
		source	$gnomadlof
		time	[timestamp]
		
		== Description ==
		
		The probability of being loss-of-function (LoF) intolerant (pLI) separates
		genes of sufficient length into LoF intolerant (pLI >= 0.9, n=3,230) or
		LoF tolerant (pLI <= 0.1, n=10,374) categories

		More info on https://gnomad.broadinstitute.org/faq
		and in https://www.nature.com/articles/nature19057
		
		== Category ==
		Annotation
	}]]
	file_write $target.opt "fields\tpLI\n"
	file mkdir $target.temp
	set tail [file tail $gnomadlof]
	wgetfile $gnomadlof $target.temp/$tail
	cg select -overwrite 1 -f {chromosome {begin=$start_position} {end=$end_position} pLI exac_pLI gene *} $target.temp/$tail $target.temp/temp.tsv
	cg select -overwrite 1 -s - $target.temp/temp.tsv $target.temp/temp2.tsv
	if {$build ne $gnomadlofbuild} {
		cg select -overwrite 1 -rf {start_position end_position} $target.temp/temp2.tsv $target.temp/temp3.tsv
		liftover_refdb $target.temp/temp3.tsv $target.zst $dest $gnomadlofbuild $build
	} else {
		cg select -overwrite 1 -rf {start_position end_position} $target.temp/temp2.tsv $target.temp2.zst
		file rename -force -- $target.temp2.zst $target.zst
	}
	cg zstindex $target.zst
	file delete -force $target.temp
}

# gnomad
# ------
set finaltarget var_${build}_gnomad.tsv.zst
set tempdir $finaltarget.temp
file mkdir $tempdir

job var_${build}_gnomad-info -targets {
	var_${build}_gnomad.tsv.info
} -vars {dest db build gnomadversion gnomadbaseurl gnomadbuild} -code {
	file_write var_${build}_gnomad.tsv.info [subst [deindent {
		= gnomAD (genome Aggregation Database) =
		
		== Download info ==
		dbname	gnomad
		version	$gnomadversion
		citation	Lek, M., Karczewski, K. J., Minikel, E. V., Samocha, K. E., Banks, E., Fennell, T., Exome Aggregation Consortium. (2016). Analysis of protein-coding genetic variation in 60,706 humans. Nature, 536(7616), 285-291. https://doi.org/10.1038/nature19057
		source	$gnomadbaseurl
		time	[timestamp]
		
		== Description ==
		The Genome Aggregation Database (gnomAD) is a resource developed by an
		international coalition of investigators, with the goal of aggregating and
		harmonizing both exome and genome sequencing data from a wide variety of
		large-scale sequencing projects, and making summary data available for the
		wider scientific community.
		
		The data set provided on this website spans 123,136 exome sequences and
		15,496 whole-genome sequences from unrelated individuals sequenced as part
		of various disease-specific and population genetic studies. The gnomAD
		Principal Investigators and groups that have contributed data to the
		current release are listed here.
		
		All data here are released for the benefit of the wider biomedical
		community, without restriction on use - see the terms of use here.

		More information on realeas 2.1 at https://macarthurlab.org/2018/10/17/gnomad-v2-1/

		== Category ==
		Annotation
	}]]
}

# multidownload does not work
# wgetfiles $url $tempdir

proc gnomadfields {file} {
	set tempfile [tempfile]
	catchchildkilled_exec zcat $file | head -10000 | cg vcf2tsv -split 1 > $tempfile
	set header [cg select -header $tempfile]
	set fields {chromosome begin end type ref alt}
	if {[llength [list_common $header {AF_popmax}]] == 1} {
		lappend fields {max_freqp=if(isnum($AF_popmax),format("%.3f",100.0*$AF_popmax),"-")}
	}
	# homfreqp is percentage of sindividuals that are homozygous -> so should divide count by (AN_popmax/2), do this by *200.0 instead of *100.0
	lappend fields {max_homfreqp=if(isnum($nhomalt_popmax),format("%.3f",200.0*$nhomalt_popmax/$AN_popmax),"-")}
	lappend fields {controls_max_freqp=if(isnum($controls_AF_popmax),format("%.3f",100.0*$controls_AF_popmax),"-")}
	lappend fields {controls_max_homfreqp=if(isnum($controls_nhomalt_popmax),format("%.3f",200.0*$controls_nhomalt_popmax/$controls_AN_popmax),"-")}
	lappend fields {non_neuro_max_freqp=if(isnum($non_neuro_AF_popmax),format("%.3f",100.0*$non_neuro_AF_popmax),"-")}
	lappend fields {non_neuro_max_homfreqp=if(isnum($non_neuro_nhomalt_popmax),format("%.3f",200.0*$non_neuro_nhomalt_popmax/$non_neuro_AN_popmax),"-")}
	lappend nh max_freqp
	foreach population {
		afr amr asj eas fin nfe sas oth nfe_nwe nfe_seu nfe_est nfe_bgr nfe_swe
		male female
	} {
		if {![inlist $header AN_$population]} continue
		lappend fields "${population}_freqp=if(def(\$AN_$population,0) < 8, \"-\", format(\"%.3f\",(100.0 * \$AC_$population)/\$AN_$population))"
		lappend fields "${population}_homfreqp=if(def(\$AN_$population,0) < 8, \"-\", format(\"%.3f\",(200.0 * \$nhomalt_$population)/\$AN_$population))"
		if {[inlist $header controls_AN_$population]} {
			lappend fields "controls_${population}_homfreqp=if(def(\$controls_AN_$population,0) < 8, \"-\", format(\"%.3f\",(200.0 * \$controls_nhomalt_$population)/\$controls_AN_$population))"
		}
		if {[inlist $header controls_AN_$population]} {
			lappend fields "controls_${population}_freqp=if(def(\$controls_AN_$population,0) < 8, \"-\", format(\"%.3f\",(100.0 * \$controls_AC_$population)/\$controls_AN_$population))"
		}
		if {[inlist $header controls_AN_${population}_male]} {
			lappend fields "controls_${population}_male_freqp=if(def(\$controls_AN_${population}_male,0) < 8, \"-\", format(\"%.3f\",(100.0 * \$controls_AC_${population}_male)/\$controls_AN_${population}_male))"
		}
		if {[inlist $header controls_AN_${population}_female]} {
			lappend fields "controls_${population}_female_freqp=if(def(\$controls_AN_${population}_female,0) < 8, \"-\", format(\"%.3f\",(100.0 * \$controls_AC_${population}_female)/\$controls_AN_${population}_female))"
		}
		if {[inlist $header non_neuro_AN_$population]} {
			lappend fields "non_neuro_${population}_freqp=if(def(\$non_neuro_AN_$population,0) < 8, \"-\", format(\"%.3f\",(100.0 * \$non_neuro_AC_$population)/\$non_neuro_AN_$population))"
		}
		if {[inlist $header faf95_${population}]} {
			lappend fields [subst {faf95_${population}_freqp=if(isnum(\$faf95_${population}), format(\"%.3f\",(100.0 * \$faf95_${population})),"-")}]
		}
		if {[inlist $header controls_faf95_${population}]} {
			lappend fields [subst {controls_faf95_${population}_freqp=if(isnum(\$controls_faf95_${population}), format(\"%.3f\",(100.0 * \$controls_faf95_${population})),"-")}]
		}
		if {[inlist $header non_neuro_faf95_${population}]} {
			lappend fields [subst {non_neuro_faf95_${population}_freqp=if(isnum(\$non_neuro_faf95_${population}), format(\"%.3f\",(100.0 * \$non_neuro_faf95_${population})),"-")}]
		}
		if {[inlist $header faf99_${population}]} {
			lappend fields [subst {faf99_${population}_freqp=if(isnum(\$faf99_${population}), format(\"%.3f\",(100.0 * \$faf99_${population})),"-")}]
		}
		if {[inlist $header controls_faf99_${population}]} {
			lappend fields [subst {controls_faf99_${population}_freqp=if(isnum(\$controls_faf99_${population}), format(\"%.3f\",(100.0 * \$controls_faf99_${population})),"-")}]
		}
		if {[inlist $header non_neuro_faf99_${population}]} {
			lappend fields [subst {non_neuro_faf99_${population}_freqp=if(isnum(\$non_neuro_faf99_${population}), format(\"%.3f\",(100.0 * \$non_neuro_faf99_${population})),"-")}]
		}
	}
	return $fields
}

set deps {}
foreach chromosome {
	1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
} {
	puts chr$chromosome
	set vcf gnomad.genomes.r$gnomadversion.sites.$chromosome.vcf.bgz
	set target $tempdir/result$chromosome.tsv.zst
	lappend deps $target
	job var_${build}_gnomad-$chromosome -targets {
		$target
	} -skip {
		$finaltarget
	} -vars {
		tempdir gnomadbaseurl vcf dest db build
	} -procs {
		gnomadfields
	} -code {
		if {![file exists $tempdir/$vcf]} {
			putslog "Downloading $vcf"
			wgetfile $gnomadbaseurl/genomes/$vcf $tempdir/$vcf
		}
		set fields [gnomadfields $tempdir/$vcf]
		putslog "Converting $vcf"
		cg vcf2tsv -split 1 -sort 0 $tempdir/$vcf | cg select --stack 1 -rc 1 -f $fields | cg collapsealleles {*}[compresspipe $target 1] > $target.temp
		file rename -force -- $target.temp $target
	}
}

job var_${build}_gnomad-final -deps $deps -targets {
	$finaltarget
	var_${build}_gnomad.tsv.opt
} -vars {tempdir dest db build finaltarget gnomadbuild} -code {
	file_write var_${build}_gnomad.tsv.opt "fields\t{max_freqp nfe_freqp}\n"
	if {$gnomadbuild eq $build} {
		exec cg cat {*}$deps {*}[compresspipe $finaltarget 12] > $tempdir/result.tsv.temp.zst
		file rename -force -- $tempdir/result.tsv.temp.zst $finaltarget
	} else {
		exec cg cat {*}$deps {*}[compresspipe $finaltarget 1] > $tempdir/result.tsv.temp.zst
		file rename -force -- [gzroot $finaltarget].info $tempdir/result.tsv.temp.info
		liftover_refdb $tempdir/result.tsv.temp.zst $finaltarget $dest $gnomadbuild $build 0
	}
	cg zindex $finaltarget
	cg index $finaltarget
	# file delete -force $tempdir
}	

set target var_${build}_gnomadex.tsv.zst
job var_${build}_gnomad_exomes -targets {
	$target
	var_${build}_gnomadex.tsv.opt
} -vars {
	gnomadexurl tempdir dest db build gnomadexbuild
} -procs {
	gnomadfields
} -code {
	file_write var_${build}_gnomadex.tsv.opt "fields\t{max_freqp nfe_freqp}\n"
	set tempdir $target.temp
	file mkdir $tempdir
	set vcf [file tail $gnomadexurl]
	wgetfile $gnomadexurl $tempdir/$vcf
	putslog "Converting $vcf"
	set fields [gnomadfields $tempdir/$vcf]
	cg vcf2tsv -split 1 -sort 0 $tempdir/$vcf | cg select --stack 1 -rc 1 -f $fields | cg collapsealleles {*}[compresspipe $target 12] > $tempdir/[file tail $target]
	if {$gnomadexbuild eq $build} {
		file rename -force -- $tempdir/[file tail $target] $target
	} else {
		file rename -force -- var_${build}_gnomadex.tsv.opt $tempdir/result.tsv.temp.info
		liftover_refdb $tempdir/[file tail $target] $target $dest $gnomadexbuild $build 0
	}
	cg zindex $target
	cg index $target
	file delete -force $tempdir
}

job var_${build}_extragnomad-final -deps {
	var_${build}_gnomad.tsv.zst
	var_${build}_gnomadex.tsv.zst
} -targets {
	extra/var_${build}_gnomad.tsv.zst
	extra/var_${build}_gnomadex.tsv.zst
	extra/var_${build}_gnomad.tsv.opt
	extra/var_${build}_gnomadex.tsv.opt
} -vars {build} -code {
	file_write extra/var_${build}_gnomad.tsv.opt "fields\t{max_homfreqp controls_max_homfreqp afr_freqp amr_freqp asj_freqp eas_freqp fin_freqp oth_freqp male_freqp female_freqp}\n"
	mklink var_${build}_gnomad.tsv.zst extra/var_${build}_gnomad.tsv.zst
	file_write extra/var_${build}_gnomadex.tsv.opt "fields\t{afr_freqp amr_freqp asj_freqp eas_freqp fin_freqp oth_freqp male_freqp female_freqp}\n"
	mklink var_${build}_gnomadex.tsv.zst extra/var_${build}_gnomadex.tsv.zst
	# file delete -force $tempdir
}

# CADD
job reg_${build}_cadd -targets {
	var_${build}_cadd.bcol
	var_${build}_cadd.bcol.bin.zst
	var_${build}_cadd.bcol.bin.zst.zsti
	var_${build}_cadd.bcol.info
} -vars {dest db build caddurl caddbuild caddversion} -code {
	set tempdir $target.temp
	file mkdir $tempdir
	set tail [file tail $caddurl]
	file_write var_${build}_cadd.bcol.info [subst [deindent {
		= CADD (Combined Annotation Dependent Depletion) =
		
		== Download info ==
		dbname	cadd
		version	$caddversion
		citation	Kircher M, Witten DM, Jain P, O'Roak BJ, Cooper GM, Shendure J. A general framework for estimating the relative pathogenicity of human genetic variants. Nat Genet. 2014 Feb 2. doi:10.1038/ng.2892 PubMed PMID:24487276
		license	non-commercial
		source	$caddurl
		time	[timestamp]
		
		== Description ==
		CADD is a tool for scoring the deleteriousness of single nucleotide variants as well as
		insertion/deletions variants in the human genome.
		While many variant annotation and scoring tools are around, most annotations 
		tend to exploit a single information type (e.g. conservation) and/or 
		are restricted in scope (e.g. to missense changes). Thus, a broadly applicable metric 
		that objectively weights and integrates diverse information is needed. 
		Combined Annotation Dependent Depletion (CADD) is a framework that 
		integrates multiple annotations into one metric by contrasting variants that 
		survived natural selection with simulated mutations.
		C-scores strongly correlate with allelic 
		diversity, pathogenicity of both coding and non-coding variants, and 
		experimentally measured regulatory effects, and also highly rank causal variants 
		within individual genome sequences. 
		More info on the CADD scores can be found on \[\[http://cadd.gs.washington.edu/home\]\]
		
		== Category ==
		Annotation
	}]]
	if {![file exists $tempdir/$tail]} {
		putslog "Downloading $tail"
		wgetfile $caddurl $tempdir/$tail
	}
	putslog "make bcol"
	file_write var_${build}_cadd.tsv.opt "fields\t{score pscore}\n"
	if {$caddbuild ne $build} {
		if {![file exists $tempdir/collapsed${build}.tsv.zst]} {
			if {![file exists $tempdir/collapsed.tsv.zst]} {
				cg select --stack 1 -hc 1 -rc 1 -f {{chrom=$Chrom} {begin = $Pos - 1} {end=$Pos} {ref=$Ref} {alt=$Alt} {score=$PHRED}} $tempdir/$tail \
					| cg collapsealleles --stack 1 {*}[compresspipe .zst 1] > $tempdir/collapsed.tsv.zst.temp
				file rename -force -- $tempdir/collapsed.tsv.zst.temp $tempdir/collapsed.tsv.zst
			}
			cg liftover --stack 1 -split 0 -s 0 $tempdir/collapsed.tsv.zst ../liftover/${caddbuild}To${build}.over.tsv {*}[compresspipe .zst] > $tempdir/collapsed${build}.tsv.temp
			file rename -force -- $tempdir/collapsed${build}.tsv.temp $tempdir/collapsed${build}.tsv.zst
			file delete $tempdir/collapsed.tsv.zst
		}
		exec cg select -s - $tempdir/collapsed${build}.tsv.zst | cg bcol make --stack 1 --precision 3 --compress 9 -t f --multicol alt --multilist A,C,T,G -p begin -c chrom $tempdir/var_${build}_cadd.bcol score
	} else {
		cg select -hc 1 -rc 1 -f {{chrom=$Chrom} {begin = $Pos - 1} {ref=$Ref} {alt=$Alt} {score=$PHRED}} $tempdir/$tail \
			| cg collapsealleles | cg select -s {chrom begin} \
			| cg bcol make --precision 3 --compress 9 -t f --multicol alt --multilist A,C,T,G -p begin -c chrom $tempdir/var_${build}_cadd.bcol score
	}
	file rename -force -- $tempdir/var_${build}_cadd.bcol.bin.zst var_${build}_cadd.bcol.bin.zst
	file rename -force -- $tempdir/var_${build}_cadd.bcol.bin.zst.zsti var_${build}_cadd.bcol.bin.zst.zsti
	file rename -force -- $tempdir/var_${build}_cadd.bcol var_${build}_cadd.bcol
	file delete -force $tempdir
}

# geneHancerRegElements and geneHancerClusteredInteractions cannot be
# downloaded (yet?) from UCSC via ftp, did this manually (tablebrowser)
# process if available in $defaultdest/downloads
if {[file exists $defaultdest/downloads/geneHancerRegElements_${build}.tsv.gz]} {
	set target reg_${build}_geneHancerRegElements.tsv.zst
	file_write [gzroot $target].opt "fields\t{score elementType evidenceSources}\n"
	cg_download_ucscinfo [gzroot $target].info $build geneHancerRegElements
	cg select -overwrite 1 -hc 1 -f {chromosome=$chrom begin=$chromStart end=$chromEnd name score elementType eliteness evidenceSources} $defaultdest/downloads/geneHancerRegElements_${build}.tsv.gz $target.temp
	exec cg select -s - $target.temp | cg regcollapse | cg zst > $target.temp2.zst
	file rename -force -- $target.temp2.zst $target
	file delete $target.temp
}
if {[file exists $defaultdest/downloads/geneHancerClusteredInteractions_${build}.tsv.gz]} {
	set target extra/reg_${build}_geneHancerClusteredInteractions.tsv.zst
	cg_download_ucscinfo [gzroot $target].info $build geneHancerClusteredInteractions
	cg select -overwrite 1 -hc 1 -f {
		chromosome=$geneHancerChrom begin=$geneHancerStart end=$geneHancerEnd geneHancerIdentifier
		score interactionscore=$value geneAssociationMethods
		geneChrom geneStart geneEnd geneName geneStrand name
	} $defaultdest/downloads/geneHancerClusteredInteractions_${build}.tsv.gz $target.temp
	cg select -overwrite 1 -s - $target.temp $target.temp2
	exec cg select -s - $target.temp | cg regcollapse | cg zst > $target.temp3.zst
	file rename -force -- $target.temp3.zst $target
	file delete $target.temp $target.temp2
}

# exome target regions 
# There does not seem to be an easily accessible source for exome target regions on the net
# here we collect the ones we can (others will copied from collected downloads after this)
# 
# The following in comments were removed because they are no longer available at the url:
#	hg19 targetseq http://tools.invitrogen.com/content/sfs/manuals/TargetSeq_exome_named_targets_hg19.bed {}
#	hg19 nextera https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_nextera/nexterarapidcapture/nexterarapidcapture_exome_targetedregions.bed {}
#	hg19 nexteraexp https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_nextera/nexterarapidcapture/nexterarapidcapture_expandedexome_targetedregions.bed {}
#	hg19 nextera_v1_2 https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_nextera/nexterarapidcapture/nexterarapidcapture_exome_targetedregions_v1.2.bed {}
foreach {srcbuild targetname url file} {
	hg19 truseq_v1_2 https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/truseq/truseq-rapid-exome-targeted-regions-manifest-v1-2-bed.zip {}
	hg19 SeqCap_EZ_v2 https://sftp.rch.cm/diagnostics/sequencing/nimblegen_annotations/ez_exome_v2/SeqCapEZ_Exome_v2.0_Design_Annotation_files.zip Design_Annotation_files/Target_Regions/SeqCap_EZ_Exome_v2.bed
	hg19 SeqCap_EZ_v3 https://sftp.rch.cm/diagnostics/sequencing/literature/nimblegen/SeqCapEZ_Exome_v3.0_Design_Annotation_files.zip SeqCapEZ_Exome_v3.0_Design_Annotation_files/SeqCap_EZ_Exome_v3_hg19_capture_targets.bed
	hg19 VCRome_V2_1 https://sftp.rch.cm/diagnostics/sequencing/nimblegen_annotations/ez_hgsc_vcrome/VCRome_2.1_design_files.zip VCRome_2_1_hg19_capture_targets.bed
	hg38 twistcore https://twistbioscience.com/sites/default/files/resources/2019-06/Twist_Exome_Target_hg38.bed {}
	hg38 twistrefseq https://twistbioscience.com/sites/default/files/resources/2019-09/Twist_Exome_RefSeq_targets_hg38.bed {}
} {
	job reg_exome_$targetname -targets {
		extra/reg_${build}_exome_$targetname.tsv
	} -vars {targetname url file dest srcbuild build} -code {
		cd ${dest}/${build}
		set fulltarget [file_absolute $target]
		file delete -force $fulltarget.temp
		file mkdir $fulltarget.temp
		cd $fulltarget.temp
		if {[file extension $url] eq ".zip"} {
			wgetfile $url temp.zip
			exec unzip temp.zip
			file delete temp.zip
			if {$file eq ""} {
				set file [lindex [glob *] 0]
			}
			set f [open $file]
			set o [open temp.bed w]
			set line [gets $f]
			puts $o $line
			while {[gets $f line] != -1} {
				if {![isint [lindex $line 1]]} break
				puts $o $line
			}
			close $o
			close $f
		} else {
			wgetfile $url temp.bed
		}
		cg bed2tsv temp.bed u$targetname.tsv
		cg select -s {chromosome begin end} u$targetname.tsv s$targetname.tsv
		cg regcollapse -o reg_${srcbuild}_exome_$targetname.tsv s$targetname.tsv
		if {$srcbuild ne $build} {
			liftover_refdb reg_${srcbuild}_exome_$targetname.tsv reg_${build}_exome_$targetname.tsv $dest $srcbuild $build
			file rename -force reg_${build}_exome_$targetname.tsv.unmapped reg_${build}_exome_$targetname-unmapped.tsv
			compress reg_${build}_exome_$targetname-unmapped.tsv reg_${build}_exome_$targetname-unmapped.tsv.zst
			file rename -force -- reg_${build}_exome_$targetname-unmapped.tsv.zst [file root $fulltarget]-unmapped.tsv.zst
		}
		compress reg_${build}_exome_$targetname.tsv reg_${build}_exome_$targetname.tsv.zst
		file delete temp.bed s$targetname.tsv u$targetname.tsv
		file rename -force -- reg_${build}_exome_$targetname.tsv.zst $fulltarget.zst
		cd ${dest}/${build}
		file delete -force $fulltarget.temp
		if {$targetname eq "SeqCap_EZ_v3"} {
			mklink $fulltarget.zst [file dir $fulltarget]/reg_${build}_exome_seqcapv3.tsv.zst
		}
	}
}


# combine twist core and refseq
cd ${dest}/${build}
job reg_exome_twistfull -checkcompressed 1 -deps {
	extra/reg_${build}_exome_twistcore.tsv extra/reg_${build}_exome_twistrefseq.tsv
} -targets {
	extra/reg_${build}_exome_twistfull.tsv
} -vars {targetname url file dest build} -code {
	set fulltarget [file_absolute $target].zst
	exec cg cat {*}$dep | cg select -s - | cg regcollapse {*}[compresspipe $fulltarget] > $fulltarget.temp
	file rename -force -- $fulltarget.temp $fulltarget
}

# copy exome target regions collected in $defaultdest/downloads to extra, or lift if needed
foreach file [glob -nocomplain $defaultdest/downloads/reg_*_exome_*.zst] {
	set tail [file tail $file]
	set newtail [join [lreplace [split $tail _] 1 1 $build] _]
	if {[file exists extra/$newtail]} {
		puts "skipping extra/$newtail: already exists"
		continue
	}
	puts "Transfering $file"
	set filebuild [lindex [split $tail _] 1]
	hardcopy $file extra/
	catch {hardcopy $file.zsti extra/}
	if {$filebuild ne $build} {
		file delete extra/$newtail
		liftover_refdb extra/$tail extra/$newtail $dest $filebuild $build
	}
}

job_wait

# todo
# cp -al /complgen/refseq/hg19/gatkres /complgen/refseqnew/hg19
# cd /complgen/refseq/backup/hg19/extra
# cp reg_hg19_conserved.tsv* reg_hg19_consnoncoding.tsv* /complgen/refseqnew/hg19/extra
# rs genome_hg19* /complgen/refseqnew/hg19/extra
# cd /complgen/refseqnew/hg19
# ln -s ../hg19-local/* .

