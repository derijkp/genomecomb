#!/bin/sh
# the next line restarts using tclsh \
exec cg source "$0" "$@"

set build hg38

# settings
set mirbasegenome hsa
set mirbaserelease 21
set mirbasebuild hg38
set dbsnpversion 150
set gencodeversion 27
set 1000g3url ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz
set 1000g3readmeurl ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/README_phase3_callset_20150220
set 1000g3build hg19
set clinvarurl ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20180128.vcf.gz
set clinvarpapuurl ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20180128_papu.vcf.gz
set kaviarurl http://s3-us-west-2.amazonaws.com/kaviar-160204-public/Kaviar-160204-Public-${build}-trim.vcf.tar
set kaviarbuild hg19
set evsurl http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.protein-hgvs-update.snps_indels.vcf.tar.gz
set evsbuild hg19
set exacurl ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/ExAC.r1.sites.vep.vcf.gz
set exacbuild hg19
set caddversion 1.5
set caddurl http://krishna.gs.washington.edu/download/CADD/v$caddversion/whole_genome_SNVs.tsv.gz
set caddbuild hg19
set dbnsfpurl ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFPv3.5a.zip
set dbnsfpbuild hg38

# arguments
if {![info exists argv]} {set argv {}}
set argv [job_init {*}$argv]
foreach {dest webcache} $argv break
if {![info exists dest]} {set dest /complgen/refseqnew}
if {![info exists webcache]} {set webcache $dest/webcache}
if {[info exists webcache]} {set env(webcache) $webcache}
set dest [file_absolute $dest]

# prepare
putslog "Installing in $dest/$build"
file mkdir ${dest}/${build}
cd ${dest}/${build}
file mkdir extra

logverbose 2
job_logdir log_jobs

# download
# ========
#

# readme
set c [file_read $genomecombdir/docs/dbdir_README.txt]
regsub {version: [0-9.]+} $c "version: 0.98.7\ntime: [lindex [timestamp] 0]" c
file_write README_dbdir.txt $c

# download genome
job genome_${build} -targets {
	genome_${build}.ifas
	genome_${build}.ifas.fai
	extra/reg_${build}_fullgenome.tsv
} -vars build -code {
	cg download_genome -alt 0 genome_${build}.ifas ${build} 2>@ stderr
	file rename -force reg_genome_${build}.tsv extra/reg_${build}_fullgenome.tsv
	cg lz4 -i 1 extra/reg_${build}_fullgenome.tsv
}

job genome_${build}_cindex -deps {
	genome_${build}.ifas
} -targets {
	genome_${build}.ssa
} -code {
	cg make_genomecindex $dep
}

job reg_${build}_sequencedgenome -deps {
	genome_${build}.ifas
} -targets {
	extra/reg_${build}_sequencedgenome.tsv.lz4
} -vars {dest build} -code {
	exec cg calcsequencedgenome --stack 1 $dep | lz4c -12 > $target.temp
	file rename -force $target.temp $target
}

# make bwa version of genome
bwarefseq_job genome_${build}.ifas

# region databases (ucsc)
# you can explicitely download info on a database using:
# cg download_ucscinfo resultfile ${build} dbname

# collapse regions
foreach db {
	cytoBand microsat oreganno rmsk simpleRepeat 
	tRNAs wgRna
	phastConsElements100way phastConsElements30way
} {
	job reg_${build}_$db -targets {
		reg_${build}_${db}.tsv
	} -vars {dest build db} -code {
		set target [gzroot $target].lz4
		cg download_ucsc $target.ucsc ${build} $db
		cg regcollapse $target.ucsc > $target.temp
		file delete $target.ucsc
		cg lz4 -i 1 $target.temp
		file rename -force $target.ucsc.info [gzroot $target].info
		file rename -force $target.temp.lz4 $target
		file rename -force $target.temp.lz4.lz4i $target.lz4i
	}
}

# join regions
foreach db {
	chainSelf dgvMerged genomicSuperDups
} {
	job reg_${build}_$db -targets {
		reg_${build}_${db}.tsv
	} -vars {dest build db} -code {
		set target [gzroot $target].lz4
		cg download_ucsc $target.ucsc ${build} $db
		cg regjoin $target.ucsc > $target.temp
		file delete $target.ucsc
		cg lz4 -i 1 $target.temp
		file rename -force $target.ucsc.info [gzroot $target].info
		file rename -force $target.temp.lz4 $target
		file rename -force $target.temp.lz4.lz4i $target.lz4i
	}
}

job reg_${build}_gwasCatalog -targets {
	reg_${build}_gwasCatalog.tsv
} -vars {build dest} -code {
	set target [gzroot $target].lz4
	cg download_ucsc $target.ucsc ${build} gwasCatalog
	cg select \
		-f {chrom start end trait pValue pubMedID name bin author pubDate journal title initSample replSample region genes riskAllele riskAlFreq pValueDesc orOrBeta ci95 platform cnv} \
		-nh {chrom start end name score pubMedID dbsnp bin author pubDate journal title initSample replSample region genes riskAllele riskAlFreq pValueDesc orOrBeta ci95 platform cnv} \
		$target.ucsc $target.temp
	cg regcollapse $target.temp > $target.temp2
	file delete $target.temp
	file delete $target.ucsc
	cg lz4 -i 1 $target.temp2
	file rename -force $target.ucsc.info [gzroot $target].info
	file rename -force $target.temp2.lz4 $target
	file rename -force $target.temp2.lz4.lz4i $target.lz4i
}

foreach db {
	rmsk simpleRepeat
} {
	job maketabix_${build}_$db -deps {
		reg_${build}_${db}.tsv
	} -targets {
		reg_${build}_${db}.tsv.gz.tbi
		reg_${build}_${db}.tsv.gz
	} -vars {build db} -code {
		cg maketabix $dep
	}
}

## 1000 genomes
job 1000g3 -targets {
	var_${build}_1000g3.tsv
	extra/var_${build}_1000g3.tsv
	extra/var_${build}_1000g3.tsv.opt
} -vars {dest 1000g3url 1000g3readmeurl 1000g3build build} -code {
	set target [gzroot $target].lz4
	if {$1000g3build eq $build} {
		cg download_1000g3 $target $1000g3url $1000g3readmeurl
	} else {
		cg download_1000g3 $target.$1000g3build.lz4 $1000g3url $1000g3readmeurl
		liftover_refdb $target.$1000g3build.lz4 $target $dest $1000g3build $build
	}
	cplinked $target extra/var_${build}_1000g3.tsv.lz4
	cplinked $target.lz4i extra/var_${build}_1000g3.tsv.lz4.lz4i
	file_write extra/var_${build}_1000g3.tsv.opt "fields\t{EUR_AF AMR_AF EAS_AF SAS_AF AFR_AF}\n"
}

# dbsnp
job dbsnp -targets {
	var_${build}_dnsnp.tsv
	var_${build}_dbsnp.tsv.opt
} -vars {dest build dbsnpversion} -code {
	set target [gzroot $target].lz4
	file_write [gzroot $target].opt "fields\t{name}\n"
	cg download_dbsnp $target ${build} snp$dbsnpversion 2>@ stderr
	cg lz4index $target
}

job dbsnpCommon -targets {
	var_${build}_dbsnpCommon.tsv
} -vars {dest build dbsnpversion} -code {
	set target [gzroot $target].lz4
	file_write [gzroot $target].opt "fields\t{freqp}\n"
	cg download_dbsnp $target ${build} snp${dbsnpversion}Common 2>@ stderr
	cg lz4index $target
}

foreach db [list dbsnp dbsnpCommon] {
	job maketabix_${build}_$db -deps {
		var_${build}_${db}.tsv
	} -targets {
		var_${build}_${db}.tsv.gz.tbi
		var_${build}_${db}.tsv.gz
	} -vars {dest build db} -code {
		cg maketabix $dep
	}
}

job clinvar -targets {
	var_${build}_clinvar.tsv
} -vars {dest build clinvarurl clinvarpapuurl} -code {
	set target [gzroot $target].lz4
	cg download_clinvar --stack 1 $target $build $clinvarurl $clinvarpapuurl
	cg lz4index $target
}

job kaviar -targets {
	var_${build}_kaviar.tsv
} -vars {dest build kaviarurl kaviarbuild} -code {
	set target [gzroot $target].lz4
	if {$kaviarbuild eq $build} {
		cg download_kaviar $target $kaviarurl
	} else {
		cg download_kaviar $target.$kaviarbuild.lz4 $kaviarurl
		liftover_refdb $target.$kaviarbuild.lz4 $target $dest $kaviarbuild $build
	}
}

# genes
foreach db [list \
	refGene knownGene wgEncodeGencodeBasicV${gencodeversion} wgEncodeGencodeCompV${gencodeversion} \
	genscan augustusGene lincRNAsTranscripts \
] {
	if {$db eq "wgEncodeGencodeCompV19"} {
		set dbname gencode
	} elseif {$db eq "wgEncodeGencodeBasicV${gencodeversion}"} {
		set dbname gencode
	} elseif {$db eq "wgEncodeGencodeCompV${gencodeversion}"} {
		set dbname cgencode
	} elseif {$db eq "lincRNAsTranscripts"} {
		set dbname lincRNA
	} else {set dbname $db}
	if {$db in "refGene lincRNAsTranscripts"} {
		set target gene_${build}_${dbname}.tsv
	} else {
		set target extra/gene_${build}_${dbname}.tsv
	}
	job gene_${build}_$dbname -targets {
		$target
		$target.gz.tbi
		$target.gz
	} -vars {dest build db} -code {
		set target [gzroot $target].lz4
		file delete $target
		cg download_genes $target $build $db
	        cg maketabix $target
		cg index $target
	}
}

set target gene_${build}_intGene.tsv
job gene_${build}_intGene -deps {
	gene_${build}_refGene.tsv
	extra/gene_${build}_gencode.tsv
	extra/gene_${build}_knownGene.tsv
} -targets {
	$target
	$target.gz
	$target.gz.tbi
} -vars {dest build db} -code {
	set target [gzroot $target].lz4
	cg intgene {*}$deps | cg lz4 > $target.temp.lz4
	file rename -force $target.temp.lz4 $target
	cg maketabix $target
	cg lz4index $target
	cg index $target
}

job reg_${build}_genes -deps {
	gene_${build}_refGene.tsv
	extra/gene_${build}_knownGene.tsv
	extra/gene_${build}_gencode.tsv
} -targets {
	extra/reg_${build}_genes.tsv
} -code {
	set target [gzroot $target].lz4
	exec cg cat -fields {chrom start end geneid} {*}$deps | cg select -s {chrom start end geneid} -f {chrom {start=$start - 2000} {end=$end + 2000} geneid} | cg regcollapse | cg lz4 > $target.temp.lz4
	file rename -force $target.temp.lz4 $target
	cg lz4index $target
}

job reg_refcoding -deps {
	gene_${build}_refGene.tsv
} -targets {
	extra/reg_${build}_refcoding.tsv
} -code {
	set target [gzroot $target].lz4
	cg gene2reg $dep | cg select -q {$type eq "CDS"} | cg select -s - | cg regjoin | cg lz4 > $target.temp.lz4
	file rename -force $target.temp.lz4 $target
	cg lz4index $target
}

job reg_exome_refGene -deps {
	gene_${build}_refGene.tsv
} -targets {
	extra/reg_${build}_exome_refGene.tsv
} -code {
	set target [gzroot $target].lz4
	cg gene2reg $dep | cg select -q {$type in "CDS UTR RNA"} | cg select -s - | cg regjoin | cg lz4 > $target.temp.lz4
	file rename -force $target.temp.lz4 $target
	cg lz4index $target
}

job reg_intcoding -deps {
	gene_${build}_intGene.tsv
} -targets {
	extra/reg_${build}_intcoding.tsv
} -code {
	set target [gzroot $target].lz4
	cg gene2reg $dep | cg select -q {$type eq "CDS"} | cg select -s - | cg regjoin | cg lz4 > $target.temp.lz4
	file rename -force $target.temp.lz4 $target
	cg lz4index $target
}

job reg_exome_intGene -deps {
	gene_${build}_intGene.tsv
} -targets {
	extra/reg_${build}_exome_intGene.tsv
} -code {
	set target [gzroot $target].lz4
	cg gene2reg $dep | cg select -q {$type in "CDS UTR RNA"} | cg select -s - | cg regjoin | cg lz4 > $target.temp.lz4
	file rename -force $target.temp.lz4 $target
	cg lz4index $target
}

file_write extra/reg_${build}_pseudoautosomal.tsv {chromosome	begin	end	name
X	10001	2781479	PAR1
X	155701383	156030895	PAR2
Y	10001	2781479	PAR1
Y	56887903	57217415	PAR2
}

job reg_${build}_phenotype -deps {
	extra/reg_${build}_genes.tsv
} -targets {
	extra/reg_${build}_phenotype.tsv
	extra/geneannot_${build}_phenotype.tsv
} -vars {dest build} -code {
	set target [gzroot $target].lz4
	# get target2, uses biomart (ensembl geneset) and clinvar for gene-phenotype correlations
	cg download_phenotype $target2 ${build}
	# to reg file
	cg geneannot2reg $dep $target2 $target.temp
	cg lz4 -i 1 $target.temp
	file rename -force $target.temp.lz4 $target
	file rename -force $target.temp.lz4.lz4i $target.lz4i
}

job reg_${build}_go -deps {
	extra/reg_${build}_genes.tsv
} -targets {
	extra/reg_${build}_go.tsv
	extra/geneannot_${build}_go.tsv
} -code {
	set target [gzroot $target].lz4
	cg download_mart $target2.lz4 hsapiens_gene_ensembl gene_ensembl_config {hgnc_symbol name_1006 namespace_1003}
	cg geneannot2reg $dep $target2.lz4 $target.temp
	cg lz4 -i 1 $target.temp
	file rename -force $target.temp.lz4 $target
	file rename -force $target.temp.lz4.lz4i $target.lz4i
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
#		cg download_biograph $target2.lz4 $id $genefile
#		cg geneannot2reg $dep $target2.lz4 $target.temp
#		cg lz4 -i 1 $target.temp
#		file rename -force $target.temp.lz4 [gzroot $target].lz4
#		file rename -force $target.temp.lz4.lz4i [gzroot $target].lz4.lz4i
#		set temp "gene ranking\n============\ngene ranking for $disease\n\nThe following scores will be given to variants in genes that"
#		append temp "0: are directly linked to the disease in the ensembl gene database or in the clinvar database\n"
#		append temp "1: have a known relation to the disease in biograph (can be association, etc.)\n"
#		append temp ">1: position in genelist sorted according to biograph prioritization (there are 18182 genes)"
#		file_write $target3 $temp
#	}
#}

# homopolymer
job reg_${build}_homopolymer -deps {
	genome_${build}.ifas
} -targets {
	reg_${build}_homopolymer.tsv.lz4
	reg_${build}_homopolymer.tsv.gz
	reg_${build}_homopolymer.tsv.gz.tbi
	reg_${build}_homopolymer.tsv.opt
} -vars {dest build db} -code {
	file_write reg_${build}_homopolymer.tsv.opt "fields\t{base size}\n"
	cg extracthomopolymers genome_${build}.ifas | cg lz4 > reg_${build}_homopolymer.tsv.temp.lz4
	file rename -force reg_${build}_homopolymer.tsv.temp.lz4 reg_${build}_homopolymer.tsv.lz4
        cg maketabix reg_${build}_homopolymer.tsv.lz4
}

# mirbase
job mir_${build}_mirbase -targets {
	mir_${build}_mirbase$mirbaserelease.tsv
	mir_${build}_mirbase$mirbaserelease.tsv.info
} -vars {mirbasegenome mirbaserelease mirbasebuild dest build db} -code {
	set target [gzroot $target].lz4
	if {$mirbasebuild ne $build} {error "error: mirbase $mirbaserelease for build $mirbasebuild (making $build)"}
	cg download_mirbase $target $mirbasegenome $mirbaserelease
}

# exome variant server
job var_${build}_evs -targets {
	extra/var_${build}_evs.tsv
	extra/var_${build}_evs.tsv.opt
	extra/var_${build}_evs.tsv.info
} -vars {dest build db evsurl evsbuild} -code {
	set target [gzroot $target].lz4
	if {$evsbuild eq $build} {
		cg download_evs $target $evsurl
	} else {
		cg download_evs $target.$evsbuild.lz4 $evsurl
		liftover_refdb $target.$evsbuild.lz4 $target $dest $evsbuild $build
	}
}

# exac
job reg_${build}_exac -targets {
	extra/var_${build}_exac.tsv
	extra/var_${build}_exac.tsv.opt
	extra/var_${build}_exac.tsv.info
} -vars {dest build db exacurl exacbuild} -code {
	set target [gzroot $target].lz4
	if {$exacbuild eq $build} {
		cg download_exac --stack 1 --verbose 2 $target $exacurl
	} else {
		cg download_exac --stack 1 --verbose 2 $target.$exacbuild.lz4 $exacurl
		liftover_refdb $target.$exacbuild.lz4 $target $dest $exacbuild $build
	}
}

# encode
foreach {jobname resultname infosrc tables} {
	enc_transcription wgEncodeRegTxn wgEncodeRegTxn {wgEncodeRegTxnCaltechRnaSeqGm12878R2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqH1hescR2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqHelas3R2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqHepg2R2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqHsmmR2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqHuvecR2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqK562R2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqNhekR2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqNhlfR2x75Il200SigPooled}
	enc_RegMarkH3k27ac wgEncodeRegMarkH3k27ac wgEncodeRegMarkH3k27ac {wgEncodeBroadHistoneGm12878H3k27acStdSig wgEncodeBroadHistoneH1hescH3k27acStdSig wgEncodeBroadHistoneHsmmH3k27acStdSig wgEncodeBroadHistoneHuvecH3k27acStdSig wgEncodeBroadHistoneK562H3k27acStdSig wgEncodeBroadHistoneNhekH3k27acStdSig wgEncodeBroadHistoneNhlfH3k27acStdSig}
	enc_RegMarkH3k4me1 wgEncodeRegMarkH3k4me1 wgEncodeRegMarkH3k4me1 {wgEncodeBroadHistoneGm12878H3k4me1StdSig wgEncodeBroadHistoneH1hescH3k4me1StdSig wgEncodeBroadHistoneHsmmH3k4me1StdSig wgEncodeBroadHistoneHuvecH3k4me1StdSig wgEncodeBroadHistoneK562H3k4me1StdSig wgEncodeBroadHistoneNhekH3k4me1StdSig wgEncodeBroadHistoneNhlfH3k4me1StdSig}
	enc_RegMarkH3k4me3	wgEncodeRegMarkH3k4me3 wgEncodeRegMarkH3k4me3 {wgEncodeBroadHistoneGm12878H3k4me3StdSig wgEncodeBroadHistoneH1hescH3k4me3StdSig wgEncodeBroadHistoneHsmmH3k4me3StdSig wgEncodeBroadHistoneHuvecH3k4me3StdSig wgEncodeBroadHistoneK562H3k4me3StdSig wgEncodeBroadHistoneNhekH3k4me3StdSig wgEncodeBroadHistoneNhlfH3k4me3StdSig}
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
			file rename -force $tempdir/reg_${build}_$resultname.tsv.temp $tempdir/reg_${build}_$resultname.tsv
		}
		cg bcol make --compress 9 -t iu -p begin -e end -c chromosome $tempdir/bcol_${build}_$resultname.bcol score < $tempdir/reg_${build}_$resultname.tsv
		file rename -force $tempdir/bcol_${build}_$resultname.bcol.bin.lz4 bcol_${build}_$resultname.bcol.bin.lz4
		file rename -force $tempdir/bcol_${build}_$resultname.bcol.bin.lz4.lz4i bcol_${build}_$resultname.bcol.bin.lz4.lz4i
		file rename -force $tempdir/bcol_${build}_$resultname.bcol bcol_${build}_$resultname.bcol
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
		file rename -force $target.temp2 $target
		file delete $target.temp
	}
}

# DNase Clusters and Txn Factor ChIP
job enc_RegDnaseClustered -targets {
	reg_${build}_wgEncodeRegDnaseClustered.tsv
	reg_${build}_wgEncodeRegDnaseClustered.tsv.info
} -vars {dest build} -code {
	cg download_ucsc $target.ucsc $build wgEncodeRegDnaseClustered
	cg regcollapse $target.ucsc > $target.temp
	file rename -force $target.ucsc.info [gzroot $target].info
	file rename -force $target.temp $target
	file delete $target.ucsc
}

# todo : wgEncodeRegDnaseWig wgEncodeRegDnase

# link local data in dir
foreach file [glob -nocomplain ../${build}-local/*] {
	catch {
		file delete extra/[file tail $file]
		cplinked $file [file tail $file]
	}
}

# dbNSFP
job var_${build}_dbnsfp -targets {
	extra/var_${build}_dbnsfp.tsv
	extra/var_${build}_dbnsfp.tsv.opt
} -vars {dest build dbnsfpurl dbnsfpbuild} -code {
	set target [gzroot $target].lz4
	cg download_dbnsfp $target $build $dbnsfpurl $dbnsfpbuild
}


job extragenome -deps {
	genome_${build}.ifas
	genome_${build}.ifas.index
	genome_${build}.ssa
} -vars build -targets {
	extra/genome_${build}.ifas extra/genome_${build}.ifas.fai extra/genome_${build}.ifas.index
	genome_${build}.fa genome_${build}.fa.fai genome_${build}.fa.index
	extra/genome_${build}.ssa
} -code {
	mklink genome_${build}.ifas extra/genome_${build}.ifas
	mklink genome_${build}.ifas.fai extra/genome_${build}.ifas.fai
	mklink genome_${build}.ifas.index extra/genome_${build}.ifas.index 
	mklink genome_${build}.ifas genome_${build}.fa
	mklink genome_${build}.ifas.fai genome_${build}.fa.fai
	mklink genome_${build}.ifas.index genome_${build}.fa.index 
	mklink genome_${build}.ssa extra/genome_${build}.ssa
}
# genome in extra
foreach file [glob genome_*] {
	catch {
		file delete extra/[file tail $file]
		mklink $file extra/[file tail $file]
	}
}

# CADD
job reg_${build}_cadd -targets {
	var_${build}_cadd.bcol
	var_${build}_cadd.bcol.bin.lz4
	var_${build}_cadd.bcol.bin.lz4.lz4i
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
		if {![file exists $tempdir/collapsed${build}.tsv.lz4]} {
			if {![file exists $tempdir/collapsed.tsv.lz4]} {
				cg select --stack 1 -hc 1 -rc 1 -f {{chrom=$Chrom} {begin = $Pos - 1} {end=$Pos} {ref=$Ref} {alt=$Alt} {score=$PHRED}} $tempdir/$tail \
					| cg collapsealleles --stack 1 | lz4c > $tempdir/collapsed.tsv.lz4.temp
				file rename $tempdir/collapsed.tsv.lz4.temp $tempdir/collapsed.tsv.lz4
			}
			cg liftover --stack 1 -split 0 -s 0 $tempdir/collapsed.tsv.lz4 ../liftover/${caddbuild}To${build}.over.tsv | lz4c > $tempdir/collapsed${build}.tsv.lz4.temp
			file rename $tempdir/collapsed${build}.tsv.lz4.temp $tempdir/collapsed${build}.tsv.lz4
			file delete $tempdir/collapsed.tsv.lz4
		}
		exec cg select -s - $tempdir/collapsed${build}.tsv.lz4 | cg bcol make --stack 1 --precision 3 --compress 9 -t f --multicol alt --multilist A,C,T,G -p begin -c chrom $tempdir/var_${build}_cadd.bcol score
	} else {
		cg select -hc 1 -rc 1 -f {{chrom=$Chrom} {begin = $Pos - 1} {ref=$Ref} {alt=$Alt} {score=$PHRED}} $tempdir/$tail \
			| cg collapsealleles | cg select -s {chrom begin} \
			| cg bcol make --precision 3 --compress 9 -t f --multicol alt --multilist A,C,T,G -p begin -c chrom $tempdir/var_${build}_cadd.bcol score
	}
	file rename -force $tempdir/var_${build}_cadd.bcol.bin.lz4 var_${build}_cadd.bcol.bin.lz4
	file rename -force $tempdir/var_${build}_cadd.bcol.bin.lz4.lz4i var_${build}_cadd.bcol.bin.lz4.lz4i
	file rename -force $tempdir/var_${build}_cadd.bcol var_${build}_cadd.bcol
	file delete -force $tempdir
}

job_wait
