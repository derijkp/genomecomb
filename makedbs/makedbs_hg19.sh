#!/bin/sh
# the next line restarts using tclsh \
exec cg source "$0" "$@"

set build hg19
set mirbasegenome hsa
set mirbaserelease 20

logverbose 2

if {![info exists argv]} {set argv {}}
set argv [job_init {*}$argv]
foreach {dest webcache} $argv break
if {![info exists dest]} {set dest /complgen/refseqnew}
if {[info exists webcache]} {set env(webcache) $webcache}
set dest [file_absolute $dest]

putslog "Installing in $dest/$build"

# download hg19
# =============
#
file mkdir ${dest}/${build}
cd ${dest}/${build}
file mkdir extra
file mkdir ${dest}/hg18

job_logdir log_jobs

# download genome
job genome_${build} -vars build -targets {genome_${build}.ifas extra/reg_${build}_fullgenome.tsv} -code {
	cg download_genome genome_${build}.ifas ${build}
	file rename -force reg_genome_${build}.tsv extra/reg_${build}_fullgenome.tsv
}

job genome_${build}_cindex -deps {genome_${build}.ifas} -targets {genome_${build}.ssa} -code {
	cg make_genomecindex $dep
}

job reg_${build}_sequencedgenome -vars {dest build} -deps {genome_${build}.ifas} -targets {extra/reg_${build}_sequencedgenome.tsv} -code {
	cg calcsequencedgenome --stack 1 $dep $target
}

# region databases (ucsc)
# you can explicitely download info on a database using:
# cg download_ucscinfo resultfile ${build} dbname

# collapse regions
foreach db {
	cytoBand evofold microsat oreganno rmsk simpleRepeat targetScanS tfbsConsSites 
	tRNAs wgRna vistaEnhancers gad
	phastConsElements46way phastConsElements46wayPlacental phastConsElements46wayPrimates
} {
	job reg_${build}_$db -targets {reg_${build}_${db}.tsv} -vars {dest build db} -code {
		cg download_ucsc $target.ucsc ${build} $db
		cg regcollapse $target.ucsc > $target.temp
		file rename -force $target.ucsc.info $target.info
		file rename -force $target.temp $target
		file delete $target.ucsc
	}
}

# join regions
foreach db {
	chainSelf dgvMerged genomicSuperDups
} {
	job reg_${build}_$db -targets {reg_${build}_${db}.tsv} -vars {dest build db} -code {
		cg download_ucsc $target.ucsc ${build} $db
		cg regjoin $target.ucsc > $target.temp
		file delete $target.ucsc
		file rename -force $target.ucsc.info $target.info
		file rename -force $target.temp $target
	}
}

job reg_${build}_gwasCatalog -vars {build dest} -targets {reg_${build}_gwasCatalog.tsv} -code {
	cg download_ucsc $target.ucsc ${build} gwasCatalog
	cg select \
		-f {chrom start end trait pValue pubMedID name bin author pubDate journal title initSample replSample region genes riskAllele riskAlFreq pValueDesc orOrBeta ci95 platform cnv} \
		-nh {chrom start end name score pubMedID dbsnp bin author pubDate journal title initSample replSample region genes riskAllele riskAlFreq pValueDesc orOrBeta ci95 platform cnv} \
		$target.ucsc $target.temp
	cg regcollapse $target.temp > $target.temp2
	file rename -force $target.ucsc.info $target.info
	file rename -force $target.temp2 $target
	file delete $target.temp
	file delete $target.ucsc
}

foreach db {
	rmsk simpleRepeat
} {
	job maketabix_${build}_$db -deps {reg_${build}_${db}.tsv} -targets {reg_${build}_${db}.tsv.gz.tbi reg_${build}_${db}.tsv.gz} -vars {build db} -code {
		cg maketabix reg_${build}_${db}.tsv
	}
}

## 1000 genomes
job 1000g3 -targets {var_hg19_1000g3.tsv extra/var_hg19_1000g3.tsv.opt} -vars {dest} -code {
	cg download_1000g3 $target hg19
	cplinked $target extra/var_hg19_1000g3.tsv
	file_write extra/var_hg19_1000g3.tsv.opt "fields\t{EUR_AF AMR_AF EAS_AF SAS_AF AFR_AF}\n"
}

# dbsnp
job dbsnp147 -targets {var_hg19_snp147.tsv var_hg19_snp147.tsv.opt} -vars {dest} -code {
	file_write $target.opt "fields\t{name}\n"
	cg download_ucsc $target.temp hg19 snp147
	cg select -f {chrom start end type ref alt name freq} $target.temp $target.temp2
	file rename -force $target.temp2 $target
	file delete $target
}

job dbsnp147Common -targets {var_hg19_snp147Common.tsv} -vars {dest} -code {
	cg download_ucsc $target.temp hg19 snp147Common
	cg select -f {chrom start end type ref alt name freq} $target.temp $target.temp2
	file rename -force $target.temp2 $target
	file delete $target
}

foreach db {
	snp147 snp147Common
} {
	job maketabix_${build}_$db -deps {var_${build}_${db}.tsv} -targets {var_${build}_${db}.tsv.gz.tbi var_${build}_${db}.tsv.gz} -vars {dest build db} -code {
		cg maketabix $dep
	}
}

job clinvar -targets {var_${build}_clinvar.tsv} -vars {dest build} -code {
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
	file delete clinvar.vcf.gz
}

job kaviar -targets {var_hg19_kaviar.tsv var_hg19_kaviar.tsv.opt} -vars {dest build} -code {
	cg download_kaviar $target hg19 http://s3-us-west-2.amazonaws.com/kaviar-160204-public/Kaviar-160204-Public-hg19-trim.vcf.tar
}

# genes
foreach db {
	refGene ensGene knownGene wgEncodeGencodeBasicV24lift37 wgEncodeGencodeCompV24lift37 genscan acembly
} {
	if {$db eq "wgEncodeGencodeCompV19"} {
		set dbname gencode
	} elseif {$db eq "wgEncodeGencodeBasicV24lift37"} {
		set dbname gencode
	} elseif {$db eq "wgEncodeGencodeCompV24lift37"} {
		set dbname cgencode
	} else {set dbname $db}
	if {$db eq "refGene"} {
		set target gene_${build}_${dbname}.tsv
	} else {
		set target extra/gene_${build}_${dbname}.tsv
	}
	job gene_${build}_$dbname -targets {$target $target.gz.tbi $target.gz} -vars {dest build db} -code {
		file delete $target
		cg download_genes $target $build $db
	        cg maketabix $target
		cg index $target
	}
}

set target gene_${build}_intGene.tsv
job gene_${build}_intGene \
-deps {gene_${build}_refGene.tsv extra/gene_${build}_gencode.tsv extra/gene_${build}_ensGene.tsv extra/gene_${build}_knownGene.tsv} \
-targets {$target $target.gz $target.gz.tbi} -vars {dest build db} -code {
	cg intgene {*}$deps > $target.temp
	file rename -force $target.temp $target
	cg maketabix $target
	cg index $target
}

job reg_${build}_genes \
-deps {gene_${build}_refGene.tsv extra/gene_${build}_ensGene.tsv extra/gene_${build}_knownGene.tsv extra/gene_${build}_gencode.tsv} \
-targets {extra/reg_${build}_genes.tsv} \
-code {
	exec cg cat -fields {chrom start end geneid} {*}$deps | cg select -s {chrom start end geneid} -f {chrom {start=$start - 2000} {end=$end + 2000} geneid} | cg regcollapse > $target.temp
	file rename -force $target.temp $target
}

job reg_${build}_phenotype -deps {extra/reg_${build}_genes.tsv} \
-targets {extra/reg_${build}_phenotype.tsv extra/geneannot_${build}_phenotype.tsv} -vars {dest build} -code {
	# get target2, uses biomart (ensembl geneset) and clinvar for gene-phenotype correlations
	cg download_phenotype $target2 ${build}
	# to reg file
	cg geneannot2reg $dep $target2 $target.temp
	file rename -force $target.temp $target
}

job reg_${build}_go -deps {extra/reg_${build}_genes.tsv} \
-targets {extra/reg_${build}_go.tsv extra/geneannot_${build}_go.tsv} -code {
	cg download_mart $target2 hsapiens_gene_ensembl gene_ensembl_config {hgnc_symbol name_1006 namespace_1003}
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
	-deps {extra/reg_${build}_genes.tsv extra/geneannot_${build}_phenotype.tsv} \
	-targets {extra/reg_${build}_biograph_$name.tsv extra/geneannot_${build}_biograph_$name.tsv extra/reg_${build}_biograph_$name.info} -code {
		set genefile [tempfile]
		cg select -q "\$phenotype_description regexp \"(?i)$pattern\"" $dep2 $genefile
		cg download_biograph $target2 $id $genefile
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
	file_write reg_${build}_homopolymer.tsv.opt "fields\t{base size}\n"
}

## mirdmg
#job reg_${build}_mirdmg -targets {mir_${build}_mirdmg.tsv mir_${build}_mirdmg.info} -vars {dest build db} -code {
#	wgetfile http://mirna.bioinf.be/public/mir_hg19_mirdmg.tsv mir_hg19_mirdmg.tsv.temp
#	file rename -force mir_hg19_mirdmg.tsv.temp mir_hg19_mirdmg.tsv
#}

# mirbase
job mir_${build}_mirbase -targets {mir_${build}_mirbase$mirbaserelease.tsv mir_${build}_mirbase$mirbaserelease.tsv.info} -vars {mirbasegenome mirbaserelease dest build db} -code {
	cg download_mirbase $target $mirbasegenome $mirbaserelease
}

# exome variant server
job var_hg19_evs -targets {extra/var_hg19_evs.tsv extra/var_hg19_evs.tsv.opt extra/var_hg19_evs.tsv.info} -vars {dest build db} -code {
	cg download_evs $target hg19 http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.protein-hgvs-update.snps_indels.vcf.tar.gz
}

# exac
job reg_hg19_exac -targets {extra/var_hg19_exac.tsv extra/var_hg19_exac.tsv.opt extra/var_hg19_exac.tsv.info} -vars {dest build db} -code {
	cg download_exac $target hg19 ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz
}

# CADD
job reg_hg19_cadd -targets {extra/var_hg19_cadd.bcol extra/var_hg19_cadd.bcol.bin.lz4 extra/var_hg19_cadd.bcol.bin.lz4.lz4i extra/var_hg19_cadd.bcol.info} -vars {dest db build} -code {
	set tempdir $target.temp
	file mkdir $tempdir
	set url http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs.tsv.gz
	set tail [file tail $url]
	file_write extra/var_${build}_cadd.bcol.info [subst [deindent {
		= CADD (Combined Annotation Dependent Depletion) =
		
		== Download info ==
		dbname	cadd
		version	1.3
		citation	Kircher M, Witten DM, Jain P, O'Roak BJ, Cooper GM, Shendure J. A general framework for estimating the relative pathogenicity of human genetic variants. Nat Genet. 2014 Feb 2. doi:10.1038/ng.2892 PubMed PMID:24487276
		license	non-commercial
		source	$url
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
		wgetfile $url $tempdir/$tail
	}
	putslog "make bcol"
	file_write extra/var_${build}_cadd.tsv.opt "fields\t{score pscore}\n"
	cg select -hc 1 -rc 1 -s {Chrom Pos Alt} \
		-f {{chrom=$Chrom} {begin = $Pos - 1} {ref=$Ref} {alt=$Alt} {score=$PHRED}} $tempdir/$tail \
		| cg collapsealleles \
		| cg bcol make --precision 3 --compress 9 -t f --multicol alt --multilist A,C,T,G -p begin -c chrom $tempdir/var_${build}_cadd.bcol score
	file rename -force $tempdir/var_${build}_cadd.bcol.bin.lz4 extra/var_${build}_cadd.bcol.bin.lz4
	file rename -force $tempdir/var_${build}_cadd.bcol.bin.lz4.lz4i extra/var_${build}_cadd.bcol.bin.lz4.lz4i
	file rename -force $tempdir/var_${build}_cadd.bcol extra/var_${build}_cadd.bcol
	file delete -force $tempdir
}

# GERP
job GERP -targets {extra/bcol_${build}_GERP.bcol extra/bcol_${build}_GERP.bcol.info} -vars {dest build tables} -code {
	set table allHg19RS_BW
	set tempdir $target.temp
	file mkdir $tempdir
	set tail [file tail $target]
	cg download_ucsc $tempdir/ucsc_${build}_$table.tsv ${build} $table
	cg ucscwb2reg -p 1 -f {} $tempdir/ucsc_${build}_$table.tsv $tempdir/reg_ucsc_${build}_$table.tsv
	cg select -s - $tempdir/reg_ucsc_${build}_$table.tsv $tempdir/reg_${build}_GERP.tsv.temp
	exec cg bcol make --precision 1 --compress 9 -t f -p begin -e end -c chromosome $tempdir/$tail score < $tempdir/reg_${build}_GERP.tsv.temp
	file rename -force $tempdir/$tail.bin.lz4 $target.bin.lz4
	file rename -force $tempdir/$tail.bin.lz4.lz4i $target.bin.lz4.lz4i
	file rename -force $tempdir/$tail $target
	file rename -force $tempdir/ucsc_${build}_$table.tsv.info $target.info
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
	job $jobname -targets {bcol_${build}_$resultname.bcol} -vars {dest build tables resultname} -code {
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
job enc_RegDnaseClustered -targets {reg_${build}_wgEncodeRegDnaseClusteredV3.tsv reg_${build}_wgEncodeRegDnaseClusteredV3.tsv.info} -vars {dest build} -code {
	cg download_ucsc $target.ucsc $build wgEncodeRegDnaseClusteredV3
	cg regcollapse $target.ucsc > $target.temp
	file rename -force $target.temp $target
	file delete $target.ucsc
	cg download_ucscinfo $target.info ${build} wgEncodeRegDnaseClusteredV3
}

job enc_RegTfbsClustered -targets {reg_${build}_wgEncodeRegTfbsClusteredV3.tsv} -vars {dest build} -code {
	cg download_ucsc $target.ucsc ${build} wgEncodeRegTfbsClusteredV3
	cg select -s - -f {chrom	start	end	name	score} $target.ucsc $target.temp
	cg regcollapse $target.temp > $target.temp2
	file rename -force $target.temp2 $target
	file rename -force $target.ucsc.info $target.info
	file delete -force $target.ucsc $target.temp
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
	cg regcollapse -o reg_hg19_exome_targetseq.tsv sreg_hg19_exome_targetseq.tsv
	file delete TargetSeq_exome_named_targets_hg19.bed sreg_hg19_exome_targetseq.tsv ureg_hg19_exome_targetseq.tsv
}
# cg select -f '* {id=NR "-" $name}' reg_hg19_exome_targetseq.tsv | less

# dbNSFP
job var_hg19_dbnsfp -targets {extra/var_hg19_dbnsfp.tsv extra/var_hg19_dbnsfp.tsv.opt} -vars {dest build} -code {
	cg download_dbnsfp $target $build ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFPv3.3a.zip
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
