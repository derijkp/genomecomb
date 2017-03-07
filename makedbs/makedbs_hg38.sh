#!/bin/sh
# the next line restarts using tclsh \
exec cg source "$0" "$@"

set build hg38
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

# download
# ========
#
file mkdir ${dest}/${build}
cd ${dest}/${build}
file mkdir extra
file mkdir ${dest}/hg19

job_logdir log_jobs

# download genome
job genome_${build} -vars build -targets {genome_${build}.ifas genome_${build}.ifas.fai extra/reg_${build}_fullgenome.tsv} -code {
	cg download_genome genome_${build}.ifas ${build} 2>@ stderr
	file rename -force reg_genome_${build}.tsv extra/reg_${build}_fullgenome.tsv
}

job genome_${build}_cindex -deps {genome_${build}.ifas} -targets {genome_${build}.ssa} -code {
	cg make_genomecindex $dep
}

job reg_${build}_sequencedgenome -vars {dest build} -deps {genome_${build}.ifas} -targets {extra/reg_${build}_sequencedgenome.tsv} -code {
	exec cg calcsequencedgenome --stack 1 $dep | lz4c -12 > $target.temp
	file rename -force $target.temp $target
}

# region databases (ucsc)
# you can explicitely download info on a database using:
# cg download_ucscinfo resultfile ${build} dbname

# collapse regions
foreach db {
	cytoBand microsat oreganno rmsk simpleRepeat 
	tRNAs wgRna
	phastConsElements100way
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
job 1000g3 -targets {$dest/hg19/var_hg19_1000g3.tsv $dest/hg19/extra/var_hg19_1000g3.tsv.opt} -vars {dest} -code {
	cg download_1000g3 $dest/hg19/$target hg19
	cplinked $target $dest/hg19/extra/var_hg19_1000g3.tsv
	file_write $dest/hg19/extra/var_hg19_1000g3.tsv.opt "fields\t{EUR_AF AMR_AF EAS_AF SAS_AF AFR_AF}\n"
}

job 1000g3_liftover -deps {${dest}/hg19/var_hg19_1000g3.tsv} -targets {var_${build}_1000g3.tsv} -vars {dest build} -code {
	cg liftover ${dest}/hg19/var_hg19_1000g3.tsv var_${build}_1000g3.tsv ${dest}/liftover/hg19ToHg38.over.tsv
}

# dbsnp
job dbsnp147 -targets {var_${build}_snp147.tsv var_${build}_snp147.tsv.opt} -vars {dest build} -code {
	file_write $target.opt "fields\t{name}\n"
	cg download_dbsnp $target ${build} snp147 2>@ stderr
}

job dbsnp147Common -targets {var_${build}_snp147Common.tsv} -vars {dest build} -code {
	cg download_dbsnp $target ${build} snp147Common 2>@ stderr
}

foreach db {
	snp147 snp147Common
} {
puts pwd=[pwd]
puts "glob=[glob var_${build}_${db}.tsv]"
	job maketabix_${build}_$db -deps {var_${build}_${db}.tsv} -targets {var_${build}_${db}.tsv.gz.tbi var_${build}_${db}.tsv.gz} -vars {dest build db} -code {
		cg maketabix $dep
	}
}

job clinvar -targets {var_${build}_clinvar.tsv} -vars {dest build} -code {
	set tempdir $target.temp
	file mkdir $tempdir
	wgetfile ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz $tempdir/clinvar.vcf.gz
	cg vcf2tsv $tempdir/clinvar.vcf.gz $tempdir/clinvar.tsv
	set f [open $tempdir/clinvar.tsv]
	set header [tsv_open $f comment]
	close $f
	if {![regexp reference=GRCh38 $comment]} {
		error "clinvar.vcf.gz is from a different reference genome version"
	}
	cg collapsealleles $tempdir/clinvar.tsv > $tempdir/clinvar_collapsed.tsv
	file_write [gzroot $target].opt "fields\t{CLNACC CLNDBN}\nheaderfields\t{clinvar_acc clinvar_disease}\n"
	file rename -force $tempdir/clinvar_collapsed.tsv $target
	file delete -force $tempdir
}

job kaviar -targets {var_${build}_kaviar.tsv} -vars {dest build} -code {
	cg download_kaviar $target ${build} http://s3-us-west-2.amazonaws.com/kaviar-160204-public/Kaviar-160204-Public-${build}-trim.vcf.tar
}

# genes
foreach db {
	refGene knownGene wgEncodeGencodeBasicV24 wgEncodeGencodeCompV24
	genscan augustusGene lincRNAsTranscripts
} {
	if {$db eq "wgEncodeGencodeCompV19"} {
		set dbname gencode
	} elseif {$db eq "wgEncodeGencodeBasicV24"} {
		set dbname gencode
	} elseif {$db eq "wgEncodeGencodeCompV24"} {
		set dbname cgencode
	} elseif {$db eq "lincRNAsTranscripts"} {
		set dbname lincRNA
	} else {set dbname $db}
	if {$db in "refGene lincRNAsTranscripts"} {
		set target gene_${build}_${dbname}.tsv
	} else {
		set target extra/gene_${build}_${dbname}.tsv
	}
	job gene_${build}_$dbname -targets {$target $target.gz.tbi $target.gz} -vars {dest build db} -code {
		cg download_genes $target $build $db
	        cg maketabix $target
		cg index $target
	}
}

set target gene_${build}_intGene.tsv
job gene_${build}_intGene \
-deps {gene_${build}_refGene.tsv extra/gene_${build}_gencode.tsv extra/gene_${build}_knownGene.tsv} \
-targets {$target $target.gz $target.gz.tbi} -vars {dest build db} -code {
	cg intgene {*}$deps > $target.temp
	file rename -force $target.temp $target
	cg maketabix $target
	cg index $target
}

job reg_${build}_genes \
-deps {gene_${build}_refGene.tsv extra/gene_${build}_knownGene.tsv extra/gene_${build}_gencode.tsv} \
-targets {extra/reg_${build}_genes.tsv} \
-code {
	exec cg cat -fields {chrom start end geneid} {*}$deps | cg select -s {chrom start end geneid} -f {chrom {start=$start - 2000} {end=$end + 2000} geneid} | cg regcollapse > $target.temp
	file rename -force $target.temp $target
}

job reg_refcoding \
-deps {gene_${build}_refGene.tsv} \
-targets {extra/reg_${build}_refcoding.tsv} \
-code {
	cg gene2reg $dep | cg select -q {$type eq "CDS"} | cg select -s - | cg regjoin > $target.temp
	file rename -force $target.temp $target
}

job reg_${build}_phenotype -deps {extra/reg_${build}_genes.tsv} \
-targets {extra/reg_${build}_phenotype.tsv extra/geneannot_${build}_phenotype.tsv} -vars {dest build} -code {
	# get target2, uses biomart (ensembl geneset) and clinvar for gene-phenotype correlations
	cg download_phenotype --stack 1 $target2 ${build}
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
#job reg_${build}_mirdmg -targets {$dest/${build}/mir_${build}_mirdmg.tsv $dest/${build}/mir_${build}_mirdmg.info} -vars {dest build db} -code {
#	wgetfile http://mirna.bioinf.be/public/mir_${build}_mirdmg.tsv mir_hg19_mirdmg.tsv.temp
#	file rename -force mir_hg19_mirdmg.tsv.temp mir_hg19_mirdmg.tsv
#}

# mirbase
job mir_${build}_mirbase -targets {mir_${build}_mirbase$mirbaserelease.tsv $dest/${build}/mir_${build}_mirbase$mirbaserelease.tsv.info} -vars {mirbasegenome mirbaserelease dest build db} -code {
	cg download_mirbase $target $mirbasegenome $mirbaserelease
}

# exome variant server
job var_${build}_evs -targets {extra/var_${build}_evs.tsv extra/var_${build}_evs.tsv.opt extra/var_${build}_evs.tsv.info} -vars {dest build db} -code {
	cg download_evs $target ${build} http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.tar.gz
}

# exac
job reg_hg19_exac -targets {${dest}/hg19/extra/var_hg19_exac.tsv ${dest}/hg19/extra/var_hg19_exac.tsv.info} -vars {dest build db} -code {
	cg download_exac $target hg19 ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz
}

job exac_${build}_liftover -deps {${dest}/hg19/extra/var_hg19_exac.tsv ${dest}/hg19/extra/var_hg19_exac.tsv.info} \
	-targets {extra/var_${build}_exac.tsv extra/var_${build}_exac.tsv.info} -vars {dest build} -code {
	file copy -force $dep2 $target2
	cg liftover -split 0 $dep $target ${dest}/liftover/hg19ToHg38.over.tsv
}

# encode
foreach {jobname resultname infosrc tables} {
	enc_transcription wgEncodeRegTxn wgEncodeRegTxn {wgEncodeRegTxnCaltechRnaSeqGm12878R2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqH1hescR2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqHelas3R2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqHepg2R2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqHsmmR2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqHuvecR2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqK562R2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqNhekR2x75Il200SigPooled wgEncodeRegTxnCaltechRnaSeqNhlfR2x75Il200SigPooled}
	enc_RegMarkH3k27ac wgEncodeRegMarkH3k27ac wgEncodeRegMarkH3k27ac {wgEncodeBroadHistoneGm12878H3k27acStdSig wgEncodeBroadHistoneH1hescH3k27acStdSig wgEncodeBroadHistoneHsmmH3k27acStdSig wgEncodeBroadHistoneHuvecH3k27acStdSig wgEncodeBroadHistoneK562H3k27acStdSig wgEncodeBroadHistoneNhekH3k27acStdSig wgEncodeBroadHistoneNhlfH3k27acStdSig}
	enc_RegMarkH3k4me1 wgEncodeRegMarkH3k4me1 wgEncodeRegMarkH3k4me1 {wgEncodeBroadHistoneGm12878H3k4me1StdSig wgEncodeBroadHistoneH1hescH3k4me1StdSig wgEncodeBroadHistoneHsmmH3k4me1StdSig wgEncodeBroadHistoneHuvecH3k4me1StdSig wgEncodeBroadHistoneK562H3k4me1StdSig wgEncodeBroadHistoneNhekH3k4me1StdSig wgEncodeBroadHistoneNhlfH3k4me1StdSig}
	enc_RegMarkH3k4me3	wgEncodeRegMarkH3k4me3 wgEncodeRegMarkH3k4me3 {wgEncodeBroadHistoneGm12878H3k4me3StdSig wgEncodeBroadHistoneH1hescH3k4me3StdSig wgEncodeBroadHistoneHsmmH3k4me3StdSig wgEncodeBroadHistoneHuvecH3k4me3StdSig wgEncodeBroadHistoneK562H3k4me3StdSig wgEncodeBroadHistoneNhekH3k4me3StdSig wgEncodeBroadHistoneNhlfH3k4me3StdSig}
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
	file rename -force $target.ucsc.info $target.info
	file rename -force $target.temp $target
	file delete $target.ucsc
}

# todo : wgEncodeRegDnaseWig wgEncodeRegDnase

# link local data in dir
foreach file [glob -nocomplain ../hg19-local/*] {
	catch {
		file delete extra/[file tail $file]
		cplinked $file [file tail $file]
	}
}

if 0 {
# extra dir
# targetseq exome
job reg_exome_targetseq -targets {extra/reg_hg19_exome_targetseq.tsv} -code {
	wgetfile http://tools.invitrogen.com/content/sfs/manuals/TargetSeq_exome_named_targets_hg19.bed $target.bed
	cg bed2sft $target.bed $target.temp
	cg select -s {chromosome begin end} $target.temp $target.temp2
	cg regcollapse -o reg_hg19_exome_targetseq.tsv $target.temp2
	file delete $target.bed $target.temp $target.temp2
}
# cg select -f '* {id=NR "-" $name}' reg_hg19_exome_targetseq.tsv | less

# dbNSFP
job var_hg19_dbnsfp -targets {extra/var_hg19_dbnsfp.tsv extra/var_hg19_dbnsfp.tsv.opt} -vars {dest} -code {
	cg download_dbnsfp $target $build ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFPv3.3a.zip
}
}


job extragenome -deps {genome_${build}.ifas genome_${build}.ifas.index genome_${build}.ssa} -vars build \
-targets {extra/genome_${build}.ifas extra/genome_${build}.ifas.index extra/genome_${build}.ssa} -code {
	mklink genome_${build}.ifas extra/genome_${build}.ifas
	mklink genome_${build}.ifas.fai extra/genome_${build}.ifas.fai
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

# compress
foreach file [jobglob *.tsv] {
	job lz4_${build}_[file tail $file] -deps {$file} -targets {$file.lz4} -vars {dest build} -code {
		cg lz4 -c 12 -i 1 $dep
	}
}

# CADD
job reg_${build}_cadd -targets {extra/var_${build}_cadd.bcol extra/var_${build}_cadd.bcol.bin.lz4 extra/var_${build}_cadd.bcol.bin.lz4.lz4.lz4i extra/var_${build}_cadd.bcol.info} -vars {dest db build} -code {
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
	if {![file exists $tempdir/collapsedhg38.tsv.lz4]} {
		if {![file exists $tempdir/collapsed.tsv.lz4]} {
			cg select --stack 1 -hc 1 -rc 1 -f {{chrom=$Chrom} {begin = $Pos - 1} {end=$Pos} {ref=$Ref} {alt=$Alt} {score=$PHRED}} $tempdir/$tail \
				| cg collapsealleles --stack 1 | lz4c > $tempdir/collapsed.tsv.lz4.temp | lz4c > $tempdir/collapsed.tsv.lz4.temp
			file rename $tempdir/collapsed.tsv.lz4.temp $tempdir/collapsed.tsv.lz4
		}
		cg liftover --stack 1 -split 0 -s 0 $tempdir/collapsed.tsv.lz4 ../liftover/hg19ToHg38.over.tsv | lz4c > $tempdir/collapsedhg38.tsv.lz4.temp
		file rename $tempdir/collapsedhg38.tsv.lz4.temp $tempdir/collapsedhg38.tsv.lz4
		file delete $tempdir/collapsed.tsv.lz4
	}
	exec cg select -s - $tempdir/collapsedhg38.tsv.lz4 | cg bcol make --stack 1 --precision 3 --compress 9 -t f --multicol alt --multilist A,C,T,G -p begin -c chrom $tempdir/var_${build}_cadd.bcol score
	file rename -force $tempdir/var_${build}_cadd.bcol.bin.lz4 extra/var_${build}_cadd.bcol.bin.lz4
	file rename -force $tempdir/var_${build}_cadd.bcol.bin.lz4.lz4i extra/var_${build}_cadd.bcol.bin.lz4.lz4i
	file rename -force $tempdir/var_${build}_cadd.bcol extra/var_${build}_cadd.bcol
	file delete -force $tempdir
}

job_wait

# todo
# cd /complgen/refseq/backup/hg19/extra
# cp reg_hg19_conserved.tsv reg_hg19_consnoncoding.tsv reg_hg19_refgene.tsv /complgen/refseq/hg19/extra
# cp reg_hg19_exome_SureSelectV4.tsv reg_hg19_exome_targetseq.tsv /complgen/refseq/hg19/extra/
# rs genome_hg19* /complgen/refseq/hg19/extra

