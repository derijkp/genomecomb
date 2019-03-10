#!/bin/sh
# the next line restarts using tclsh \
exec cg source "$0" "$@"

set build hg18
set mirbasegenome hsa
set mirbaserelease 20

logverbose 2

if {![info exists argv]} {set argv {}}
set argv [job_init {*}$argv]
foreach {dest webcache} $argv break
if {![info exists dest]} {set dest /complgen/refseqnew}
if {![info exists webcache]} {set webcache $dest/webcache}
if {[info exists webcache]} {set env(webcache) $webcache}
set dest [file_absolute $dest]

putslog "Installing in $dest/$build"

# download hg18
# =============
#
file mkdir ${dest}/${build}
cd ${dest}/${build}
file mkdir extra
file mkdir ${dest}/hg19

job_logdir log_jobs

# download genome
job genome_${build} -vars build -targets {genome_${build}.ifas genome_${build}.ifas.fai extra/reg_${build}_fullgenome.tsv} -code {
	cg download_genome genome_${build}.ifas ${build}
	file rename -force reg_genome_${build}.tsv extra/reg_${build}_fullgenome.tsv
}

job genome_${build}_cindex -deps {genome_${build}.ifas} -targets {genome_${build}.ssa} -code {
	cg make_genomecindex $dep
}

job reg_${build}_sequencedgenome -vars {build} -deps {genome_${build}.ifas} -targets {extra/reg_${build}_sequencedgenome.tsv} -code {
	exec cg calcsequencedgenome --stack 1 $dep {*}[compresspipe $target 12] > $target.temp
	file rename -force $target.temp $target
}

# region databases (ucsc)
# you can explicitely download info on a database using:
# cg download_ucscinfo resultfile ${build} dbname

# collapse regions
foreach db {
	cytoBand evofold gwasCatalog microsat oreganno rmsk simpleRepeat targetScanS tfbsConsSites 
	tRNAs wgRna vistaEnhancers gad
	phastConsElements28way phastConsElements28wayPlacMammal phastConsElements44way
} {
	job reg_${build}_$db -targets {reg_${build}_${db}.tsv} -vars {build db} -code {
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
	job reg_${build}_$db -targets {reg_${build}_${db}.tsv} -vars {build db} -code {
		cg download_ucsc $target.ucsc ${build} $db
		cg regjoin $target.ucsc > $target.temp
		file delete $target.ucsc
		file rename -force $target.ucsc.info $target.info
		file rename -force $target.temp $target
	}
}

# direct
foreach db {
	firstEF rnaGene
} {
	job reg_${build}_$db -targets {reg_${build}_${db}.tsv} -vars {build db} -code {
		cg download_ucsc reg_${build}_${db}.tsv ${build} $db
	}
}

job reg_${build}_gwasCatalog -vars {build} -deps {ucsc_${build}_gwasCatalog.tsv} -targets {reg_${build}_gwasCatalog.tsv} -code {
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
job 1000g3 -targets {../hg19/var_hg19_1000g3.tsv ../hg19/extra/var_hg19_1000g3.tsv.opt} -vars {dest} -code {
	cg download_1000g3 $target hg19
	cplinked $target ../hg19/extra/var_hg19_1000g3.tsv
	file_write ../hg19/extra/var_hg19_1000g3.tsv.opt "fields\t{EUR_AF AMR_AF EAS_AF SAS_AF AFR_AF}\n"
}

job 1000g3_liftover -deps {../hg19/var_hg19_1000g3.tsv} -targets {var_${build}_1000g3.tsv} -vars {dest build} -code {
	cg liftover $dep $target ${dest}/liftover/hg19ToHg18.over.tsv
}

# dbsnp
job dbsnp130 -targets {var_${build}_snp130.tsv} -vars {dest build} -code {
	cg download_dbsnp $target ${build} snp130 2>@ stderr
}

# dbsnp147 (hg19)
job dbsnp147 -targets {../hg19/var_hg19_snp147.tsv ../hg19/var_hg19_snp147.tsv.opt} -vars {dest} -code {
	file_write $target.opt "fields\t{name}\n"
	cg download_dbsnp $target hg19 snp147 2>@ stderr
}

# dbsnp147Common (hg19)
job dbsnp147Common -targets {../hg19/var_hg19_snp147Common.tsv} -vars {dest} -code {
	cg download_dbsnp $target hg19 snp147Common 2>@ stderr
}

# liftover
job dbsnp147lift -deps {../hg19/var_hg19_snp147.tsv} -targets {var_${build}_snp147lift.tsv} -vars {dest build} -code {
	cg liftover -split 0 $dep $target ${dest}/liftover/hg19ToHg18.over.tsv
	file delete $target.unmapped
}

job dbsnp147Commonlift -deps {../hg19/var_hg19_snp147Common.tsv} -targets {var_${build}_snp147Commonlift.tsv} -vars {dest build} -code {
	cg liftover -split 0 $dep $target ../liftover/hg19ToHg18.over.tsv
	file delete $target.unmapped
}

foreach db {
	snp130 snp147lift snp147Commonlift
} {
	job maketabix_${build}_$db -deps {var_${build}_${db}.tsv} -targets {reg_${build}_${db}.tsv.gz.tbi reg_${build}_${db}.tsv.gz} -vars {dest build db} -code {
		cg maketabix $dep
	}
}

# genes
foreach db {
	refGene ensGene knownGene wgEncodeGencodeManualV3 wgEncodeGencodeAutoV3 genscan augustusAbinitio acembly
} {
	if {$db eq "wgEncodeGencodeCompV19"} {
		set dbname gencode
	} elseif {$db eq "wgEncodeGencodeManualV3"} {
		set dbname gencode
	} elseif {$db eq "wgEncodeGencodeAutoV3"} {
		set dbname gencodea
	} else {set dbname $db}
	if {$db eq "refGene"} {
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
-deps {gene_${build}_refGene.tsv extra/gene_${build}_gencode.tsv extra/gene_${build}_ensGene.tsv extra/gene_${build}_knownGene.tsv} \
-targets {$target $target.gz $target.gz.tbi} -vars {dest build db} -code {
	cg intgene {*}$deps > $target.temp
	file rename -force $target.temp $target
	cg maketabix $target
	cg index $target
}

job reg_${build}_genes -targets {extra/reg_${build}_genes.tsv} \
-deps {gene_${build}_refGene.tsv extra/gene_${build}_ensGene.tsv extra/gene_${build}_knownGene.tsv extra/gene_${build}_gencode.tsv extra/gene_${build}_gencodea.tsv} \
-code {
	exec cg cat -fields {chrom start end geneid} {*}$deps | cg select -s {chrom start end geneid} | cg regcollapse > $target.temp
	file rename -force $target.temp $target
}

job reg_refcoding \
-deps {gene_${build}_refGene.tsv} \
-targets {extra/reg_${build}_refcoding.tsv} \
-code {
	cg gene2reg $dep | cg select -q {$type eq "CDS"} | cg select -s - | cg regjoin > $target.temp
	file rename -force $target.temp $target
}

# homopolymer
job reg_${build}_homopolymer -deps {genome_${build}.ifas} -targets {reg_${build}_homopolymer.tsv reg_${build}_homopolymer.tsv.gz reg_${build}_homopolymer.tsv.gz.tbi reg_${build}_homopolymer.tsv.opt} -vars {dest build db} -code {
	cg extracthomopolymers genome_${build}.ifas > reg_${build}_homopolymer.tsv.temp
	file rename -force reg_${build}_homopolymer.tsv.temp reg_${build}_homopolymer.tsv
        cg maketabix reg_${build}_homopolymer.tsv
	file_write reg_${build}_homopolymer.tsv.opt "fields\t{base size}\n"
}

# mirbase (hg19)
job mir_${build}_mirbase -targets {../hg19/mir_hg19_mirbase$mirbaserelease.tsv ../hg19/mir_hg19_mirbase$mirbaserelease.tsv.info} -vars {mirbasegenome mirbaserelease dest build db} -code {
	cg download_mirbase $target $mirbasegenome $mirbaserelease
}

# mirbase hg18 liftover
job reg_hg18_mirbase_liftover -deps {../hg19/mir_hg19_mirbase$mirbaserelease.tsv ../hg19/mir_hg19_mirbase$mirbaserelease.tsv.info} \
-targets {mir_${build}_mirbase$mirbaserelease.tsv mir_${build}_mirbase$mirbaserelease.tsv.opt mir_${build}_mirbase$mirbaserelease.tsv.info} -vars {dest build db mirbaserelease} -code {
	file_write $target.opt "fields\t{ID}\n"
	file copy -force ../hg19/mir_hg19_mirbase$mirbaserelease.tsv.info $target.info
	cg liftover $dep $target ../liftover/hg19ToHg18.over.tsv
}

# encode
# enc_transcription
job enc_transcription -targets {reg_${build}_wgEncodelogRnaSeq.tsv} -vars {dest build} -code {
	set tempdir $target.temp
	file mkdir $tempdir
	set tables {wgEncodeCaltechRnaSeqRawSignalRep1Gm12878CellLongpolyaBb12x75 wgEncodeCaltechRnaSeqRawSignalRep1H1hescCellPapBb2R2x75 wgEncodeCaltechRnaSeqRawSignalRep2Hepg2CellPapBb2R2x75 wgEncodeCaltechRnaSeqRawSignalRep1HuvecCellPapBb2R2x75 wgEncodeCaltechRnaSeqRawSignalRep1K562CellLongpolyaBb12x75 wgEncodeCaltechRnaSeqRawSignalRep1NhekCellPapBb2R2x75}
	foreach table $tables {
		if {[file exists $tempdir/ucsc_${build}_$table.tsv]} continue
		cg download_ucsc $tempdir/ucsc_${build}_$table.tsv ${build} $table
	}
	set todo {}
	foreach table $tables {
		cg ucscwiggle2reg -n 0.01 -p 1 -f {log10($value)} $tempdir/ucsc_${build}_$table.tsv $tempdir/
		file delete $tempdir/ucsc_${build}_$table.tsv
		lappend todo $tempdir/reg_ucsc_${build}_$table.tsv
	}
	cg regcollapse -o reg_${build}_wgEncodelogRnaSeq.tsv {*}$todo
	file delete -force $tempdir
}

# make enc_transcription info file
job enc_transcription_info -targets {reg_${build}_wgEncodeCaltechRnaSeq.tsv.info} -vars {dest build} -code {
	cg download_ucscinfo $target.temp ${build} wgEncodeCaltechRnaSeq
	set c [file_read $target.temp]
	set c [string_change $c [list {== Description ==} [deindent [subst {
		== Agregation info ==
		This file combines the data from 6 cell lines.
		score contains the logarithm of the highest score, with a precision of 1 digit after the decimal point.
		num contains the number of lines for which the original value is >= 0.01 (thus log >= -2.0)
		
		== Description ==
	}]]]]
	file_write $target.temp2 $c
	file rename -force $target.temp2 $target
	file delete $target.temp
}

foreach {jobname resultname infosrc tables} {
	enc_H3K4Me1 wgEncodeH3k4me1 wgEncodeRegMarkEnhH3k4me1 {wgEncodeBroadChipSeqSignalGm12878H3k4me1 wgEncodeBroadChipSeqSignalH1hescH3k4me1 wgEncodeBroadChipSeqSignalHmecH3k4me1 wgEncodeBroadChipSeqSignalHsmmH3k4me1 wgEncodeBroadChipSeqSignalHuvecH3k4me1 wgEncodeBroadChipSeqSignalK562H3k4me1 wgEncodeBroadChipSeqSignalNhekH3k4me1 wgEncodeBroadChipSeqSignalNhlfH3k4me1}
	enc_H3K27Ac wgEncodeH3k27ac wgEncodeRegMarkEnhH3k27ac {wgEncodeBroadChipSeqSignalGm12878H3k27ac wgEncodeBroadChipSeqSignalHepg2H3k27ac wgEncodeBroadChipSeqSignalHmecH3k27ac wgEncodeBroadChipSeqSignalHsmmH3k27ac wgEncodeBroadChipSeqSignalHuvecH3k27ac wgEncodeBroadChipSeqSignalK562H3k27ac wgEncodeBroadChipSeqSignalNhekH3k27ac wgEncodeBroadChipSeqSignalNhlfH3k27ac}
	enc_H3K4Me3	wgEncodeH3k4me3 wgEncodeRegMarkPromoter   {wgEncodeBroadChipSeqSignalGm12878H3k4me3 wgEncodeBroadChipSeqSignalH1hescH3k4me3 wgEncodeBroadChipSeqSignalHepg2H3k4me3 wgEncodeBroadChipSeqSignalHmecH3k4me3 wgEncodeBroadChipSeqSignalHsmmH3k4me3 wgEncodeBroadChipSeqSignalHuvecH3k4me3 wgEncodeBroadChipSeqSignalK562H3k4me3 wgEncodeBroadChipSeqSignalNhekH3k4me3 wgEncodeBroadChipSeqSignalNhlfH3k4me3}
	
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
			cg download_ucsc $tempdir/ucsc_${build}_$table.tsv $build $table 2>@ stderr
			cg ucscwiggle2reg -n 10 -p 0 -f {5*(($value+4)/5)} $tempdir/ucsc_${build}_$table.tsv $tempdir/reg_ucsc_${build}_$table.tsv
		}
		if {![file exists $tempdir/reg_${build}_$resultname.tsv]} {
			cg regcollapse -o $tempdir/reg_${build}_$resultname.tsv.temp {*}$todo
			file rename -force $tempdir/reg_${build}_$resultname.tsv.temp $tempdir/reg_${build}_$resultname.tsv
		}
		cg bcol make --compress 9 -t iu -p begin -e end -c chromosome $tempdir/bcol_${build}_$resultname.bcol score < $tempdir/reg_${build}_$resultname.tsv
		file rename -force $tempdir/bcol_${build}_$resultname.bcol.bin.zst bcol_${build}_$resultname.bcol.bin.zst
		file rename -force $tempdir/bcol_${build}_$resultname.bcol.bin.zst.zsti bcol_${build}_$resultname.bcol.bin.zst.zsti
		file rename -force $tempdir/bcol_${build}_$resultname.bcol bcol_${build}_$resultname.bcol
		file delete -force $tempdir
	}
	# make info file
	job ${jobname}_info -targets {bcol_${build}_$resultname.tsv.info} -vars {dest build infosrc tables} -code {
		cg download_ucscinfo $target.temp ${build} $infosrc
		set c [file_read $target.temp]
		set c [string_change $c [list {== Description ==} [string trim [subst {
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
job enc_RegDnaseClustered -targets {reg_${build}_wgEncodeRegDnaseClustered.tsv reg_${build}_wgEncodeRegDnaseClustered.info} -vars {dest build} -code {
	cg download_ucsc $target.ucsc $build wgEncodeRegDnaseClusteredV3
	cg regcollapse $target.ucsc > $target.temp
	file rename -force $target.temp $target
	file delete $target.ucsc
	cg download_ucscinfo $target.info ${build} wgEncodeRegDnaseClusteredV3
}

job enc_RegTfbsClustered -targets {reg_${build}_wgEncodeRegTfbsClustered.tsv reg_${build}_wgEncodeRegTfbsClustered.info} -vars {dest build} -code {
	cg download_ucsc $target.ucsc ${build} wgEncodeRegTfbsClustered
	cg select -s - -f {chrom	start	end	name	score} $target.ucsc $target.temp
	cg regcollapse $target.temp > $target.temp2
	file rename -force $target.temp2 $target
	file rename -force $target.ucsc.info $target.info
	file delete -force $target.ucsc $target.temp
}

# link local data in dir
catch {
	cplinked {*}[glob ../hg18-local/*] .
}

# extra

# dbNSFPzip
job var_${build}_dbnsfp -targets {extra/var_${build}_dbnsfp.tsv extra/var_${build}_dbnsfp.tsv.opt} -vars {dest build} -code {
	cg download_dbnsfp $target $build ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFPv3.3a.zip 2>@ stderr
}

# compress
foreach file [jobglob *.tsv] {
	job zst_${build}_[file tail $file] -deps {$file} -targets {$file.zst} -vars {dest build} -code {
		cg zst -c 12 -i 1 $dep
	}
}

job_wait

# todo
# cd /complgen/backup2013-10/hg18/
# cp reg_hg18_conserved.tsv reg_hg18_consnoncoding.tsv reg_hg18_refgene.tsv /complgen/hg18/extra
