#!/bin/sh
# the next line restarts using tclsh \
exec cg source "$0" "$@"

# Script to download and prepare all reference databases (including genome) for hg38
# call using
# makedbs_hg38.sh ?options? ?dest? ?webcache?
# where 
# * dest: the directory where the reference databases will be installed (default /complgen/refseqnew/hg38)
# * webcache: directory that will be used to cache downloads (default /complgen/refseqnew/webcache)
# options:
# -d distr: allow distributed processing on 
#           e.g. a grid engine cluster with -d sge, 
#           or locally on multiple cores with e.g. -d 8 for 8 cores
# If the command is interrupted, you can just restart using the same
# command, and it will continue from where it was stopped.

# settings
# ========

# basic settings (for makerefdb_job)
# --------------
set build hg38
set defaultdest /complgen/refseqnew

set genomeurl http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
set par {chromosome	begin	end	name
X	10001	2781479	PAR1
X	155701383	156030895	PAR2
Y	10001	2781479	PAR1
Y	56887903	57217415	PAR2
}
set organelles {chromosome
chrM
}
set dbsnpversion 155
# set refSeqFuncElemsurl https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/109.20200522/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz
set refSeqFuncElemsurl https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/110/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz
set mirbase hsa-22.1:hg38
set regionsdb_collapse {
	cytoBand microsat oreganno rmsk simpleRepeat 
	tRNAs wgRna
	phastConsElements100way phastConsElements30way
}
set regionsdb_join {
	chainSelf dgvMerged genomicSuperDups
}

set gencodeversion 45
# list with geneset name (first word) and one or more of the following keywords
# int : include in intGene
# extra : place in the extra dir instead of in base annotation dir
# reg : make a region file from it (in extra)
set genesdb [list \
	{refGene int reg} \
	{ncbiRefSeq extra int reg} \
	[list wgEncodeGencodeBasicV${gencodeversion} gencode extra int reg] \
	{knownGene extra int reg} \
	[list wgEncodeGencodeCompV${gencodeversion} cgencode extra] \
	{genscan extra} \
	{augustusGene extra} \
	{lincRNAsTranscripts lincRNA} \
]

# extra settings
# --------------
set 1000g3url http://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz
set 1000g3readmeurl http://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/README_phase3_callset_20150220
set 1000g3build hg19
set clinvarurl https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20240407.vcf.gz
set clinvarpapuurl https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20240407_papu.vcf.gz
#set kaviarurl http://s3-us-west-2.amazonaws.com/kaviar-160204-public/Kaviar-160204-Public-${build}-trim.vcf.tar
#set kaviarbuild hg19
set evsurl http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.protein-hgvs-update.snps_indels.vcf.tar.gz
set evsbuild hg19
set exacurl ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/ExAC.r1.sites.vep.vcf.gz
set exacbuild hg19
set caddversion 1.7
set caddurl http://krishna.gs.washington.edu/download/CADD/v$caddversion/GRCh38/whole_genome_SNVs.tsv.gz
set caddbuild hg38
set gnomadbuild hg38
set gnomadversion 4.1
set gnomadbaseurl https://storage.googleapis.com/gcp-public-data--gnomad/release/$gnomadversion
set gnomadexbuild hg38
set gnomadexversion 4.1
set gnomadexurl https://storage.googleapis.com/gcp-public-data--gnomad/release/$gnomadexversion/vcf/exomes
set gnomadlofbuild hg38
set gnomadlofversion 4.1
set gnomadlof https://storage.googleapis.com/gcp-public-data--gnomad/release/$gnomadlofversion/constraint/gnomad.v$gnomadlofversion.constraint_metrics.tsv
set gnomadsvbuild hg38
set gnomadsvversion 4.1
set gnomadsv https://storage.googleapis.com/gcp-public-data--gnomad/release/$gnomadsvversion/genome_sv/gnomad.v$gnomadsvversion.sv.sites.vcf.gz
# from https://github.com/quinlan-lab/ccrhtml
set ccrversion 2.20180420
set ccrurl https://s3.us-east-2.amazonaws.com/ccrs/ccrs/ccrs.autosomes.v${ccrversion}.bed.gz
set ccrbuild hg19
# set dbnsfpurl ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFP4.3a.zip
set dbnsfpurl https://dbnsfp.s3.amazonaws.com/dbNSFP4.5a.zip
set dbnsfpbuild hg38
#set gtfurl https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ensGene.gtf.gz
#set gtffile genes_hg38_ensGene.gtf.gz
#set gencodegtfurl ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz
#set gencodegtffile extra/gene_hg38_gencode.v39.gtf
set transcriptsurl http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz
set transcriptsgtf extra/gene_hg38_gencode.v45.gtf
set dm_alphamissense_url https://storage.googleapis.com/dm_alphamissense/AlphaMissense_isoforms_hg38.tsv.gz
set dm_alphamissense_canonical_url https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz

# prepare
# =======

if {![info exists argv]} {set argv {}}

# keep actual command line used for log
set cmdline [clean_cmdline cg [info script] {*}$args]

# arguments, start job system
logverbose 2
set argv [job_init {*}$argv]
cg_options makedbs_hg38.sh argv {
} {dest webcache} 0 2

# foreach {dest webcache} $argv break

if {![info exists dest] || $dest eq ""} {set dest $defaultdest}
if {![info exists webcache] || $webcache eq ""} {set webcache $dest/webcache}
if {[info exists webcache]} {set env(webcache) $webcache}
set dest [file_absolute $dest]

# prepare
putslog "Installing in $dest/$build"
file mkdir ${dest}/${build}
cd ${dest}/${build}
file mkdir extra

# set 
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
	-pseudoautosomal $par \
	-organelles $organelles \
	-regionsdb_collapse $regionsdb_collapse \
	-regionsdb_join $regionsdb_join \
	-dbsnp $dbsnpversion \
	-refSeqFuncElemsurl $refSeqFuncElemsurl \
	-genesdb $genesdb \
	-mirbase $mirbase \
	-transcriptsurl $transcriptsurl \
	-transcriptsgtf $transcriptsgtf \
	-webcache $webcache \
	$dest/$build

# rest after this is hg38 specific code
# -------------------------------------

set target extra/reg_${ref}_distrg.tsv
job reg_${ref}_nolowgene -deps {
	extra/gene_${ref}_cgencode.tsv.zst
	extra/reg_${ref}_nolowgene200k.tsv.zst
} -targets {
	$target
} -vars {ref} -code {
	proc sizes {list} {
		set list [split [string trim $list] \n]
		set list [list_subindex $list 0]
		foreach region $list {
			if {$region eq ""} continue
			puts "$region\t[format %.2f [expr {([lindex [split $region -] end] - [lindex [split $region -] end-1])/1000000.0}]]"
		}
	}
	
	proc finerregions {nogenes problemregions cutoff} {
		foreach problemregion $problemregions {
			foreach {chr prevbegin finalend} [split $problemregion -] break
			set data [split [cg select -sh /dev/null \
				-q "region(\"$problemregion\") and \$end - \$begin > $cutoff" \
				$nogenes] \n]
			set regions {}
			list_foreach {c begin end} $data {
				set mid [expr {($end+$begin)/2}]
				lappend regions $chr-$prevbegin-$mid
				set prevbegin $mid
			}
			lappend regions $chr-$prevbegin-$finalend
			foreach region $regions {
				puts "$region\t[format %.2f [expr {([lindex [split $region -] end] - [lindex [split $region -] end-1])/1000000.0}]]"
			}
		}
		return $regions
	}
	set problemregions {
		chr11-54985028-79789552 chr19-0-25900983
		chr21-5371491-11824755
		chr17-24847945-52275265 chr16-0-9238611
		chr19-25900983-56983846
		chr1-144101052-157398672
		chr3-43888438-61406513
		chr9-124145453-138394717
		chr12-47540720-58325029
		chr15-22326183-37706401
		chr16-66095566-75984075
	}
	# sizes [join [bsort $problemregions] \n]
	set regions [distrreg_regs g5000000 [pwd] s 0]
	# These regions are to large still from the 200k file, remove and replace with smaller
	set regions [list_lremove $regions $problemregions]
	lappend regions {*}{
		chr1-144101052-147845873
		chr1-147845873-152478129
		chr1-152478129-157398672
		chr3-43888438-48126125
		chr3-48126125-55120733
		chr3-55120733-61406513
		chr9-124145453-126754848
		chr9-126754848-132120962
		chr9-132120962-135408505
		chr9-135408505-138394717
		chr11-54985028-57948734
			chr11-57948734-62009676
			chr11-62009676-67845906
			chr11-67845906-69232857
			chr11-69232857-71317516
		chr11-71317516-72482417
		chr11-72482417-76282961
		chr11-76282961-79789552
		chr12-47540720-51553014
		chr12-51553014-54758975
		chr12-54758975-58325029
		chr15-22326183-26330477
		chr15-26330477-31278860
		chr15-31278860-37706401
		chr16-0-4148279
		chr16-4148279-9238611
		chr16-66095566-69515678
		chr16-69515678-75984075
		chr17-24847945-27710211
		chr17-27710211-34421252
		chr17-24847945-34421252
		chr17-34421252-36794009
		chr17-36794009-46430772
		chr17-46430772-52275265
		chr19-0-5398114
		chr19-5398114-13682392
		chr19-13682392-20866608
		chr19-20866608-25900983
		chr19-25900983-30152779
		chr19-30152779-34119393
		chr19-34119393-40142915
		chr19-40142915-47927180
		chr19-47927180-56983846
		chr21-5371491-7569202
		chr21-7569202-9450211
		chr21-9450211-11824755
	}
	# set nogenes /complgen/refseqnew/hg38/extra/reg_hg38_nogenes.tsv.zst
	# sizes {}
	# finerregions $nogenes $problemregions 100000
	set o [open $target.temp w]
	puts $o [join {chromosome begin end} \t]
	foreach region [bsort $regions] {
		append region --
		puts $o [join [lrange [split $region -] 0 2] \t]
	}
	close $o
	file rename -force $target.temp $target
	# check
	cg regsubtract extra/reg_${ref}_fullgenome.tsv.zst $target | cg covered
}

#job gtf -targets {
#	$gtffile
#} -vars {gtfurl gtffile} -code {
#	wgetfile $gtfurl $gtffile.temp
#	file rename $gtffile.temp $gtffile
#}

job sniffles_trf -optional 1 -targets {
	extra/sniffles_hg38.trf.bed
} -vars {} -code {
	wgetfile https://github.com/fritzsedlazeck/Sniffles/raw/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed extra/sniffles_hg38.trf.bed.temp
	exec cg bed2tsv extra/sniffles_hg38.trf.bed.temp \
		| cg select -f {{chromosome=if($chromosome regexp "^chr",$chromosome,"chr$chromosome")} begin end} \
		| cg select -s chromosome \
		| cg tsv2bed > extra/sniffles_hg38.trf.bed.temp2
	file rename -force extra/sniffles_hg38.trf.bed.temp2 extra/sniffles_hg38.trf.bed
	file delete extra/sniffles_hg38.trf.bed.temp
}

set target extra/collapsed[file tail $transcriptsgtf]
job collapsedgencodegtf -deps {
	$transcriptsgtf
} -targets {
	$target
} -vars {transcriptsgtf} -code {
	set tempdir [tempdir]
	wgetfile https://raw.githubusercontent.com/broadinstitute/gtex-pipeline/TOPMed_RNAseq_v2/gene_model/collapse_annotation.py $tempdir/collapse_annotation.py
	exec [findpython3] $tempdir/collapse_annotation.py $transcriptsgtf $target.temp
	file rename -force $target.temp $target
}

job reg_${build}_gwasCatalog -targets {
	reg_${build}_gwasCatalog.tsv.zst
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

## 1000 genomes
job 1000g3 -targets {
	var_${build}_1000g3.tsv.zst
	extra/var_${build}_1000g3.tsv.zst
	extra/var_${build}_1000g3.tsv.opt
} -vars {dest 1000g3url 1000g3readmeurl 1000g3build build} -code {
	set target [gzroot $target].zst
	if {$1000g3build eq $build} {
		cg download_1000g3 $target $1000g3url $1000g3readmeurl
	} else {
		cg download_1000g3 $target.$1000g3build.zst $1000g3url $1000g3readmeurl
		liftover_refdb $target.$1000g3build.zst $target $dest $1000g3build $build
		file delete $target.$1000g3build.zst
		mv $target.$1000g3build.info [gzroot $target].info
	}
	cg zindex $target
	cplinked $target extra/var_${build}_1000g3.tsv.zst
	cplinked $target.zsti extra/var_${build}_1000g3.tsv.zst.zsti
	file_write extra/var_${build}_1000g3.tsv.opt "fields\t{EUR_AF AMR_AF EAS_AF SAS_AF AFR_AF}\n"
}

job clinvar -targets {
	var_${build}_clinvar.tsv.zst
} -vars {dest build clinvarurl clinvarpapuurl} -code {
	set target [gzroot $target].zst
	cg download_clinvar --stack 1 $target $build $clinvarurl $clinvarpapuurl
	cg zindex $target
}

if {[info exists kaviarurl]} {
	job kaviar -targets {
		var_${build}_kaviar.tsv.zst
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
	extra/reg_${build}_genes.tsv.zst
} -targets {
	extra/reg_${build}_phenotype.tsv.zst
	extra/geneannot_${build}_phenotype.tsv.zst
} -vars {dest build} -code {
	set target [gzroot $target].zst
	# get target2, uses biomart (ensembl geneset) and clinvar for gene-phenotype correlations
	cg download_phenotype $target2 ${build}
	# to reg file
	cg geneannot2reg $dep $target2 $target.temp
	compress $target.temp $target
}

job reg_${build}_go -deps {
	extra/reg_${build}_genes.tsv.zst
} -targets {
	extra/reg_${build}_go.tsv.zst
	extra/geneannot_${build}_go.tsv.zst
} -code {
	cg download_mart $target2 hsapiens_gene_ensembl gene_ensembl_config {external_gene_name go_id name_1006 namespace_1003 go_linkage_type}
	cg geneannot2reg $dep $target2 $target.temp
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

## exome variant server
#job var_${build}_evs -targets {
#	extra/var_${build}_evs.tsv.zt
#	extra/var_${build}_evs.tsv.opt
#	extra/var_${build}_evs.tsv.info
#} -vars {dest build db evsurl evsbuild} -code {
#	set target [gzroot $target].zst
#	if {$evsbuild eq $build} {
#		cg download_evs $target $evsurl
#	} else {
#		cg download_evs $target.$evsbuild.zst $evsurl
#		liftover_refdb $target.$evsbuild.zst $target $dest $evsbuild $build
#	}
#}

# exac
job var_${build}_exac -targets {
	extra/var_${build}_exac.tsv.zst
	extra/var_${build}_exac.tsv.opt
	extra/var_${build}_exac.tsv.info
} -vars {dest build db exacurl exacbuild} -code {
	set target [gzroot $target].zst
	if {$exacbuild eq $build} {
		cg download_exac --stack 1 --verbose 2 $target $exacurl
	} else {
		cg download_exac --stack 1 --verbose 2 $target.$exacbuild.zst $exacurl
		liftover_refdb $target.$exacbuild.zst $target $dest $exacbuild $build 0
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
	cg download_ucsc $target.ucsc $build wgEncodeRegDnaseClustered
	cg regcollapse $target.ucsc > $target.temp
	file rename -force -- $target.ucsc.info [gzroot $target].info
	file rename -force -- $target.temp $target
	file delete $target.ucsc
}

# todo : wgEncodeRegDnaseWig wgEncodeRegDnase

## link local data in dir
#foreach file [glob -nocomplain ../${build}-local/*] {
#	catch {
#		file delete extra/[file tail $file]
#		cplinked $file [file tail $file]
#	}
#}

# dbNSFP
job var_${build}_dbnsfp -targets {
	extra/var_${build}_dbnsfp.tsv.zst
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
	reg_${build}_ccr.tsv.zst
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

# alphamissense
# -------------
job alphamissense -deps {
} -vars {
	build dm_alphamissense_url dest
} -targets {
	var_${build}_alphamissense.tsv.zst
} -code {
	file_write [gzroot $target].info [subst [deindent {
		= alphamissense (alphafold2 based prediction of pathogenic missense variants) =
		adds am_pathogenicity (higher values are predicted to be more pathogenic) and am_class
		
		== Download info ==
		dbname	alphamissense
		version	no version given
		citation	Jun Cheng et al., Accurate proteome-wide missense variant effect prediction with AlphaMissense.Science381,eadg7492(2023).DOI:10.1126/science.adg7492
		license	CC BY-NC-SAL cite, non-commercial research use only
		source	$dm_alphamissense_url
		time	[timestamp]
		
		== Description ==
		
		More info on aphamissense in https://www.science.org/doi/10.1126/science.adg7492
		license and desciption of files in https://storage.cloud.google.com/dm_alphamissense/README.pdf
		
		Predictions for all possible missense variants for 60k non-canonical transcript isoforms
		(hg38, GENCODE V32). This file has transcript_id but no UniProt accession numbers.
		Predictions for non-canonical isoforms were not thoroughly evaluated and should be
		used with caution.
		
		== Category ==
		Annotation
	}]]
	file_write $target.opt "fields\t{am_pathogenicity am_class}\n"
	file mkdir $target.temp
	set tail [file tail $dm_alphamissense_url]
	wgetfile $dm_alphamissense_url $target.temp/$tail
	cg select -overwrite 1 -hc 1 -f {
		chromosome=$CHROM {begin=$POS - 1} end=$POS type="snp" ref=$REF alt=$ALT
		transcript_id protein_variant am_pathogenicity am_class
	} $target.temp/$tail | cg select -s - | cg collapsealleles | cg zst > $target.temp.zst
	file rename -force -- $target.temp.zst $target
	cg zstindex $target
}

job alphamissensecanonical -deps {
} -vars {
	build dm_alphamissense_canonical_url dest
} -targets {
	extra/var_${build}_alphamissensecanonical.tsv.zst
} -code {
	file_write [gzroot $target].info [subst [deindent {
		= alphamissense (alphafold2 based prediction of pathogenic missense variants) =
		adds am_pathogenicity (higher values are predicted to be more pathogenic) and am_class
		
		== Download info ==
		dbname	alphamissensecanonical
		version	no version given
		citation	Jun Cheng et al., Accurate proteome-wide missense variant effect prediction with AlphaMissense.Science381,eadg7492(2023).DOI:10.1126/science.adg7492
		license	CC BY-NC-SAL cite, non-commercial research use only
		source	$dm_alphamissense_canonical_url
		time	[timestamp]
		
		== Description ==
		
		More info on aphamissense in https://www.science.org/doi/10.1126/science.adg7492
		license and desciption of files in https://storage.cloud.google.com/dm_alphamissense/README.pdf
		
		Predictions for all possible single nucleotide missense variants (71M) from 19k human
		protein-coding genes (canonical transcripts).
		
		== Category ==
		Annotation
	}]]
	file_write $target.opt "fields\t{am_pathogenicity am_class}\n"
	file mkdir $target.temp
	set tail [file tail $dm_alphamissense_canonical_url]
	file mkdir $target.temp
	wgetfile $dm_alphamissense_canonical_url $target.temp/$tail
	cg select -overwrite 1 -hc 1 -f {
		chromosome=$CHROM {begin=$POS - 1} end=$POS ref=$REF alt=$ALT
		uniprot_id transcript_id protein_variant am_pathogenicity am_class
	} $target.temp/$tail | cg select -s - | cg collapsealleles | cg zst > $target.temp.zst
	file rename -force -- $target.temp.zst $target
	cg zstindex $target
	# clean
	file delete -force $target.temp
}

# gnomad lof
# ----------
job lofgnomad -deps {
} -vars {
	build gnomadlof gnomadlofbuild gnomadlofversion dest
} -targets {
	reg_${build}_lofgnomad.tsv.zst
} -code {
	file_write [gzroot $target].info [subst [deindent {
		= lofgnomad (loss of function tolerence) =
		
		== Download info ==
		dbname	lofgnomad
		version	$gnomadlofversion
		citation	Karczewski KJ, Francioli LC, Tiao G, Cummings BB, Alföldi J, Wang Q, Collins RL, Laricchia KM, Ganna A, Birnbaum DP, Gauthier LD, Brand H, Solomonson M, Watts NA, Rhodes D, Singer-Berk M, England EM, Seaby EG, Kosmicki JA, Walters RK, Tashman K, Farjoun Y, Banks E, Poterba T, Wang A, Seed C, Whiffin N, Chong JX, Samocha KE, Pierce-Hoffman E, Zappala Z, O'Donnell-Luria AH, Minikel EV, Weisburd B, Lek M, Ware JS, Vittal C, Armean IM, Bergelson L, Cibulskis K, Connolly KM, Covarrubias M, Donnelly S, Ferriera S, Gabriel S, Gentry J, Gupta N, Jeandet T, Kaplan D, Llanwarne C, Munshi R, Novod S, Petrillo N, Roazen D, Ruano-Rubio V, Saltzman A, Schleicher M, Soto J, Tibbetts K, Tolonen C, Wade G, Talkowski ME; Genome Aggregation Database Consortium; Neale BM, Daly MJ, MacArthur DG. The mutational constraint spectrum quantified from variation in 141,456 humans. Nature. 2020 May;581(7809):434-443. doi: 10.1038/s41586-020-2308-7. Epub 2020 May 27. Erratum in: Nature. 2021 Feb;590(7846):E53. Erratum in: Nature. 2021 Sep;597(7874):E3-E4. PMID: 32461654; PMCID: PMC7334197.
		citation	Lek M., Karczewski K., Exome Aggregation Consortium. Analysis of protein-coding genetic variation in 60,706 humans. Nature volume 536, pages 285-291 (2016)
		license	cite
		source	$gnomadlof
		time	[timestamp]
		
		== Description ==
		
		This region annotation assigns two gene constraint metrics (pLI and LOEUF) to all
		variants in the region covered by the gene (including intronic, etc.)

		pLI is the the probability of being loss-of-function (LoF) intolerant. It separates
		genes of sufficient length into LoF intolerant (pLI >= 0.9) or LoF tolerant (pLI <= 0.1) categories

		Low LOEUF (loss-of-function observed/expected upper bound fraction) values are indicative of strong intolerance.
		It is suggested to use LOEUF < 0.6 if a hard threshold is needed

		As different metrics are available for each gene (per transcript), the one given here is the maximum pLI
		and lowest LOEUF.

		More info on https://gnomad.broadinstitute.org/help/constraint
		and in https://www.nature.com/articles/nature19057
		
		== Category ==
		Annotation
	}]]
	file_write [gzroot $target].opt "fields\tpLI\tLOEUF\n"
	file mkdir $target.temp
	set tail [file tail $gnomadlof]
	wgetfile $gnomadlof $target.temp/$tail
	unset -nocomplain pLIa
	unset -nocomplain LOEUFa
	set f [gzopen $target.temp/$tail]
	set header [tsv_open $f]
	set poss [list_cor $header {gene gene_id lof.pLI lof.oe_ci.upper}]
	while 1 {
		if {[gets $f line] == -1} break
		foreach {gene geneid pLI LOEUF} [list_sub [split $line \t] $poss] break
		if {$pLI ne "NA" && (![info exists pLIa($gene)] || $pLIa($gene) < $pLI)} {
			set pLIa($gene) $pLI
		}
		if {$LOEUF ne "NA" && (![info exists LOEUFa($gene)] || $LOEUFa($gene) > $LOEUF)} {
			set LOEUFa($gene) $LOEUF
		}
	}
	gzclose $f
	set f [gzopen [gzfile gene_${build}_intGene.tsv]]
	set header [tsv_open $f]
	set poss [list_cor $header {chrom start end gene geneid name2}]
	unset -nocomplain genea
	set genes {}
	while 1 {
		if {[gets $f line] == -1} break
		foreach {chrom start end gene geneid name2} [list_sub [split $line \t] $poss] break
		if {[regexp {_alt$} $chrom] || [regexp {_fix$} $chrom] || [regexp {_random$} $chrom]} continue
		if {![info exists pLIa($gene)]} {
			set gene $name2
			if {![info exists pLIa($gene)]} {
				set gene $name2
				if {![info exists pLIa($gene)]} {
					set gene $geneid
					if {![info exists pLIa($gene)]} continue
				}
			}
		}
		if {![info exists genea($gene)]} {
			set genea($gene) [list $chrom $start $end $gene]
			lappend genes $gene
		} else {
			if {$chrom != [lindex $genea($gene) 0]} {
				error "gene $gene on 2 chromosomes ($chrom and [lindex $genea($gene) 0])"
			}
			if {$start < [lindex $genea($gene) 1]} {
				lset genea($gene) 1 $start
			}
			if {$end > [lindex $genea($gene) 2]} {
				lset genea($gene) 2 $end
			}
		}
		
	}
	gzclose $f
	set o [wgzopen $target.temp/reg_${build}_lofgnomad.tsv.zst]
	puts $o [join {chromosome begin end gene pLI LOEUF} \t]
	foreach gene $genes {
		set resultline $genea($gene)
		lappend resultline $pLIa($gene) $LOEUFa($gene)
		puts $o [join $resultline \t]
	}
	gzclose $o
	cg select -s - $target.temp/reg_${build}_lofgnomad.tsv.zst $target.temp/reg_${build}_lofgnomad2.tsv.zst
	file rename -force $target.temp/reg_${build}_lofgnomad2.tsv.zst reg_${build}_lofgnomad.tsv.zst
	cg zstindex reg_${build}_lofgnomad.tsv.zst
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
	# homfreqp is percentage of sindividuals that are homozygous -> so should divide count by (AN_popmax/2), do this by *200.0 instead of *100.0
	foreach population {
		afr amr asj eas fin nfe sas oth nfe_nwe nfe_seu nfe_est nfe_bgr nfe_swe
		male female
	} {
		if {![inlist $header AN_$population]} continue
		if {[lsearch $header AC_${population}] == -1} {
			error "field AC_${population} not found in $file"
		}
		if {[lsearch $header nhomalt_${population}] == -1} {
			error "field AC_${population} not found in $file"
		}
		lappend fields "${population}_freqp=if(def(\$AN_$population,0) < 8, \"-\", format(\"%.3f\",(100.0 * \$AC_${population})/\$AN_$population))"
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
	set maxpops {afr amr eas nfe sas}
	lappend fields "max_freqp=lmaxd(\$[join $maxpops {_freqp,$}]_freqp,\"-\")"
	lappend fields "max_homfreqp=lmaxd(\$[join $maxpops {_homfreqp,$}]_homfreqp,\"-\")"
	return $fields
}

set deps {}
foreach chromosome {
	1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
} {
	puts chr$chromosome
	set vcf gnomad.genomes.v$gnomadversion.sites.chr$chromosome.vcf.bgz
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
			wgetfile $gnomadbaseurl/vcf/genomes/$vcf $tempdir/$vcf
		}
		set fields [gnomadfields $tempdir/$vcf]
		putslog "Converting $vcf"
		cg vcf2tsv -split 1 -sort 0 $tempdir/$vcf | cg select --stack 1 -rc 1 -f $fields | cg collapsealleles {*}[compresspipe $target 1] > $target.temp
		file rename -- $target.temp $target
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
		file rename -- [gzroot $finaltarget].info $tempdir/result.tsv.temp.info
		liftover_refdb $tempdir/result.tsv.temp.zst $finaltarget $dest $gnomadbuild $build 0
	}
	cg zindex $finaltarget
	cg index $finaltarget
	# file delete -force $tempdir
}

# gnomadex
# --------
set finaltarget var_${build}_gnomadex.tsv.zst
set tempdir $finaltarget.temp
file mkdir $tempdir

set deps {}
foreach chromosome {
	1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
} {
	puts chr$chromosome
	set vcf gnomad.exomes.v$gnomadversion.sites.chr$chromosome.vcf.bgz
	set target $tempdir/result$chromosome.tsv.zst
	lappend deps $target
	job var_${build}_gnomadex-$chromosome -targets {
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
			wgetfile $gnomadbaseurl/vcf/exomes/$vcf $tempdir/$vcf
		}
		set fields [gnomadfields $tempdir/$vcf]
		putslog "Converting $vcf"
		cg vcf2tsv -split 1 -sort 0 $tempdir/$vcf | cg select --stack 1 -rc 1 -f $fields | cg collapsealleles {*}[compresspipe $target 1] > $target.temp
		file rename -- $target.temp $target
	}
}

job var_${build}_gnomadex-final -deps $deps -targets {
	$finaltarget
	var_${build}_gnomad.tsv.opt
} -vars {tempdir dest db build finaltarget gnomadbuild} -code {
	file_write var_${build}_gnomad.tsv.opt "fields\t{max_freqp nfe_freqp}\n"
	if {$gnomadbuild eq $build} {
		exec cg cat {*}$deps {*}[compresspipe $finaltarget 12] > $tempdir/result.tsv.temp.zst
		file rename -force -- $tempdir/result.tsv.temp.zst $finaltarget
	} else {
		exec cg cat {*}$deps {*}[compresspipe $finaltarget 1] > $tempdir/result.tsv.temp.zst
		file rename -- [gzroot $finaltarget].info $tempdir/result.tsv.temp.info
		liftover_refdb $tempdir/result.tsv.temp.zst $finaltarget $dest $gnomadbuild $build 0
	}
	cg zindex $finaltarget
	cg index $finaltarget
	# file delete -force $tempdir
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

job svgnomad -deps {
} -targets {
	sv_${build}_svgnomad.tsv.zst
} -vars {
	build gnomadsv gnomadsvbuild gnomadsvversion dest
} -procs {
	gnomadfields
} -code {
	file_write [gzroot $target].info [subst [deindent {
		= svgnomad (loss of function tolerence) =
		
		== Download info ==
		dbname	svgnomad
		version	$gnomadsvversion
		citation	Collins RL, Brand H, Karczewski KJ, Zhao X, Alföldi J, Francioli LC, Khera AV, Lowther C, Gauthier LD, Wang H, Watts NA, Solomonson M, O'Donnell-Luria A, Baumann A, Munshi R, Walker M, Whelan CW, Huang Y, Brookings T, Sharpe T, Stone MR, Valkanas E, Fu J, Tiao G, Laricchia KM, Ruano-Rubio V, Stevens C, Gupta N, Cusick C, Margolin L; Genome Aggregation Database Production Team; Genome Aggregation Database Consortium; Taylor KD, Lin HJ, Rich SS, Post WS, Chen YI, Rotter JI, Nusbaum C, Philippakis A, Lander E, Gabriel S, Neale BM, Kathiresan S, Daly MJ, Banks E, MacArthur DG, Talkowski ME. A structural variation reference for medical and population genetics. Nature. 2020 May;581(7809):444-451. doi: 10.1038/s41586-020-2287-8. Epub 2020 May 27. Erratum in: Nature. 2021 Feb;590(7846):E55. PMID: 32461652; PMCID: PMC7334194
		license	cite
		source	$gnomadsv
		time	[timestamp]
		
		== Description ==
		
		gnomad structural variants
		
		More info on https://gnomad.broadinstitute.org/faq
		and in https://doi.org/10.1038/s41586-020-2287-8
		
		== Category ==
		Annotation
	}]]
	file_write [gzroot $target].opt "fields\t{max_freqp eur_freqp}\n"
	file mkdir $target.temp
	set tail [file tail $gnomadsv]
	wgetfile $gnomadsv $target.temp/$tail
	set root [file root [gzroot $tail]]
	#
	# convert to tsv
	set temptsv $target.temp/$root.tsv.zst
	file delete $target.temp/$root.tsv.zst
	cg vcf2tsv -split 0 $target.temp/$tail $temptsv
	#
	# change fields
	set temptarget $target.temp/$root-2.tsv.zst
	set header [cg select -header $temptsv]
	file delete $temptarget
	catch {gzclose $o} ; catch {gzclose $f}
	set o [wgzopen $temptarget]
	set f [gzopen $temptsv]
	while {[gets $f line] != -1} {
		if {[regexp {^#fields} $line]} break
		if {![regexp {^#} $line]} break
		puts $o $line
	}
	puts $o [deindent {
		#fields	table
		#fields	field	number	type	description	source
		#fields	chromosome	1	String	Chromosome/Contig	var
		#fields	begin	1	Integer	Begin of feature (0 based - half open)	var
		#fields	end	1	Integer	End of feature (0 based - half open)	var
		#fields	type	1	String	Type of feature (snp,del,ins,...)	var
		#fields	ref	1	String	Reference sequence, can be a number for large features	var
		#fields	alt	1	String	Alternative sequence, can be a number for large features	var
		#fields	name	1	String	name of feature	var
		#fields	CPX_INTERVALS	.	String	Genomic intervals constituting complex variant.	info
		#fields	CPX_TYPE	1	String	Class of complex variant.	info
	}]
	set fields {chromosome begin end type ref alt CPX_INTERVALS CPX_TYPE name}
	foreach population {
		afr amr eas fin mid nfe rmi sas
		XY XX
	} {
		set upopulation [string toupper $population]
		if {![inlist $header AN_${population}]} continue
		set POP [string toupper $population]
		lappend fields "${population}_an=\$AN_${population}"
		puts $o "\#fields\t${population}_an\t1\tInteger\tTotal number of $POP alleles genotyped (biallelic sites only)"
		lappend fields "${population}_ac=\$AC_${population}"
		puts $o "\#fields\t${population}_ac\t1\tInteger\tNumber of non-reference $POP observed genotyped (for biallelic sites)"
		lappend fields "${population}_freqp=if(def(\$AN_${population},0) < 8, \"-\", vformat(\"%.4f\",(100.0 @* vdef(\$AF_${population},0))))"
		puts $o "\#fields\t${population}_freqp\tA\tFloat\t$POP allele frequency in percent (for biallelic sites)"
		if {[inlist $header FREQ_HOMALT_${population}]} {
			lappend fields "${population}_homfreqp=if(def(\$AN_${population},0) < 8, \"-\", vformat(\"%.4f\",(100.0 @* vdef(\$FREQ_HOMALT_${population},0))))"
			puts $o "\#fields\t${population}_homfreqp\tA\tFloat\t$POP homozygous alternate genotype frequency in percent (biallelic sites only)"
		}
	}
	lappend fields "max_freqp=vformat(\"%.4f\",(100.0 @* vdef(\$GRPMAX_AF,0)))"
	puts $o "\#fields\t${population}_homfreqp\tA\tFloat\t$POP homozygous alternate genotype frequency in percent (biallelic sites only)"
	while {[regexp {^#} $line]} {
		if {![regexp {^#fields} $line]} {
			puts $o $line
		}
		if {[gets $f line] == -1} break
	}
	gzclose $f
	exec cg select -rc 1 -f $fields $temptsv >@ $o
	gzclose $o
	file delete $target
	if {$build ne $gnomadsvbuild} {
		liftover_refdb $temptarget $target $dest $gnomadsvbuild $build 0
	} else {
		file rename -force -- $temptarget $target
	}
	cg zstindex $target
	file delete -force $target.temp
}

# GATK resources
job gatkres_${build} -targets {
	gatkres
} -vars {
	dest build
} -code {
	mkdir ${dest}/${build}/gatkres.temp
	cd ${dest}/${build}/gatkres.temp
	foreach file {
		1000G_phase1.snps.high_confidence.hg38.vcf.gz
		1000G_omni2.5.hg38.vcf.gz
		hapmap_3.3.hg38.vcf.gz
		Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
		Homo_sapiens_assembly38.dbsnp138.vcf
	} {
		wgetfile https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/$file
		catch {wgetfile https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/$file.tbi}
	}
	cg gzip Homo_sapiens_assembly38.dbsnp138.vcf
	cg maketabix Homo_sapiens_assembly38.dbsnp138.vcf.gz
	file rename ${dest}/${build}/gatkres.temp ${dest}/${build}/gatkres
	cd ${dest}/${build}
}
cd ${dest}/${build}

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
				file rename -- $tempdir/collapsed.tsv.zst.temp $tempdir/collapsed.tsv.zst
			}
			cg liftover --stack 1 -split 0 -s 0 $tempdir/collapsed.tsv.zst ../liftover/${caddbuild}To${build}.over.tsv {*}[compresspipe .zst] > $tempdir/collapsed${build}.tsv.temp
			file rename -- $tempdir/collapsed${build}.tsv.temp $tempdir/collapsed${build}.tsv.zst
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
	puts "Transfering $file"
	set tail [file tail $file]
	set newtail [join [lreplace [split $tail _] 1 1 $build] _]
	if {[file exists extra/$newtail]} {
		puts "skipping extra/$newtail: already exists"
		continue
	}
	set filebuild [lindex [split $tail _] 1]
	hardcopy $file extra/
	catch {hardcopy $file.zsti extra/}
	if {$filebuild ne $build} {
		file delete extra/$newtail
		liftover_refdb extra/$tail extra/$newtail $dest $filebuild $build
	}
}

job_wait

