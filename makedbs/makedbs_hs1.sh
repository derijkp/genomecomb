#!/bin/sh
# the next line restarts using tclsh \
exec cg source "$0" "$@"

# Script to download and prepare all reference databases (including genome) for hs1
# call using
# makedbs_hs1.sh ?options? ?dest? ?webcache?
# where 
# * dest: the directory where the reference databases will be installed (default /complgen/refseqnew/hs1)
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
set build hs1
set defaultdest /complgen/refseqnew

set genomeurl {}
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
set mirbase hsa-22.1:hg38
set regionsdb_collapse {
	cytoBand microsat oreganno rmsk simpleRepeat 
	tRNAs wgRna
	phastConsElements100way phastConsElements30way
}
set regionsdb_join {
	chainSelf dgvMerged genomicSuperDups
}

set gencodeversion 42
# list with geneset name (first word) and one or more of the following keywords
# int : include in intGene
# extra : place in the extra dir instead of in base annotation dir
# reg : make a region file from it (in extra)
set genesdb [list \
	{hub_567047_ncbiRefSeqCurated int reg} \
	{hub_567047_ncbiRefSeq extra int reg} \
	{hub_567047_catLiftOffGenesV1 int} \
]

# extra settings
# --------------
set 1000g3url http://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz
set 1000g3readmeurl http://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/README_phase3_callset_20150220
set 1000g3build hg19
set clinvarurl https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20221119.vcf.gz
set clinvarpapuurl https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20221119_papu.vcf.gz
#set kaviarurl http://s3-us-west-2.amazonaws.com/kaviar-160204-public/Kaviar-160204-Public-${build}-trim.vcf.tar
#set kaviarbuild hg19
set evsurl http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.protein-hgvs-update.snps_indels.vcf.tar.gz
set evsbuild hg19
set exacurl ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/ExAC.r1.sites.vep.vcf.gz
set exacbuild hg19
set caddversion 1.6
set caddurl http://krishna.gs.washington.edu/download/CADD/v$caddversion/GRCh38/whole_genome_SNVs.tsv.gz
set caddbuild hg38
set gnomadbuild hg38
set gnomadversion 3.1.2
set gnomadbaseurl https://storage.googleapis.com/gcp-public-data--gnomad/release/$gnomadversion/vcf
set gnomadexbuild hg19
set gnomadexversion 2.1.1
set gnomadexurl https://storage.googleapis.com/gcp-public-data--gnomad/release/$gnomadexversion/vcf/exomes/gnomad.exomes.r$gnomadexversion.sites.vcf.bgz
set gnomadlofbuild hg19
set gnomadlofversion 2.1.1
set gnomadlof https://storage.googleapis.com/gcp-public-data--gnomad/release/$gnomadlofversion/constraint/gnomad.v$gnomadlofversion.lof_metrics.by_gene.txt.bgz
set gnomadsvbuild hg19
set gnomadsvversion 2.1
set gnomadsv https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz
set ccrversion 2.20180420
set ccrurl https://s3.us-east-2.amazonaws.com/ccrs/ccrs/ccrs.autosomes.v${ccrversion}.bed.gz
set ccrbuild hg19
# set dbnsfpurl ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFP4.3a.zip
set dbnsfpurl https://dbnsfp.s3.amazonaws.com/dbNSFP4.3a.zip
set dbnsfpbuild hg38
set transcriptsurl http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.annotation.gtf.gz
set transcriptsgtf extra/gene_hs1_gencode.v42.gtf

# prepare
# =======

if {![info exists argv]} {set argv {}}

# keep actual command line used for log
set cmdline [clean_cmdline cg [info script] {*}$args]

# arguments, start job system
logverbose 2
set argv [job_init {*}$argv]
cg_options makedbs_hs1.sh argv {
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
	-genesdb $genesdb \
	-mirbase $mirbase \
	-webcache $webcache \
	$dest/$build

proc bed2genetsv {}

wgetfile http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/bigBedToBed
exec chmod u+x ./bigBedToBed

set genesets {}
foreach {dir db} {
	catLiftOffGenesV1 catLiftOffGenesV1
	ncbiRefSeq ncbiRefSeqCurated
	ncbiRefSeq ncbiRefSeq
} {
	lappend genesets gene_hs1_$db.tsv.zst
	job gene_$db -vars {
		dir db
	} -targets {
		gene_hs1_$db.tsv.zst
	} -code {
		wgetfile http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/bigBedToBed
		exec chmod u+x ./bigBedToBed
		wgetfile https://hgdownload.soe.ucsc.edu/gbdb/hs1/$dir/$db.bb
		exec ./bigBedToBed $db.bb $db.bed
		cg bed2genetsv $db.bed temp_gene_hs1_$db.tsv.zst
		cg select -s - temp_gene_hs1_$db.tsv.zst temp2_gene_hs1_$db.tsv.zst
		file rename -force temp2_gene_hs1_$db.tsv.zst gene_hs1_$db.tsv.zst
		file delete $db.bb $db.bed temp_gene_hs1_$db.tsv.zst
	}
}

set target gene_${build}_intGene.tsv.zst
job intgene -deps $genesets -targets {
	$target
} -vars {
	genesets
} -code {
	cg intgene {*}$genesets {*}[compresspipe $target 12] > $target.temp
	file rename $target.temp $target
}

set target var_${build}_clinvar.tsv.zst
job clinvar -targets {
	$target
} -vars {
} -code {
	wgetfile https://hgdownload.soe.ucsc.edu/gbdb/hs1/clinVar20220313/chm13v2.0_ClinVar20220313.vcf.gz
	cg vcf2tsv chm13v2.0_ClinVar20220313.vcf.gz temp_$target
	file rename temp_$target $target
	file delete chm13v2.0_ClinVar20220313.vcf.gz
}

set target reg_${build}_rmsk.tsv.zst
job rmsk -targets {
	$target
} -vars {
} -code {
	wgetfile https://hgdownload.soe.ucsc.edu/gbdb/hs1/t2tRepeatMasker/chm13v2.0_rmsk.bb
	exec ./bigBedToBed chm13v2.0_rmsk.bb chm13v2.0_rmsk.bed
	cg bed2tsv chm13v2.0_rmsk.bed temp_reg_hs1_rmsk.tsv.zst
	cg select -s - temp_reg_hs1_rmsk.tsv.zst temp2_reg_hs1_rmsk.tsv.zst
	file rename -force temp2_reg_hs1_rmsk.tsv.zst reg_hs1_rmsk.tsv.zst
	file delete chm13v2.0_rmsk.bed temp_reg_hs1_rmsk.tsv.zst
}

set target var_${build}_snp.tsv.zst
job snp -targets {
	$target
} -vars {
} -code {
	wgetfile https://hgdownload.soe.ucsc.edu/gbdb/hs1/dbSNP155/chm13v2.0_dbSNPv155.vcf.gz
	cg vcf2tsv chm13v2.0_dbSNPv155.vcf.gz temp_$target
	file rename temp_$target $target
	file delete chm13v2.0_dbSNPv155.vcf.gz
}

job_wait

