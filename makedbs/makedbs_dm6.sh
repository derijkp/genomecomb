#!/bin/sh
# the next line restarts using tclsh \
exec cg source "$0" "$@"

# Script to download and prepare all reference databases (including genome) for dm6
# call using
# makedbs_dm6.sh ?options? ?dest? ?webcache?
# where 
# * dest: the directory where the reference databases will be installed (default /complgen/refseqnew/dm6)
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
set build dm6
set defaultdest /complgen/refseqnew
set par {}
set genomeurl {}
set dbsnpversion {}
set refSeqFuncElemsurl https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Drosophila_melanogaster/all_assembly_versions/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff.gz
set mirbase dme-22.1:dm6
set regionsdb_collapse {
	cytoBand rmsk simpleRepeat microsat
}
set regionsdb_join {}

# list with geneset name (first word) and one or more of the following keywords
# int : include in intGene
# extra : place in the extra dir instead of in base annotation dir
# reg : make a region file from it (in extra)
set genesdb [list \
	{refGene int reg} \
	{ncbiRefSeq extra int reg} \
	{ensGene extra int reg} \
	{genscan extra} \
	{augustusGene extra} \
]
set downloads {
	https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/genes/dm6.ensGene.gtf.gz gene_dm6_ensGene.gtf.gz
}
set transcriptsurl https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/genes/dm6.ensGene.gtf.gz
set transcriptsgtf extra/gene_dm6_ensGene.gtf

# extra settings
# --------------

# prepare
# =======

if {![info exists argv]} {set argv {}}

# keep actual command line used for log
set cmdline [clean_cmdline cg [info script] {*}$args]

# arguments, start job system
logverbose 2
set argv [job_init {*}$argv]

cg_options makedbs_mm10.sh argv {
} {dest webcache} 0 2

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
	-dbsnp $dbsnpversion \
	-regionsdb_collapse $regionsdb_collapse \
	-regionsdb_join $regionsdb_join \
	-refSeqFuncElemsurl $refSeqFuncElemsurl \
	-genesdb $genesdb \
	-mirbase $mirbase \
	-transcriptsurl $transcriptsurl \
	-transcriptsgtf $transcriptsgtf \
	-webcache $webcache \
	$dest/$build

# rest after this is dm6 specific code
# ------------------------------------

job_wait
