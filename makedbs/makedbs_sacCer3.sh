#!/bin/sh
# the next line restarts using tclsh \
exec cg source "$0" "$@"

# Script to download and prepare all reference databases (including genome) for sacCer3
# call using
# makedbs_sacCer3.sh ?options? ?dest? ?webcache?
# where 
# * dest: the directory where the reference databases will be installed (default /complgen/refseqnew/sacCer3)
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
set build sacCer3
set defaultdest /complgen/refseqnew

set genomeurl {}
set par {}
set dbsnpversion {}
# set refSeqFuncElemsurl https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/109.20200522/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz
set refSeqFuncElemsurl https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/110/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz
set mirbase {}
set regionsdb_collapse {
	cytoBand microsat rmsk simpleRepeat oreganno
	sgdOther
	phastConsElements7way
}
set regionsdb_join {
}

# list with geneset name (first word) and one or more of the following keywords
# int : include in intGene
# extra : place in the extra dir instead of in base annotation dir
# reg : make a region file from it (in extra)
set genesdb [list \
	{refGene int reg} \
	{ncbiRefSeq extra int reg} \
	{sgdGene extra reg} \
	{ensGene extra int reg} \
	{augustusGene extra} \
]

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
cg_options makedbs_sacCer3.sh argv {
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
	-regionsdb_collapse $regionsdb_collapse \
	-regionsdb_join $regionsdb_join \
	-dbsnp $dbsnpversion \
	-refSeqFuncElemsurl $refSeqFuncElemsurl \
	-genesdb $genesdb \
	-mirbase $mirbase \
	-webcache $webcache \
	$dest/$build

# rest after this is sacCer3 specific code
# -------------------------------------

job_wait

