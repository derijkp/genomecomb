#!/bin/sh
# the next line restarts using tclsh \
exec cg source "$0" "$@"

# settings
# ========

set defaultdest /complgen/refseqnew
set build susScr3
set par {}
set dbsnpversion 138
set mirbase ssc-21:hg38
set genesdb [list {refGene int} \
	{ensGene int} \
	{genscan extra} \
]

# prepare
# =======

# keep actual command line used for log
set cmdline "[list cd [pwd]] \; [list [info script] {*}$argv]"

# arguments, start job system
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

# set 
logverbose 2
set_job_logdir log_jobs
job_logfile ${dest}/${build}/log_makedbs_${build} ${dest}/${build} $cmdline

# download
# ========
#
# actual downloads are done using makerefdb_job, which does the same as "cg makerefdb ..."
# but jobs are run under the previously started jobmanager and job settings

makerefdb_job \
	-genomeurl $genomeurl \
	-refSeqFuncElemsurl $refSeqFuncElemsurl \
	-genesdb $genesdb \
	-pseudoautosomal $par \
	-dbsnp $dbsnpversion \
	-mirbase $mirbase \
	-webcache $webcache \
	$dest/$build >@ stdout 2>@ stderr

job_wait
