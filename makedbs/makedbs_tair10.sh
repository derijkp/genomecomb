#!/bin/sh
# the next line restarts using tclsh \
exec cg source "$0" "$@"

# Script to download and prepare all reference databases (including genome) for mm10
# call using
# makedbs_mm10.sh ?options? ?dest? ?webcache?
# where 
# * dest: the directory where the reference databases will be installed (default /complgen/refseqnew/mm10)
# * webcache: directory that will be used to cache downloads (default /complgen/refseqnew/webcache)
# options:
# -d distr: allow distributed processing on 
#           e.g. a grid engine cluster with -d sge, 
#           or locally on multiple cores with e.g. -d 8 for 8 cores
# If the command is interrupted, you can just restart using the same
# command, and it will continue from where it was stopped.

# settings
# ========

# basic settings
# --------------

set build tair10
set defaultdest /complgen/refseqnew
# should first try Arabidopsis_thaliana.TAIR10.dna_sm.primary_assembly.fa.gz
# toplevel here because there are no alt/hap
set source https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-56/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz
set genomeseq genome_tair10.ifas

set organelles {chromosome
Mt
Pt
}


# prepare
# =======

if {![info exists argv]} {set argv {}}

# keep actual command line used for log
set cmdline [clean_cmdline cg [info script] {*}$args]

# arguments, start job system
logverbose 2
set argv [job_init {*}$argv]
cg_options makedbs_tair10.sh argv {
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

exec wget $source >@ stdout 2>@ stderr
cg fas2ifas Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz $genomeseq
file delete Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz
file_write $genomeseq.info [subst [deindent {
	= Genome (build $build) =
	
	== Download info ==
	dbname	genome
	version	[timestamp]
	license	free
	source	$source
	time	[timestamp]
	
	== Description ==
	The genome sequence downloaded from source, and chromosome sorted according 
	to a natural sort (for use in genomecomb)
	
	== Category ==
	Genome
}]]

# make samtools index
exec samtools faidx $genomeseq
#
set extra [file dir $genomeseq]/extra
mkdir $extra

file_write extra/reg_${build}_organelles.tsv $organelles

set rfile $extra/reg_${build}_fullgenome.tsv
putslog "Making $rfile"
set data [file_read $genomeseq.fai]
set o [open $rfile.temp w]
puts $o chromosome\tbegin\tend
list_foreach {chromosome len} [lrange [split [string trim $data] \n] 0 end] {
	if {$chromosome eq ""} continue
	puts $o $chromosome\t0\t$len
}
close $o
file rename -force -- $rfile.temp $rfile

job genome_${build}_forcram -deps {
	genome_${build}.ifas
} -targets {
	genome_${build}.ifas.forcram
} -code {
	cg fasta2cramref $dep $target
}

job reg_${build}_sequencedgenome -deps {
	genome_${build}.ifas
} -targets {
	extra/reg_${build}_sequencedgenome.tsv.zst
} -vars {dest build} -code {
	exec cg calcsequencedgenome --stack 1 $dep {*}[compresspipe $target 12] > $target.temp
	file rename -force -- $target.temp $target
}

distrreg_nolowgene $dbdir 100000 100k

# make minimap2 versions of genome
refseq_minimap2_job genome_${build}.ifas sr
refseq_minimap2_job genome_${build}.ifas map-ont
refseq_minimap2_job genome_${build}.ifas splice

# homopolymer
job reg_${build}_homopolymer -deps {
	genome_${build}.ifas
} -targets {
	reg_${build}_homopolymer.tsv.zst
	reg_${build}_homopolymer.tsv.gz
	reg_${build}_homopolymer.tsv.gz.tbi
	reg_${build}_homopolymer.tsv.opt
} -vars {dest build db} -code {
	set target reg_${build}_homopolymer.tsv.zst
	file_write [gzroot $target].opt "fields\t{base size}\n"
	cg extracthomopolymers genome_${build}.ifas {*}[compresspipe $target 12] > $target.temp
	file rename -force -- $target.temp $target
        cg maketabix $target
}

job_wait
