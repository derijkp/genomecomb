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

# basic settings (for makerefdb_job)
# --------------

set build mm10
set defaultdest /complgen/refseqnew

# basic settings
# use default (ucsc goldenPath) genome source
set genomeurl {}
set par {chromosome	begin	end	name
X	169969759	170931299	PAR1
Y	90745845	91644698	PAR1
}
set organelles {chromosome
chrM
}
set dbsnpversion {}
set refSeqFuncElemsurl https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/annotation_releases/108.20200622/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.gff.gz
set mirbase mmu-22.1:mm10
set regionsdb_collapse {
	cytoBand oreganno cpgIslandExt tRNAs
	microsat rmsk simpleRepeat 
	multiz60way phastConsElements60way phastConsElements60wayPlacental phastConsElements60wayEuarchontoGlires
}
set regionsdb_join {
	genomicSuperDups
}

set gencodeversion 25
# list with geneset name (first word) and one or more of the following keywords
# int : include in intGene
# extra : place in the extra dir instead of in base annotation dir
# reg : make a region file from it (in extra)
set genesdb [list \
	{refGene int reg} \
	[list wgEncodeGencodeBasicVM${gencodeversion} gencode extra int reg] \
	{knownGene extra int reg} \
	[list wgEncodeGencodeCompVM${gencodeversion} cgencode extra] \
	{genscan extra} \
]

# extra settings
# --------------
set transcriptsurl http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
set transcriptsgtf extra/gene_mm10_gencode.vM25.gtf

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
	-organelles $organelles \
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

# rest after this is mm10 specific code
# -------------------------------------

set evasnpurl https://hgdownload.soe.ucsc.edu/gbdb/$build/bbi/evaSnp5.bb
job var_${build}evaSnp -deps {
} -targets {
	var_${build}_evaSnp.tsv.zst
} -vars {
	build evasnpurl
} -code {
	if {![file exists ${build}_evaSnp5.bb.bed]} {
		cg_download_bigbed $evasnpurl ${build}_evaSnp5.bb.temp.bed
		file rename ${build}_evaSnp5.bb.temp.bed ${build}_evaSnp5.bb.bed
	}
	set f [open ${build}_evaSnp5.bb.bed]
	set o [wgzopen ${build}_evaSnp5.tsv.temp.zst]
	puts $o [join {chromosome begin end type ref alt name} \t]
	while {[gets $f line] != -1} {
		foreach {chromosome begin end ref alt name} [list_sub [split $line \t] {0 1 2 9 10 3}] break
		# putsvars chromosome begin end ref alt name
		set l1 [string length $ref]
		set l2 [string length $alt]
		if {$l1 > 1} {
			if {$l2 > 1} {
				set type sub
			} elseif {[string index $ref 0] ne [string index $alt 0]} {
				set type sub
			} else {
				set type del
				set ref [string range $ref 1 end]
				set alt [string range $alt 1 end]
			}
		} elseif {$l2 > 1} {
			if {[string index $ref 0] ne [string index $alt 0]} {
				set type sub
			} else {
				set type ins
				set ref [string range $ref 1 end]
				set alt [string range $alt 1 end]
			}
		} else {
			set type snp
		}
		puts $o $chromosome\t$begin\t$end\t$type\t$ref\t$alt\t$name
	}
	gzclose $o
	close $f
	cg select -s - ${build}_evaSnp5.tsv.temp.zst ${build}_evaSnp5.tsv.temp2.zst
	file rename -force ${build}_evaSnp5.tsv.temp2.zst	var_${build}_evaSnp.tsv.zst
	file delete ${build}_evaSnp5.bb.bed ${build}_evaSnp5.tsv.temp.zst
}

job_wait
