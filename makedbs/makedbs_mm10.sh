#!/bin/sh
# the next line restarts using tclsh \
exec cg source "$0" "$@"

set build mm10
set defaultdest /complgen/refseqnew

# settings
set mirbasegenome mmu
set mirbaserelease 22
set mirbasebuild mm10
set dbsnpversion 142
set gencodeversion 18

# arguments
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

logverbose 2
job_logdir log_jobs

# download
# ========
#

# readme
set c [file_read $genomecombdir/docs/dbdir_README.txt]
regsub {version: [0-9.]+} $c "version: 0.99\ntime: [lindex [timestamp] 0]" c
file_write README_dbdir.txt $c

# download genome
job genome_${build} -targets {
	genome_${build}.ifas
	genome_${build}.ifas.fai
	extra/reg_${build}_fullgenome.tsv
} -vars build -code {
	cg download_genome -alt 0 genome_${build}.ifas ${build} 2>@ stderr
	file rename -force reg_genome_${build}.tsv extra/reg_${build}_fullgenome.tsv
	cg lz4 -i 1 extra/reg_${build}_fullgenome.tsv
}

job genome_${build}_cindex -deps {
	genome_${build}.ifas
} -targets {
	genome_${build}.ssa
} -code {
	cg make_genomecindex $dep
}

job reg_${build}_sequencedgenome -deps {
	genome_${build}.ifas
} -targets {
	extra/reg_${build}_sequencedgenome.tsv.lz4
} -vars {dest build} -code {
	exec cg calcsequencedgenome --stack 1 $dep | lz4c -12 > $target.temp
	file rename -force $target.temp $target
}

# make bwa version of genome
bwarefseq_job genome_${build}.ifas

# region databases (ucsc)
# you can explicitely download info on a database using:
# cg download_ucscinfo resultfile ${build} dbname



# collapse regions
foreach db {
	cytoBand oreganno cpgIslandExt tRNAs
	microsat rmsk simpleRepeat 
	multiz60way phastConsElements60way phastConsElements60wayPlacental phastConsElements60wayEuarchontoGlires
} {
	job reg_${build}_$db -targets {
		reg_${build}_${db}.tsv
	} -vars {dest build db} -code {
		set target [gzroot $target].lz4
		cg download_ucsc $target.ucsc ${build} $db
		cg regcollapse $target.ucsc > $target.temp
		file delete $target.ucsc
		cg lz4 -i 1 $target.temp
		file rename -force $target.ucsc.info [gzroot $target].info
		file rename -force $target.temp.lz4 $target
		file rename -force $target.temp.lz4.lz4i $target.lz4i
	}
}

# join regions
foreach db {
	genomicSuperDups
} {
	job reg_${build}_$db -targets {
		reg_${build}_${db}.tsv
	} -vars {dest build db} -code {
		set target [gzroot $target].lz4
		cg download_ucsc $target.ucsc ${build} $db
		cg regjoin $target.ucsc > $target.temp
		file delete $target.ucsc
		cg lz4 -i 1 $target.temp
		file rename -force $target.ucsc.info [gzroot $target].info
		file rename -force $target.temp.lz4 $target
		file rename -force $target.temp.lz4.lz4i $target.lz4i
	}
}

foreach db {
	rmsk simpleRepeat
} {
	job maketabix_${build}_$db -deps {
		reg_${build}_${db}.tsv
	} -targets {
		reg_${build}_${db}.tsv.gz.tbi
		reg_${build}_${db}.tsv.gz
	} -vars {build db} -code {
		cg maketabix $dep
	}
}

# dbsnp
job dbsnp -targets {
	var_${build}_dbsnp.tsv
	var_${build}_dbsnp.tsv.opt
} -vars {dest build dbsnpversion} -code {
	set target [gzroot $target].lz4
	file_write [gzroot $target].opt "fields\t{name}\n"
	cg download_dbsnp $target ${build} snp$dbsnpversion 2>@ stderr
	cg lz4index $target
}

job dbsnpCommon -targets {
	var_${build}_dbsnpCommon.tsv
} -vars {dest build dbsnpversion} -code {
	set target [gzroot $target].lz4
	file_write [gzroot $target].opt "fields\t{freqp}\n"
	cg download_dbsnp $target ${build} snp${dbsnpversion}Common 2>@ stderr
	cg lz4index $target
}

foreach db [list dbsnp dbsnpCommon] {
	job maketabix_${build}_$db -deps {
		var_${build}_${db}.tsv
	} -targets {
		var_${build}_${db}.tsv.gz.tbi
		var_${build}_${db}.tsv.gz
	} -vars {dest build db} -code {
		cg maketabix $dep
	}
}

# genes
foreach db [list \
	refGene wgEncodeGencodeBasicVM${gencodeversion} wgEncodeGencodeCompVM${gencodeversion} \
	genscan \
] {
	if {$db eq "wgEncodeGencodeCompV19"} {
		set dbname gencode
	} elseif {$db eq "wgEncodeGencodeBasicVM${gencodeversion}"} {
		set dbname gencode
	} elseif {$db eq "wgEncodeGencodeCompVM${gencodeversion}"} {
		set dbname cgencode
	} elseif {$db eq "lincRNAsTranscripts"} {
		set dbname lincRNA
	} else {set dbname $db}
	if {$db in "refGene lincRNAsTranscripts"} {
		set target gene_${build}_${dbname}.tsv
	} else {
		set target extra/gene_${build}_${dbname}.tsv
	}
	job gene_${build}_$dbname -targets {
		$target
		$target.gz.tbi
		$target.gz
	} -vars {dest build db} -code {
		set target [gzroot $target].lz4
		file delete $target
		cg download_genes $target $build $db
	        cg maketabix $target
		cg index $target
	}
}

set target gene_${build}_intGene.tsv
job gene_${build}_intGene -deps {
	gene_${build}_refGene.tsv
	extra/gene_${build}_gencode.tsv
} -targets {
	$target
	$target.gz
	$target.gz.tbi
} -vars {dest build db} -code {
	set target [gzroot $target].lz4
	cg intgene {*}$deps | cg lz4 > $target.temp.lz4
	file rename -force $target.temp.lz4 $target
	cg maketabix $target
	cg lz4index $target
	cg index $target
}

job reg_${build}_genes -deps {
	gene_${build}_refGene.tsv
} -targets {
	extra/reg_${build}_genes.tsv
} -code {
	set target [gzroot $target].lz4
	exec cg cat -fields {chrom start end geneid} {*}$deps | cg select -s {chrom start end geneid} -f {chrom {start=$start - 2000} {end=$end + 2000} geneid} | cg regcollapse | cg lz4 > $target.temp.lz4
	file rename -force $target.temp.lz4 $target
	cg lz4index $target
}

job reg_refcoding -deps {
	gene_${build}_refGene.tsv
} -targets {
	extra/reg_${build}_refcoding.tsv
} -code {
	set target [gzroot $target].lz4
	cg gene2reg $dep | cg select -q {$type eq "CDS"} | cg select -s - | cg regjoin | cg lz4 > $target.temp.lz4
	file rename -force $target.temp.lz4 $target
	cg lz4index $target
}

job reg_exome_refGene -deps {
	gene_${build}_refGene.tsv
} -targets {
	extra/reg_${build}_exome_refGene.tsv
} -code {
	set target [gzroot $target].lz4
	cg gene2reg $dep | cg select -q {$type in "CDS UTR RNA"} | cg select -s - | cg regjoin | cg lz4 > $target.temp.lz4
	file rename -force $target.temp.lz4 $target
	cg lz4index $target
}

job reg_intcoding -deps {
	gene_${build}_intGene.tsv
} -targets {
	extra/reg_${build}_intcoding.tsv
} -code {
	set target [gzroot $target].lz4
	cg gene2reg $dep | cg select -q {$type eq "CDS"} | cg select -s - | cg regjoin | cg lz4 > $target.temp.lz4
	file rename -force $target.temp.lz4 $target
	cg lz4index $target
}

job reg_exome_intGene -deps {
	gene_${build}_intGene.tsv
} -targets {
	extra/reg_${build}_exome_intGene.tsv
} -code {
	set target [gzroot $target].lz4
	cg gene2reg $dep | cg select -q {$type in "CDS UTR RNA"} | cg select -s - | cg regjoin | cg lz4 > $target.temp.lz4
	file rename -force $target.temp.lz4 $target
	cg lz4index $target
}

file_write extra/reg_${build}_pseudoautosomal.tsv {chromosome	begin	end	name
X	169969759	170931299	PAR1
Y	90745845	91644698	PAR1
}

# homopolymer
job reg_${build}_homopolymer -deps {
	genome_${build}.ifas
} -targets {
	reg_${build}_homopolymer.tsv.lz4
	reg_${build}_homopolymer.tsv.gz
	reg_${build}_homopolymer.tsv.gz.tbi
	reg_${build}_homopolymer.tsv.opt
} -vars {dest build db} -code {
	file_write reg_${build}_homopolymer.tsv.opt "fields\t{base size}\n"
	cg extracthomopolymers genome_${build}.ifas | cg lz4 > reg_${build}_homopolymer.tsv.temp.lz4
	file rename -force reg_${build}_homopolymer.tsv.temp.lz4 reg_${build}_homopolymer.tsv.lz4
        cg maketabix reg_${build}_homopolymer.tsv.lz4
}

# mirbase
job mir_${build}_mirbase -targets {
	mir_${build}_mirbase$mirbaserelease.tsv
	mir_${build}_mirbase$mirbaserelease.tsv.info
} -vars {mirbasegenome mirbaserelease mirbasebuild dest build db} -code {
	set target [gzroot $target].lz4
	if {$mirbasebuild ne $build} {error "error: mirbase $mirbaserelease for build $mirbasebuild (making $build)"}
	cg download_mirbase $target $mirbasegenome $mirbaserelease
}

# link local data in dir
foreach file [glob -nocomplain ../${build}-local/*] {
	catch {
		file delete extra/[file tail $file]
		cplinked $file [file tail $file]
	}
}

job extragenome -deps {
	genome_${build}.ifas
	genome_${build}.ifas.index
	genome_${build}.ssa
} -vars build -targets {
	extra/genome_${build}.ifas extra/genome_${build}.ifas.fai extra/genome_${build}.ifas.index
	genome_${build}.fa genome_${build}.fa.fai genome_${build}.fa.index
	extra/genome_${build}.ssa
} -code {
	mklink genome_${build}.ifas extra/genome_${build}.ifas
	mklink genome_${build}.ifas.fai extra/genome_${build}.ifas.fai
	mklink genome_${build}.ifas.index extra/genome_${build}.ifas.index 
	mklink genome_${build}.ifas genome_${build}.fa
	mklink genome_${build}.ifas.fai genome_${build}.fa.fai
	mklink genome_${build}.ifas.index genome_${build}.fa.index 
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
