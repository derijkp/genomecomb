#!/bin/sh
# the next line restarts using tclsh \
exec cg source "$0" "$@"

set build susScr3
set mirbasegenome ssc
set mirbaserelease 21
if {![info exists argv]} {set argv {}}

logverbose 2

if {![info exists argv]} {set argv {}}
set argv [job_init {*}$argv]
foreach {dest webcache} $argv break
if {![info exists dest]} {set dest /complgen/refseqnew}
if {![info exists webcache]} {set webcache $dest/webcache}
if {[info exists webcache]} {set env(webcache) $webcache}
set dest [file_absolute $dest]

putslog "Installing in $dest/$build"

# download susScr3
# ================
#
file mkdir ${dest}/${build}
cd ${dest}/${build}
file mkdir extra

job_logdir log_jobs

# download genome
job genome_${build} -vars build -targets {genome_${build}.ifas extra/reg_${build}_fullgenome.tsv} -code {
	cg downloadgenome ${build} genome_${build}.ifas
	file rename -force -- reg_genome_${build}.tsv extra/reg_${build}_fullgenome.tsv
	cg zst -i 1 extra/reg_${build}_fullgenome.tsv
}

set ifasfile genome_${build}.ifas
job genome_${build}_cindex -deps {$ifasfile} -targets {genome_${build}.ssa} -code {
	cg make_genomecindex $dep
}

job reg_${build}_sequencedgenome -vars {dest build} -deps {genome_${build}.ifas} -targets {extra/reg_${build}_sequencedgenome.tsv} -code {
	exec cg calcsequencedgenome --stack 1 $dep {*}[compresspipe $target 12] > $target.temp
	file rename -force -- $target.temp $target
}

# make bwa version of genome
bwarefseq_job genome_${build}.ifas

# region databases (ucsc)
# you can explicitely download info on a database using:
# cg download_ucscinfo resultfile ${build} dbname

# collapse regions
foreach db {
	cytoBandIdeo rmsk simpleRepeat 
} {
	job reg_${build}_$db -targets {
		reg_${build}_${db}.tsv
	} -vars {dest build db} -code {
		set target [gzroot $target].zst
		cg download_ucsc $target.ucsc ${build} $db
		cg regcollapse $target.ucsc > $target.temp
		file delete $target.ucsc
		file rename -force -- $target.ucsc.info [gzroot $target].info
		compress $target.temp $target
	}
}

# join regions
foreach db {
} {
	job reg_${build}_$db -targets {
		reg_${build}_${db}.tsv
	} -vars {dest build db} -code {
		set target [gzroot $target].zst
		cg download_ucsc $target.ucsc ${build} $db
		cg regjoin $target.ucsc > $target.temp
		file delete $target.ucsc
		file rename -force -- $target.ucsc.info [gzroot $target].info
		compress $target.temp $target
	}
}

# other databases
foreach db {refLink} {
	job other_${build}_$db -vars {dest build db} -targets {other_${build}_${db}.tsv} -code {
		cg download_ucsc $target ${build} $db
	}
}

foreach db {
	rmsk simpleRepeat
} {
	job maketabix_${build}_$db -deps {reg_${build}_${db}.tsv} -targets {reg_${build}_${db}.tsv.gz.tbi reg_${build}_${db}.tsv.gz} -vars {build db} -code {
		cg maketabix reg_${build}_${db}.tsv
	}
}

# dbsnp
job dbsnp138 -targets {var_${build}_snp138.tsv} -vars {dest build} -code {
	cg download_dbsnp $target ${build} snp138 2>@ stderr
}

foreach db {
	snp138
} {
	job maketabix_${build}_$db -deps {var_${build}_${db}.tsv} -targets {var_${build}_${db}.tsv.gz.tbi var_${build}_${db}.tsv.gz} -vars {dest build db} -code {
		cg maketabix var_${build}_${db}.tsv
	}
}

# genes
foreach db {
	refGene ensGene genscan
} {
	set dbname $db
	if {$db eq "refGene"} {
		set target gene_${build}_${dbname}.tsv
	} else {
		set target extra/gene_${build}_${dbname}.tsv
	}
	job gene_${build}_$db -targets {$target $target.gz.tbi $target.gz} -vars {dest build db} -code {
		cg download_genes $target $build $db
	        cg maketabix $target
		cg index $target
	}
}

set target gene_${build}_intGene.tsv
job gene_${build}_intGene \
-deps {gene_${build}_refGene.tsv extra/gene_${build}_gencode.tsv extra/gene_${build}_ensGene.tsv} \
-targets {$target $target.gz $target.gz.tbi} -vars {dest build db} -code {
	cg intgene {*}$deps > $target.temp
	file rename -force -- $target.temp $target
	cg maketabix $target
	cg index $target
}

# homopolymer
job reg_${build}_homopolymer -deps {genome_${build}.ifas} -targets {reg_${build}_homopolymer.tsv reg_${build}_homopolymer.tsv.gz reg_${build}_homopolymer.tsv.gz.tbi reg_${build}_homopolymer.tsv.opt} -vars {dest build db} -code {
	cg extracthomopolymers genome_${build}.ifas > reg_${build}_homopolymer.tsv.temp
	file rename -force -- reg_${build}_homopolymer.tsv.temp reg_${build}_homopolymer.tsv
        cg maketabix reg_${build}_homopolymer.tsv
	file_write reg_${build}_homopolymer.tsv.opt "fields\t{base size}\n"
}

# mirbase
job reg_${build}_mirbase -targets {mir_${build}_mirbase$mirbaserelease.tsv mir_${build}_mirbase.info} -vars {mirbasegenome mirbaserelease dest build db} -code {
	set organism $mirbasegenome
	cg downloadmirbase mir_${build}_mirbase$mirbaserelease.tsv $organism $release
}

job extragenome -deps {genome_${build}.ifas genome_${build}.ifas.index genome_${build}.ssa} -vars build \
-targets {
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
catch {
	foreach file [glob genome_*] {
		file delete extra/[file tail $file]
		mklink $file extra/[file tail $file]
	}
}

# compress
foreach file [jobglob *.tsv] {
	job zst_${build}_[file tail $file] -deps {$file} -targets {$file.zst} -vars {dest build} -code {
		cg zst -c 12 -i 1 $dep
	}
}

job_wait
