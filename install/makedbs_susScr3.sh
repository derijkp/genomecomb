#!/bin/sh
# the next line restarts using tclsh \
exec cg source "$0" "$@"

set build susScr3
if {![info exists argv]} {set argv {}}
set argv [job_init {*}$argv]
if {[llength $argv]} {
	set dest [lindex $argv 0]
} else {
	set dest /complgen/refseq/
}


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
	file rename -force reg_genome_${build}.tsv extra/reg_${build}_fullgenome.tsv
}

set ifasfile genome_${build}.ifas
job genome_${build}_cindex -deps {$ifasfile} -targets {genome_${build}.ssa} -code {
	cg make_genomecindex $dep
}

job reg_${build}_sequencedgenome -vars {dest build} -deps {genome_${build}.ifas} -targets {extra/reg_${build}_sequencedgenome.tsv} -code {
	cg calcsequencedgenome $dest ${build}
}

# region databases (ucsc)
# you can explicitely download info on the databases using:
# cg downloaddbinfo ${dest} ${build} simpleRepeat microsat rmsk genomicSuperDups chainSelf

# collapse regions
foreach db {
	cytoBandIdeo rmsk simpleRepeat 
} {
	job reg_${build}_$db -targets {reg_${build}_${db}.tsv} -vars {dest build db} -code {
		cg downloaddb $dest ${build} $db
		cg collapseoverlap ucsc_${build}_${db}.tsv
		file delete ucsc_${build}_${db}.tsv
	}
}

# join regions
foreach db {
} {
	job reg_${build}_$db -targets {reg_${build}_${db}.tsv} -vars {dest build db} -code {
		cg downloaddb $dest ${build} $db
		cg regjoin ucsc_${build}_${db}.tsv > reg_${build}_${db}.tsv.temp
		file delete ucsc_${build}_${db}.tsv
		file rename -force reg_${build}_${db}.tsv.temp reg_${build}_${db}.tsv
	}
}

# other databases
foreach db {refLink} {
	job other_${build}_$db -vars {dest build db} -targets {other_${build}_${db}.tsv} -code {
		cg downloaddb $dest ${build} $db
		file rename -force ucsc_${build}_${db}.tsv other_${build}_${db}.tsv
	}
}

foreach db {
	rmsk simpleRepeat
} {
	job maketabix_${build}_$db -deps {reg_${build}_${db}.tsv} -targets {reg_${build}_${db}.tsv.gz.tbi reg_${build}_${db}.tsv.gz} -vars {build db} -code {
		cg maketabix reg_${build}_${db}.tsv
		exec gunzip -c reg_${build}_${db}.tsv.gz > reg_${build}_${db}.tsv
	}
}

# dbsnp
job dbsnp138 -targets {${dest}/$build/var_${build}_snp138.tsv} -vars {dest build} -code {
	cd $dest/$build
	cg downloaddb ${dest} $build snp138
}

foreach db {
	snp138
} {
	job maketabix_${build}_$db -deps {var_${build}_${db}.tsv} -targets {var_${build}_${db}.tsv.gz.tbi var_${build}_${db}.tsv.gz} -vars {dest build db} -code {
		cg maketabix var_${build}_${db}.tsv
		# exec gunzip -c var_${build}_${db}.tsv.gz > var_${build}_${db}.tsv
		cg select -f {chrom start end type ref alt name freq} ${dest}/${build}/var_${build}_${db}.tsv.gz ${dest}/${build}/var_${build}_${db}.tsv
	}
}

# genes
foreach db {
	refGene ensGene genscan
} {
	job gene_${build}_$db -targets {gene_${build}_${db}.tsv gene_${build}_${db}.tsv.gz.tbi gene_${build}_${db}.tsv.gz} -vars {dest build db} -code {
	        cg downloaddb ${dest} ${build} $db
	        cg select -s - ucsc_${build}_${db}.tsv gene_${build}_${db}.tsv.temp
		file rename -force gene_${build}_${db}.tsv.temp gene_${build}_${db}.tsv
	        file rename -force reg_${build}_${db}.info gene_${build}_${db}.info
	        cg maketabix gene_${build}_${db}.tsv
	        exec gunzip -c gene_${build}_${db}.tsv.gz > gene_${build}_${db}.tsv.temp
		file rename -force gene_${build}_${db}.tsv.temp gene_${build}_${db}.tsv
		file delete ucsc_${build}_${db}.tsv
	}
}

# homopolymer
job reg_${build}_homopolymer -deps {genome_${build}.ifas} -targets {reg_${build}_homopolymer.tsv reg_${build}_homopolymer.tsv.gz reg_${build}_homopolymer.tsv.gz.tbi reg_${build}_homopolymer.tsv.opt} -vars {dest build db} -code {
	cg extracthomopolymers genome_${build}.ifas > reg_${build}_homopolymer.tsv.temp
	file rename -force reg_${build}_homopolymer.tsv.temp reg_${build}_homopolymer.tsv
        cg maketabix reg_${build}_homopolymer.tsv
        gunzip reg_${build}_homopolymer.tsv.gz reg_${build}_homopolymer.tsv
	file_write reg_${build}_homopolymer.tsv.opt "fields\t{base size}\n"
}

# mirbase
job reg_${build}_mirbase -targets {$dest/${build}/reg_${build}_mirbase.tsv $dest/${build}/reg_${build}_mirbase.tsv.opt $dest/${build}/reg_${build}_mirbase.info} -vars {dest build db} -code {
	set organism ssc
	cd $dest/${build}
	file_write $dest/${build}/reg_${build}_mirbase.tsv.opt "fields\t{ID}\n"
	exec -ignorestderr wget -c --tries=45 --directory-prefix=${dest}/tmp/${build} ftp://mirbase.org/pub/mirbase/20/genomes/$organism.gff2
	cg gff2sft ${dest}/tmp/${build}/$organism.gff2 ${dest}/tmp/${build}/reg_${build}_mirbase.tsv.temp
	cg select -s - ${dest}/tmp/${build}/reg_${build}_mirbase.tsv.temp ${dest}/tmp/${build}/reg_${build}_mirbase.tsv.temp2
	file rename -force ${dest}/tmp/${build}/reg_${build}_mirbase.tsv.temp2 reg_${build}_mirbase.tsv
	exec -ignorestderr wget -c ftp://mirbase.org/pub/mirbase/20/README
	file rename -force README reg_${build}_mirbase.info
}

job extragenome -deps {genome_${build}.ifas genome_${build}.ifas.index genome_${build}.ssa} -vars build \
-targets {extra/genome_${build}.ifas extra/genome_${build}.ifas.index extra/genome_${build}.ssa} -code {
	mklink genome_${build}.ifas extra/genome_${build}.ifas
	mklink genome_${build}.ifas.index extra/genome_${build}.ifas.index 
	mklink genome_${build}.ssa extra/genome_${build}.ssa
}
# genome in extra
catch {
	foreach file [glob genome_*] {
		file delete extra/[file tail $file]
		mklink $file extra/[file tail $file]
	}
}

job_wait
