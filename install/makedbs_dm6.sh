#!/bin/sh
# the next line restarts using tclsh \
exec cg source "$0" "$@"

set build dm6
set dbs_reg_collapse {cytoBand rmsk simpleRepeat microsat cons27way}
set dbs_reg_join {}
set dbs_other {refLink}
set dbs_tabix {rmsk simpleRepeat}
set dbs_var {}
set dbs_gene {refGene extra/ensGene extra/genscan}
set mirbasegenome dme
set mirbaserelease 21

logverbose 2

if {![info exists argv]} {set argv {}}
set argv [job_init {*}$argv]
foreach {dest webcache} $argv break
if {![info exists dest]} {set dest /complgen/refseqnew}
if {[info exists webcache]} {set env(webcache) $webcache}

putslog "Creating dir $dest/$build"

# download genome
# ===============
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
foreach db $dbs_reg_collapse
	cytoBandIdeo rmsk simpleRepeat 
} {
	job reg_${build}_$db -targets {reg_${build}_${db}.tsv} -vars {dest build db} -code {
		cg download_ucsc $target.ucsc ${build} $db
		cg regcollapse $target.ucsc > $target.temp
		file rename -force $target.ucsc.info $target.info
		file rename -force $target.temp $target
		file delete $target.ucsc
	}
}

# join regions
foreach db $dbs_reg_join {
	job reg_${build}_$db -targets {reg_${build}_${db}.tsv} -vars {dest build db} -code {
		cg download_ucsc $target.ucsc ${build} $db
		cg regjoin $target.ucsc > $target.temp
		file delete $target.ucsc
		file rename -force $target.ucsc.info $target.info
		file rename -force $target.temp $target
	}
}

# other databases
foreach db $dbs_other {
	job other_${build}_$db -vars {dest build db} -targets {other_${build}_${db}.tsv} -code {
		cg download_ucsc other_${build}_${db}.tsv ${build} $db
	}
}

foreach db $dbs_tabix {
	job maketabix_${build}_$db -deps {reg_${build}_${db}.tsv} -targets {reg_${build}_${db}.tsv.gz.tbi reg_${build}_${db}.tsv.gz} -vars {build db} -code {
		cg maketabix reg_${build}_${db}.tsv
	}
}

# var dbs
foreach db $dbs_var {
	job $db -targets {var_${build}_${db}.tsv} -vars {dest build db} -code {
		cg download_ucsc $target.temp $build $db
		cg select -f {chrom start end type ref alt name freq} $target.temp $target.temp2
		file rename -force $target.temp2 $target
		file delete $target.tmp
	}
	job maketabix_${build}_$db -deps {var_${build}_${db}.tsv} -targets {var_${build}_${db}.tsv.gz.tbi var_${build}_${db}.tsv.gz} -vars {dest build db} -code {
		cg maketabix var_${build}_${db}.tsv
	}
}

# genes
foreach db $dbs_gene {
	set dbname [file tail $db]
	set dir [file dir $db]
	set target $dir/gene_${build}_${dbname}.tsv
	job gene_${build}_$dbname -targets {$target $target.gz.tbi $target.gz} -vars {dest build dbname} -code {
		cg download_genes $target $build $dbname
	        cg maketabix $target
		cg index $target
	}
}

# homopolymer
job reg_${build}_homopolymer -deps {genome_${build}.ifas} -targets {reg_${build}_homopolymer.tsv reg_${build}_homopolymer.tsv.gz reg_${build}_homopolymer.tsv.gz.tbi reg_${build}_homopolymer.tsv.opt} -vars {dest build db} -code {
	cg extracthomopolymers genome_${build}.ifas > reg_${build}_homopolymer.tsv.temp
	file rename -force reg_${build}_homopolymer.tsv.temp reg_${build}_homopolymer.tsv
        cg maketabix reg_${build}_homopolymer.tsv
	file_write reg_${build}_homopolymer.tsv.opt "fields\t{base size}\n"
}

# mirbase
job mir_${build}_mirbase -targets {$dest/${build}/mir_${build}_mirbase$mirbaserelease.tsv $dest/${build}/mir_${build}_mirbase.info} -vars {mirbasegenome mirbaserelease dest build db} -code {
	set resultfile $dest/${build}/mir_${build}_mirbase$mirbaserelease.tsv
	cg downloadmirbase $resultfile.temp $mirbasegenome $mirbaserelease
	file rename -force $resultfile.temp $resultfile
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
