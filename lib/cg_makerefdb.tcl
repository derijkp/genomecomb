proc makerefdb_job {args} {
	upvar job_logdir job_logdir
	set cmdline [clean_cmdline cg makerefdb {*}$args]

	# default settings
	# ----------------

	set webcache {}
	set genomeurl {}
	set regionsdb_collapse {
		cytoBand microsat oreganno rmsk simpleRepeat 
		tRNAs wgRna
		phastConsElements100way phastConsElements30way
	}
	set regionsdb_join {
		chainSelf dgvMerged genomicSuperDups
	}
	set genesdb [list {refGene int} \
		{ensGene int} \
		{knownGene extra int} \
		{genscan extra} \
		{augustusGene extra} \
		{lincRNAsTranscripts lincRNA} \
	]
	set dbsnp {}
	set refSeqFuncElemsurl {}
	set mirbase {}
	set pseudoautosomal {}
	set transcriptsurl {}
	set transcriptsgtf {}

	# options
	# -------
	cg_options makerefdb args {
		-genomeurl {
			set genomeurl $value
		}
		-pseudoautosomal {
			set pseudoautosomal $value
		}
		-genesdb {
			set genesdb $value
		}
		-dbsnp {
			set dbsnp $value
		}
		-refSeqFuncElemsurl {
			set refSeqFuncElemsurl $value
		}
		-mirbase {
			set mirbase $value
		}
		-regionsdb_collapse {
			set regionsdb_collapse $value
		}
		-regionsdb_join {
			set regionsdb_join $value
		}
		-transcriptsurl {
			set transcriptsurl $value
		}
		-transcriptsgtf {
			set transcriptsgtf $value
		}
		-webcache {
			set webcache [file_absolute $value]
		}
	} {dbdir} 1 1
	set dbdir [file_absolute $dbdir]
	set build [file tail $dbdir]
	set dest [file dir $dbdir]
	if {$webcache eq ""} {
		set webcache $dest/webcache
	}
	set env(webcache) $webcache

	# prepare
	# -------
	putslog "Installing in $dest/$build"
	file mkdir ${dest}/${build}
	cd ${dest}/${build}
	file mkdir extra

	# logfile
	# -------
	job_logfile $dbdir/log_makerefdb_$build $dbdir $cmdline {*}[versions samtools]

	# download
	# ========
	#

	# readme
	set c [file_read $::appdir/help/dbdir_README.txt]
	regsub {version: [0-9.]+} $c "version: [lindex [timestamp] 0]\ntime: [lindex [timestamp] 0]" c
	file_write README_dbdir.txt $c

	set refseq genome_${build}.ifas

	# download genome
	job genome_${build} -targets {
		genome_${build}.ifas
		genome_${build}.ifas.fai
		extra/reg_${build}_fullgenome.tsv
		genome_${build}.ifas.index
	} -vars {
		build genomeurl
	} -code {
		cg download_genome -alt 0 -url $genomeurl genome_${build}.ifas ${build} 2>@ stderr
		file rename -force -- reg_genome_${build}.tsv extra/reg_${build}_fullgenome.tsv
		cg zst -i 1 extra/reg_${build}_fullgenome.tsv
	}

	job genome_${build}_cindex -deps {
		genome_${build}.ifas
	} -targets {
		genome_${build}.ssa
	} -code {
		cg make_genomecindex $dep
	}

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

	# make bwa version of genome
	refseq_bwa_job genome_${build}.ifas

	# make ngmlr version of genome
	refseq_ngmlr_job genome_${build}.ifas ont

	# make minimap2 versions of genome
	refseq_minimap2_job genome_${build}.ifas sr
	refseq_minimap2_job genome_${build}.ifas map-ont
	refseq_minimap2_job genome_${build}.ifas splice

	job extragenome -deps {
		genome_${build}.ifas
		genome_${build}.ifas.index
		genome_${build}.ssa
	} -vars build -targets {
		extra/genome_${build}.ifas extra/genome_${build}.ifas.fai extra/genome_${build}.ifas.index
		genome_${build}.fa genome_${build}.fa.fai genome_${build}.fa.index
		extra/genome_${build}.ssa
	} -code {
		mklink -matchtime 0 genome_${build}.ifas extra/genome_${build}.ifas
		mklink -matchtime 0 genome_${build}.ifas.fai extra/genome_${build}.ifas.fai
		mklink -matchtime 0 genome_${build}.ifas.index extra/genome_${build}.ifas.index 
		mklink -matchtime 0 genome_${build}.ifas genome_${build}.fa
		mklink -matchtime 0 genome_${build}.ifas.fai genome_${build}.fa.fai
		mklink -matchtime 0 genome_${build}.ifas.index genome_${build}.fa.index 
		mklink -matchtime 0 genome_${build}.ssa extra/genome_${build}.ssa
	}

	job genome_${build}_dict -deps {
		genome_${build}.fa
	} -targets {
		genome_${build}.dict
	} -code {
		exec samtools dict -o $target.temp $dep
		file rename $target.temp $target
	}

	# genome in extra
	foreach file [jobglob genome_*] {
		set target extra/[file tail $file]
		file delete $target
		mklink $file $target
	}

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

	# region databases (ucsc)
	# -----------------------

	# you can explicitely download info on a database using:
	# cg download_ucscinfo resultfile ${build} dbname

	# collapse regions
	foreach db $regionsdb_collapse {
		if {![file exists [gzfile reg_${build}_${db}.tsv]] && ![cg_check_ucsc $build $db]} {
			putslog "skipping $db, track not found for $build"
			continue
		}
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

	foreach db {
		rmsk simpleRepeat
	} {
		job maketabix_${build}_$db -optional 1 -deps {
			reg_${build}_${db}.tsv
		} -targets {
			reg_${build}_${db}.tsv.gz.tbi
			reg_${build}_${db}.tsv.gz
		} -vars {build db} -code {
			cg maketabix $dep
		}
	}

	# join regions
	foreach db $regionsdb_join {
		if {![file exists [gzfile reg_${build}_${db}.tsv]] && ![cg_check_ucsc $build $db]} {
			putslog "skipping $db, track not found for $build"
			continue
		}
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

	# refSeqFuncElemsurl
	if {$refSeqFuncElemsurl ne ""} {
		job reg_${build}_refSeqFuncElems -targets {
			reg_${build}_refSeqFuncElems.tsv.zst
		} -vars {
			dest build refSeqFuncElemsurl
		} -code {
			set target [gzroot $target].zst
			file mkdir $target.temp
			set tail [file tail $refSeqFuncElemsurl]
			wgetfile $refSeqFuncElemsurl $target.temp/$tail
			cg gff2tsv $target.temp/$tail $target.temp/reg.tsv
			set temp [string trim [cg select -q {$type eq "region" and $chromosome regexp "^NC_0" and $attr_chromosome ne ""} -g {chromosome * attr_chromosome *} $target.temp/reg.tsv]]
			set temp [lrange [split $temp \n] 1 end]
			set chrcode {chromosome=if(}
			list_foreach {chrid chr} $temp {
				append chrcode "\$chromosome eq \"$chrid\",\"$chr\","
			}
			append chrcode "\"?\")"
			set fields [list $chrcode {*}[list_common {
				begin end type strand source phase ID Dbxref gbkey Note experiment function regulatory_class standard_name
			} [cg select -h $target.temp/reg.tsv]]]
		#	set fields {chromosome begin end type strand source phase ID Dbxref gbkey Note experiment function regulatory_class standard_name}
			cg select -stack 1 -v 2 -overwrite 1 -q {$source eq "RefSeqFE" and $type ne "biological_region" and $chromosome regexp "^NC_0"} -f $fields $target.temp/reg.tsv $target.temp/refseqfe.tsv
	
			cg select -s - -overwrite 1 $target.temp/refseqfe.tsv $target.temp/srefseqfe.tsv
			cg regcollapse $target.temp/srefseqfe.tsv | cg zst -compressionlevel 11 > $target.temp/result.tsv.zst
			# opt and info
			file_write [gzroot $target].opt "fields\t{type}\n"
			set version [timestamp]
			set temp [exec head -20 $target.temp/srefseqfe.tsv]
			regexp {Annotation Release ([0-9.]+)} $temp temp version
			file_write [gzroot $target].info [subst [deindent {
				= RefSeqFE (Refseq functional elements) =
				
				== Download info ==
				dbname	RefSeqFE
				version	$version
				source	$refSeqFuncElemsurl
				time	[timestamp]
				
				== Description ==
				
				More info on https://www.ncbi.nlm.nih.gov/refseq/functionalelements/
				
				== Category ==
				Annotation
			}]]
			file rename -force $target.temp/result.tsv.zst $target
			file delete -force $target.temp
		}
	}

	# genes
	# -----
	set intdeps {}
	set regdeps {}
	foreach line $genesdb {
		set db [lindex $line 0]
		set dbname $db
		set extra 0; set int 0 ; set reg 0
		foreach el [lrange $line 1 end] {
			if {$el eq "extra"} {
				set extra 1
			} elseif {$el eq "int"} {
				set int 1
			} elseif {$el eq "reg"} {
				set reg 1
			} else {
				set dbname $el
			}
		}
		if {$extra} {
			set target extra/gene_${build}_${dbname}.tsv
		} else {
			set target gene_${build}_${dbname}.tsv
		}
		if {![file exists [gzfile $target]] && ![cg_check_ucsc $build $db]} {
			putslog "skipping $db, track not found for $build"
			continue
		}
		if {$int} {
			lappend intdeps $target.zst
		}
		if {$reg} {
			lappend regdeps $target.zst
		}
		job gene_${build}_$dbname -targets {
			$target.zst
			$target.gz.tbi
			$target.gz
		} -vars {dest build db} -code {
			set target [gzroot $target].zst
			file delete $target
			cg download_genes $target $build $db
		        cg maketabix $target
			cg index $target
		}
	}

	set target gene_${build}_intGene.tsv
	job gene_${build}_intGene -deps $intdeps -targets {
		$target.zst
		$target.gz
		$target.gz.tbi
	} -vars {dest build db} -code {
		set target [gzroot $target].zst
		cg intgene {*}$deps {*}[compresspipe $target 12] > $target.temp
		file rename -force -- $target.temp $target
		cg maketabix $target
		cg zindex $target
		cg index $target
	}

	job reg_${build}_genes -deps $regdeps -targets {
		extra/reg_${build}_genes.tsv
	} -code {
		set target [gzroot $target].zst
		exec cg cat -fields {chrom start end geneid} {*}$deps | cg select -s {chrom start end geneid} -f {chrom {start=$start - 2000} {end=$end + 2000} geneid} | cg regcollapse {*}[compresspipe $target 12] > $target.temp
		file rename -force -- $target.temp $target
		cg zindex $target
	}

	set dep gene_${build}_refGene.tsv
	if {![jobfileexists $dep]} {
		set dep gene_${build}_ncbiRefSeq.tsv
	}
	if {![jobfileexists $dep]} {
		set dep extra/gene_${build}_ncbiRefSeq.tsv
	}
	job reg_refcoding -deps {
		$dep
	} -targets {
		extra/reg_${build}_refcoding.tsv
	} -vars {
		build
	} -code {
		mkdbs_write_info $target Regions {
			Coding regions extracted from refGene
		} source $dep build $build
		set target [gzroot $target].zst
		cg gene2reg $dep | cg select -q {$type eq "CDS"} | cg select -s - | cg regjoin {*}[compresspipe $target 12] > $target.temp
		file rename -force -- $target.temp $target
		cg zindex $target
	}

	job reg_exome_refGene -deps {
		$dep
	} -targets {
		extra/reg_${build}_exome_refGene.tsv
	} -vars {
		build
	} -code {
		mkdbs_write_info $target Regions {
			Exome regions (CDS, UTR and RNA) extracted from refGene
		} source $dep build $build
		set target [gzroot $target].zst
		cg gene2reg $dep | cg select -q {$type in "CDS UTR RNA"} | cg select -s - | cg regjoin {*}[compresspipe $target 12] > $target.temp
		file rename -force -- $target.temp $target
		cg zindex $target
	}

	job reg_intcoding -deps {
		gene_${build}_intGene.tsv
	} -targets {
		extra/reg_${build}_intcoding.tsv
	} -vars {
		build
	} -code {
		mkdbs_write_info $target Regions {
			Coding regions extracted from intGene
		} source $dep build $build
		set target [gzroot $target].zst
		cg gene2reg $dep | cg select -q {$type eq "CDS"} | cg select -s - | cg regjoin {*}[compresspipe $target 12] > $target.temp
		file rename -force -- $target.temp $target
		cg zindex $target
	}

	job reg_exome_intGene -deps {
		gene_${build}_intGene.tsv
	} -targets {
		extra/reg_${build}_exome_intGene.tsv
	} -vars {
		build
	} -code {
		mkdbs_write_info $target Regions {
			Exome regions (CDS, UTR and RNA) extracted from intGene
		} source $dep build $build
		set target [gzroot $target].zst
		cg gene2reg $dep | cg select -q {$type in "CDS UTR RNA"} | cg select -s - | cg regjoin {*}[compresspipe $target 12] > $target.temp
		file rename -force -- $target.temp $target
		cg zindex $target
	}

	if {$pseudoautosomal ne ""} {
		file_write extra/reg_${build}_pseudoautosomal.tsv $pseudoautosomal
	}

	# mirbase
	if {$mirbase ne ""} {
		foreach {mirbasegenome mirbaserelease mirbasebuild} [split $mirbase -:] break
		job mir_${build}_mirbase -targets {
			mir_${build}_mirbase$mirbaserelease.tsv
			mir_${build}_mirbase$mirbaserelease.tsv.info
		} -vars {mirbasegenome mirbaserelease mirbasebuild dest build db} -code {
			set target [gzroot $target].zst
			if {$mirbasebuild ne $build} {error "error: mirbase $mirbaserelease for build $mirbasebuild (making $build)"}
			cg download_mirbase $target $mirbasegenome $mirbaserelease
		}
	}
	
	# dbsnp
	if {$dbsnp ne ""} {
		job dbsnp -targets {
			var_${build}_dbsnp.tsv.zst
			var_${build}_dbsnp.tsv.opt
			var_${build}_dbsnp.tsv.gz
			var_${build}_dbsnp.tsv.gz.tbi
		} -vars {dest build dbsnp} -code {
			set target [gzroot $target].zst
			file_write [gzroot $target].opt "fields\t{name}\n"
			cg download_dbsnp $target ${build} snp$dbsnp 2>@ stderr
			cg zindex $target
			cg maketabix $target
		}
	
		job dbsnpCommon -optional 1 -targets {
			var_${build}_dbsnpCommon.tsv.zst
			var_${build}_dbsnpCommon.tsv.gz
			var_${build}_dbsnpCommon.tsv.gz.tbi
		} -vars {dest build dbsnp} -code {
			set target [gzroot $target].zst
			file_write [gzroot $target].opt "fields\t{freqp}\n"
			cg download_dbsnp $target ${build} snp${dbsnp}Common 2>@ stderr
			cg zindex $target
			cg maketabix $target
		}
	}

	set fullgenome [list [jobgzfile $dbdir/extra/reg_*_fullgenome.tsv]]
	set rmskfile [jobgzfile $dbdir/reg_*_rmsk.tsv]
	job norep100000file -optional 1 -deps {
		$fullgenome
		$rmskfile
	} -targets {
		$dbdir/extra/reg_${build}_norep100000.tsv.zst
	} -vars {
		dbdir
	} -code {
		distrreg_norep100000file $dbdir
	}

	set genefile [jobgzfile \
		$dbdir/gene_*_intGene.tsv.zst \
		$dbdir/extra/gene_*_gencode.tsv.zst \
		$dbdir/extra/gene_*.tsv.zst \
		$dbdir/gene_*_refGene.tsv.zst \
	]
	job nolowgene -optional 1 -deps {
		$fullgenome
		$genefile
	} -targets {
		$dbdir/extra/reg_${build}_nolowgene200k.tsv.zst
	} -vars {
		dbdir
	} -code {
		distrreg_nolowgene $dbdir
	}

	if {$transcriptsgtf ne ""} {
		if {$transcriptsurl eq ""} {error "-transcriptsurl empty"}
		job gtf -targets {
			$transcriptsgtf
		} -vars {transcriptsurl transcriptsgtf} -code {
			set ext [file ext $transcriptsurl]
			wgetfile $transcriptsurl $transcriptsgtf.temp$ext
			if {[gziscompressed $transcriptsgtf.temp$ext] && ![gziscompressed $transcriptsgtf]} {
				cg unzip $transcriptsgtf.temp$ext
				file rename $transcriptsgtf.temp $transcriptsgtf
			} else {
				file rename $transcriptsgtf.temp$ext $transcriptsgtf
			}
		}
	
		# make star versions of the genome
		refseq_star_job $refseq $transcriptsgtf 2
	}
}

proc cg_makerefdb {args} {
	set args [job_init {*}$args]
#	putslog pwd:\ [pwd]
#	putslog cmd:\ [list cg makerefdb {*}$args]
#	putslog jobargs:\ [job_curargs]
#	if {[llength $args] < 1} {
#		errorformat makerefdb
#	}
	makerefdb_job {*}$args
	job_wait
}

