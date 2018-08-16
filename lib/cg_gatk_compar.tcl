proc gatk_compar_job args {
	upvar job_logdir job_logdir
	set distrreg chr
	set dbdir {}
	set usecombinegvcfs 0
	set batchsize 50
	set maxmem 10
	set vqsr {}
	set newqual true
	set gatkres {}
	set cmdline [list cg gatk_compar {*}$args]
	cg_options gatk_compar args {
		-dbdir {
			set dbdir $value
		}
		-distrreg {
			set distrreg $value
		}
		-maxmem {
			set maxmem $value
		}
		-newqual {
			if {[true $value]} {
				set newqual true
			} else {
				set newqual false
			}
		}
		-batchsize {
			set batchsize $value
		}
		-gatkres {
			set gatkres [file_absolute $value]
		}
		-vqsr {
			set vqsr $value
		}
		-usecombinegvcfs {
			set usecombinegvcfs $value
		}
	} {resultvcf} 2 ... {
		combines gvcfs into joint variant called vcf using gatk tools
	}
	if {$gatkres eq ""} {
		set gatkres $dbdir/gatkres
	}
	if {$vqsr ne ""} {
		# Variant (Quality Score) Recalibration resources and settings
		set indelannotations {QD DP FS SOR MQRankSum ReadPosRankSum InbreedingCoeff}
		set snpannotations {QD DP FS SOR MQ MQRankSum ReadPosRankSum InbreedingCoeff}
		set indelreslist {}
		set iresfound 0
		foreach {pattern res} {
			Mills_and_1000G_gold_standard.indels*.vcf mills,known=false,training=true,truth=true,prior=12.0
			dbsnp*.vcf dbsnp,known=true,training=false,truth=false,prior=2.0
		} {
			set file [jobglob1 $gatkres/$pattern]
			if {[file exists $file]} {
				set iresfound 1
			}
		}
		set sresfound 0
		set snpreslist {}
		foreach {pattern res} {
			hapmap*.sites.vcf hapmap,known=false,training=true,truth=true,prior=15.0
			*omni*.sites.vcf omni,known=false,training=true,truth=true,prior=12.0
			1000G_phase1.snps.high_confidence.*.sites.vcf 1000G,known=false,training=true,truth=false,prior=10.0
			dbsnp*.vcf dbsnp,known=true,training=false,truth=false,prior=2.0
		} {
			set file [jobglob1 $gatkres/$pattern]
			if {[file exists $file]} {
				set sresfound 1
			}
		}
		if {!$iresfound} {
			puts stderr "No INDEL recalibration resources found. (use -gatkres option?)"
			set vqsr ""
		}
		if {!$sresfound} {
			puts stderr "No SNP recalibration resources found. (use -gatkres option?)"
			set vqsr ""
		}
	}
	set resultvcf [file_absolute $resultvcf]
	set resulttsv [file root [gzroot $resultvcf]].tsv.lz4
	set name [file root [file tail [gzroot $resultvcf]]]
	job_logdir [file dir $resultvcf]/log_jobs
	job_logfile [file dir $resultvcf]/gatk_compar_$name [file dir $resultvcf] $cmdline \
		{*}[versions gatk]
	set dbdir [dbdir $dbdir]
	if {![info exists gatkres]} {
		set gatkres $dbdir/gatkres
	}
	set refseq [lindex [glob $dbdir/genome_*.ifas] 0]
	set gatkrefseq [gatk_refseq_job $refseq]
	if {$usecombinegvcfs} {
		set filesarg {}
		set deps [list $gatkrefseq]
		foreach file $args {
			if {[file extension $file] eq ".gz" && ![file exists $file.tbi]} {
				exec tabix -p vcf $file
			}
			lappend filesarg --variant $file
			lappend deps $file $file.tbi
		}
		job gatkbp_combine-$name -mem [expr {$maxmem + 4}]g -cores 1 \
		 -deps $deps \
		 -targets {$resultvcf} \
		 -vars {region gatkrefseq resultvcf samplemapfile maxmem filesarg} \
		 -code {
			set combinedgvcf [file root [gzroot $resultvcf]]-combined.gvcf
			gatkexec [list -XX:ParallelGCThreads=1 -d64 -Xms${maxmem}g -Xmx${maxmem}g] CombineGVCFs \
				-R $gatkrefseq \
				{*}$filesarg \
				-O $combinedgvcf
			#
			# to vcf
			gatkexec [list -XX:ParallelGCThreads=1 -d64 -Xms${maxmem}g -Xmx${maxmem}g] GenotypeGVCFs \
				-R $gatkrefseq \
				-V $combinedgvcf \
				-O $resultvcf.temp[file extension $resultvcf] \
				-G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation
			file delete $combinedgvcf
			file rename $resultvcf.temp[file extension $resultvcf] $resultvcf
		}
	} else {
		set regions [distrreg_regs $distrreg $refseq]
		# sample mapping file
		set samplemapfile $resultvcf.samplemap
		set o [open $samplemapfile w]
		set deps [list $gatkrefseq]
		foreach file $args {
			if {[file extension $file] eq ".gz"} {
				job tabix-[file tail $file] -deps $file -targets $file.tbi -code {
					exec tabix -p vcf $dep
				}
			}
			set analysis [file_analysis $file]
			puts $o $analysis\t[file_absolute $file]
			lappend deps $file $file.tbi
		}
		close $o
		# run
		set vcffiles {}
		foreach region $regions {
			job gatkbp_combine-$name-$region -mem [expr {$maxmem + 4}]g -cores 1 \
			 -skip $resultvcf \
			 -deps $deps \
			 -targets {$resultvcf.$region} \
			 -vars {region gatkrefseq resultvcf samplemapfile batchsize maxmem newqual} \
			 -code {
				# combine using GenomicsDBImport
				set tempfile [tempfile]
				file delete $tempfile
				gatkexec [list -XX:ParallelGCThreads=1 -d64 -Xms${maxmem}g -Xmx${maxmem}g] GenomicsDBImport \
					--genomicsdb-workspace-path $tempfile \
					--sample-name-map $samplemapfile \
					--overwrite-existing-genomicsdb-workspace true \
					--TMP_DIR [scratchdir] \
					--reader-threads 2 \
					--batch-size $batchsize \
					--intervals $region >@ stdout 2>@stderr
				#
				# to vcf
				gatkexec [list -XX:ParallelGCThreads=1 -d64 -Xms${maxmem}g -Xmx${maxmem}g] GenotypeGVCFs \
					-R $gatkrefseq \
					-V gendb://$tempfile \
					-O $resultvcf.$region.temp \
					-G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \
					--verbosity ERROR \
					--TMP_DIR [scratchdir] \
					-new-qual $newqual \
					--only-output-calls-starting-in-intervals \
					--intervals $region >@ stdout 2>@stderr
				file delete -force $tempfile
				file rename $resultvcf.$region.temp $resultvcf.$region
				file delete $resultvcf.$region.temp.idx
			}
			lappend vcffiles $resultvcf.$region
		}
		set deps $vcffiles ; lappend deps $samplemapfile
		job gatkbp_combine-$name -deps $deps -targets {$resultvcf} \
		 -rmtargets $deps \
		 -vars {vcffiles resultvcf} \
		 -code {
			set tempfile $resultvcf.temp[file extension $resultvcf]
			cg_vcfcat -o $tempfile {*}$vcffiles
			file rename $tempfile $resultvcf
			file delete {*}$deps
		}
	}
	gatk_IndexFeatureFile_job $resultvcf
	#
	if {$vqsr ne ""} {
		gatk_vqsr_job -dbdir $dbdir -gatkres $gatkres $resultvcf $vqsr
	}
	#
	# convert to tsv
	set mincoverage 8
	set mingenoqual 25
	set split 1
	job gatkbp_combine-$name -deps {
		$resultvcf
	} -targets {
		$resulttsv
	} -vars {
		resultvcf resulttsv split mincoverage mingenoqual
	} -code {
		#set fields {chromosome begin end type ref alt quality alleleSeq1 alleleSeq2}
		#lappend fields [subst {sequenced=if(\$genoqual < $mingenoqual || \$coverage < $mincoverage,"u","v")}]
		#lappend fields [subst {zyg=if(\$genoqual < $mingenoqual || \$coverage < $mincoverage,"u",\$zyg)}]
		#lappend fields *
		#exec cg vcf2tsv -split $split -removefields {
		#	name filter AN AC AF AA ExcessHet InbreedingCoeff MLEAC MLEAF NDA RPA RU STR
		#} tmp/combined.vcf | cg select -f $fields > tmp/combined.tsv
		#
		exec cg vcf2tsv -split $split -removefields {
			name filter AN AC AF AA ExcessHet InbreedingCoeff MLEAC MLEAF NDA RPA RU STR
		} $resultvcf $resulttsv
	}
	
}

proc cg_gatk_compar {args} {
	set args [job_init {*}$args]
	gatk_compar_job {*}$args
	job_wait
}

