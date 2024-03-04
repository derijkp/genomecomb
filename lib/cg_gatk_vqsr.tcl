proc gatk_vqsr_job args {
	upvar job_logdir job_logdir
	set dbdir {}
	set usecombinegvcfs 0
	set batchsize 50
	set maxmem 10
	set newqual true
	set cmdline [clean_cmdline cg gatk_vqsr {*}$args]
	cg_options gatkbp_combine args {
		-dbdir {
			set dbdir $value
		}
		-gatkres {
			set gatkres [file_absolute $value]
		}
		-maxmem {
			set maxmem $value
		}
	} {rawvcf resultvcf} 2 2 {
		runs variant quality score recalibration on vcf
	}
	set dbdir [dbdir $dbdir]
	if {![info exists gatkres]} {
		set gatkres $dbdir/gatkres
	}
	set refseq [lindex [glob $dbdir/genome_*.ifas] 0]
	set gatkrefseq [gatk_refseq_job $refseq]
	set rawvcf [file_absolute $rawvcf]
	set resultvcf [file_absolute $resultvcf]
	if {![info exists job_logdir]} {
		set_job_logdir [file dir $resultvcf]/log_jobs
	}
	job_logfile [file dir $resultvcf]/gatk_vqsr_[file tail $resultvcf] [file dir $resultvcf] $cmdline \
		{*}[versions gatk]
	#
	# Variant (Quality Score) Recalibration resources and settings
	set indelannotations {QD DP FS SOR MQRankSum ReadPosRankSum InbreedingCoeff}
	set snpannotations {QD DP FS SOR MQ MQRankSum ReadPosRankSum InbreedingCoeff}
	set indelreslist {}
	foreach {pattern res} {
		Mills_and_1000G_gold_standard.indels*.vcf mills,known=false,training=true,truth=true,prior=12.0
		*dbsnp*.vcf dbsnp,known=true,training=false,truth=false,prior=2.0
	} {
		set file [jobglob1 -checkcompressed 1 $gatkres/$pattern]
		if {[file exists $file]} {
			lappend indelreslist -resource $res:$file
		}
	}
	set snpreslist {}
	foreach {pattern res} {
		hapmap*.sites.vcf hapmap,known=false,training=true,truth=true,prior=15.0
		*omni*.sites.vcf omni,known=false,training=true,truth=true,prior=12.0
		1000G_phase1.snps.high_confidence.*.sites.vcf 1000G,known=false,training=true,truth=false,prior=10.0
		*dbsnp*.vcf dbsnp,known=true,training=false,truth=false,prior=2.0
	} {
		set file [jobglob1 -checkcompressed 1 $gatkres/$pattern]
		if {[file exists $file]} {
			lappend snpreslist -resource $res:$file
		}
	}
	if {![llength $indelreslist]} {
		error "No INDEL recalibration resources found. (use -gatkres option?)"
	}
	if {![llength $snpreslist]} {
		error "No SNP recalibration resources found. (use -gatkres option?)"
	}
	# Variant (Quality Score) Recalibration
	# build indel model
	job [job_relfile2name gatkvqsr- $resultvcf].indel.recal -skip $resultvcf -deps {
		$rawvcf $gatkrefseq
	} -targets {
		$resultvcf.indel.recal $resultvcf.indel.recal.idx $resultvcf.indel.tranches
	} -vars {
		rawvcf gatkrefseq indelreslist indelannotations resultvcf maxmem
	} -code {
		gatkexec [list -XX:ParallelGCThreads=1 -Xms1g -Xmx${maxmem}g] VariantRecalibrator \
			-R $gatkrefseq \
			-V $rawvcf \
			-O $resultvcf.indel.recal.temp \
			--tranches-file $resultvcf.indel.tranches.temp \
			{*}$indelreslist \
			-an {*}[join $indelannotations " -an "] \
			-mode INDEL \
			-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
			--max-gaussians 4 \
			>@ stdout 2>@stderr
			file rename -- $resultvcf.indel.recal.temp $resultvcf.indel.recal
			file rename -- $resultvcf.indel.recal.temp.idx $resultvcf.indel.recal.idx
			file rename -- $resultvcf.indel.tranches.temp $resultvcf.indel.tranches
	}
	# gatk_IndexFeatureFile_job $resultvcf.indel.recal
	# apply indel model
	job [job_relfile2name gatkvqsr- $resultvcf.indelscalibrated] -skip $resultvcf -deps {
		$rawvcf $gatkrefseq $resultvcf.indel.recal $resultvcf.indel.recal.idx $resultvcf.indel.tranches
	} -targets {
		$rawvcf.indelscalibrated
	} -vars {
		rawvcf gatkrefseq resultvcf maxmem
	} -code {
		gatkexec [list -XX:ParallelGCThreads=1 -Xms1g -Xmx${maxmem}g] ApplyVQSR \
			-R $gatkrefseq \
			-V $rawvcf \
			-O $rawvcf.indelscalibrated.temp \
			--truth-sensitivity-filter-level 99.0 \
			--recal-file $resultvcf.indel.recal \
			--tranches-file $resultvcf.indel.tranches \
			-mode INDEL \
			>@ stdout 2>@stderr
			file rename -- $rawvcf.indelscalibrated.temp $rawvcf.indelscalibrated
	}
	#
	# build snp model
	job [job_relfile2name gatkvqsr- $resultvcf.snp.recal] -skip $resultvcf -deps {
		$rawvcf $gatkrefseq
	} -targets {
		$resultvcf.snp.recal $resultvcf.snp.recal.idx $resultvcf.snp.tranches
	} -vars {
		rawvcf gatkrefseq snpreslist snpannotations resultvcf maxmem
	} -code {
		gatkexec [list -XX:ParallelGCThreads=1 -Xms1g -Xmx${maxmem}g] VariantRecalibrator \
			-R $gatkrefseq \
			-V $rawvcf \
			-O $resultvcf.snp.recal.temp \
			--tranches-file $resultvcf.snp.tranches.temp \
			{*}$snpreslist \
			-an {*}[join $snpannotations " -an "] \
			-mode SNP \
			-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
			>@ stdout 2>@stderr
		file rename -- $resultvcf.snp.recal.temp $resultvcf.snp.recal
		file rename -- $resultvcf.snp.recal.temp.idx $resultvcf.snp.recal.idx
		file rename -- $resultvcf.snp.tranches.temp $resultvcf.snp.tranches
	}
	# gatk_IndexFeatureFile_job $resultvcf.snp.recal
	# apply snp model
	job [job_relfile2name gatkvqsr- $resultvcf] -skip $resultvcf -deps {
		$rawvcf.indelscalibrated $gatkrefseq $resultvcf.snp.recal $resultvcf.snp.recal.idx $resultvcf.snp.tranches
	} -targets {
		$resultvcf
	} -vars {
		rawvcf gatkrefseq resultvcf maxmem
	} -code {
		gatkexec [list -XX:ParallelGCThreads=1 -Xms1g -Xmx${maxmem}g] ApplyVQSR \
			-R $gatkrefseq \
			-V $rawvcf.indelscalibrated \
			-O $resultvcf.temp[file extension $resultvcf] \
			--truth-sensitivity-filter-level 99.0 \
			--recal-file $resultvcf.snp.recal \
			--tranches-file $resultvcf.snp.tranches \
			-mode SNP \
			>@ stdout 2>@stderr
		file rename -- $resultvcf.temp[file extension $resultvcf] $resultvcf
	}
}

proc cg_gatk_vqsr {args} {
	set args [job_init {*}$args]
	gatk_vqsr_job {*}$args
	job_wait
}


