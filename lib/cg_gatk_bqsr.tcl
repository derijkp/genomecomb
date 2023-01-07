# info why not to use bqsr:
# Impact of post-alignment processing in variant discovery from whole exome data
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5048557/

proc gatk_bqsr_job args {
	upvar job_logdir job_logdir
	set dbdir {}
	set maxmem 4
	set distrreg 0
	set cmdline [clean_cmdline cg gatk_bqsr {*}$args]
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
		-distrreg {
			set distrreg [distrreg_checkvalue $value]
		}
	} {srcbam resultbam} 2 2 {
		perform Base (Quality Score) Recalibration on bam
	}
	set dbdir [dbdir $dbdir]
	if {![info exists gatkres]} {
		set gatkres $dbdir/gatkres
	}
	set refseq [lindex [glob $dbdir/genome_*.ifas] 0]
	set gatkrefseq [gatk_refseq_job $refseq]
	set srcbam [file_absolute $srcbam]
	set resultbam [file_absolute $resultbam]
	if {![info exists job_logdir]} {
		set_job_logdir [file dir $resultbam]/log_jobs
	}
	job_logfile [file dir $resultbam]/gatk_bqsr_[file tail $resultbam] [file dir $resultbam] $cmdline \
		{*}[versions gatk]
	#
	# Base (Quality Score) Recalibration resources and settings
	set knownsitesopt {}
	set deps [list $srcbam $gatkrefseq]
	foreach pattern {*dbsnp*.vcf hapmap*vcf} {
		set knownvcf [gzfile $dbdir/gatkres/$pattern]
		if {[file exists $knownvcf]} {
			lappend knownsitesopt --known-sites $knownvcf
			lappend deps $knownvcf
		}
	}
	if {![llength $knownsitesopt]} {
		error "knownsites file no found: searched $dbdir/gatkres"
	}
	# Base (Quality Score) Recalibration
	# build model
	if {$distrreg eq "0"} {
		job [job_relfile2name gatkbqsr- $resultbam.table] \
		-skip $resultbam \
		-mem ${maxmem}G \
		-deps $deps \
		-targets {
			$resultbam.bqsr.table
		} -vars {
			srcbam gatkrefseq resultbam knownsitesopt maxmem
		} -code {
			gatkexec [list -XX:ParallelGCThreads=1 -Xms1g -Xmx${maxmem}g] BaseRecalibrator \
				{*}$knownsitesopt \
				-R $gatkrefseq \
				-I $srcbam \
				-O $resultbam.bqsr.table.temp \
				>@ stdout 2>@stderr
			file rename -- $resultbam.bqsr.table.temp $resultbam.bqsr.table
		}
	} else {
		set workdir [workdir $resultbam.bqsr.table]
		file mkdir $workdir
		set regions [list_remove [distrreg_regs $distrreg $refseq] unaligned]
		set list {}
		foreach region $regions {
			set target $workdir/bqsr.table.$region
			lappend list $target
			set lopts {}
			foreach r [distrreg_reg2gatklist $refseq $region] {
				lappend lopts -L $r
			}
			job [job_relfile2name gatkbqsr- $resultbam.table.$region] \
			-skip $resultbam \
			-mem ${maxmem}G \
			-deps $deps \
			-targets {
				$target
			} -vars {
				srcbam gatkrefseq lopts resultbam knownsitesopt maxmem region
			} -code {
				gatkexec [list -XX:ParallelGCThreads=1 -Xms1g -Xmx${maxmem}g] BaseRecalibrator \
					{*}$knownsitesopt \
					-R $gatkrefseq \
					{*}$lopts \
					-I $srcbam \
					-O $target.temp \
					>@ stdout 2>@stderr
				file rename -- $target.temp $target
			}
		}
		job [job_relfile2name gatkbqsr- $resultbam.table] \
		-skip $resultbam \
		-mem ${maxmem}G \
		-deps $list \
		-targets {
			$resultbam.bqsr.table
		} -vars {
			resultbam maxmem workdir
		} -code {
			set inputopts {}
			foreach file $deps {
				lappend inputopts --input $file
			}
			gatkexec [list -XX:ParallelGCThreads=1 -Xms1g -Xmx${maxmem}g] GatherBQSRReports \
				{*}$inputopts \
				-O $workdir/result.bqsr.table \
				>@ stdout 2>@stderr
			file rename -- $workdir/result.bqsr.table $resultbam.bqsr.table
		}
	}
	# apply model
	job [job_relfile2name gatkbqsr- $resultbam] \
	-skip $resultbam \
	-mem ${maxmem}G -deps {
		$srcbam $resultbam.bqsr.table
	} -targets {
		$resultbam
	} -vars {
		srcbam gatkrefseq resultbam maxmem
	} -code {
		gatkexec [list -XX:ParallelGCThreads=1 -Xms1g -Xmx${maxmem}g] ApplyBQSR \
			-R $gatkrefseq \
			-I $srcbam \
			--bqsr-recal-file $resultbam.bqsr.table \
			-O $resultbam.temp[file extension $resultbam] \
			>@ stdout 2>@stderr
		file rename -- $resultbam.temp[file extension $resultbam] $resultbam
	}
}

proc cg_gatk_bqsr {args} {
	set args [job_init {*}$args]
	gatk_bqsr_job {*}$args
	job_wait
}



