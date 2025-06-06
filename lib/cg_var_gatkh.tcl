proc validate_var_gatkh {refseq distrreg datatype} {
	# seperate because command is not the same as cmd
	if {[version gatk] eq "?"} {
		error "command not available, make sure gatk is installed, e.g. using \"cg install gatk\""
	}
}
proc var_gatkh_tools {} {
	return {gatkh}
}

proc sreg_gatkh_job {job varallfile resultfile {mincoverage 8} {mingenoqual 25} {skips {}}} {
	upvar job_logdir job_logdir
	job $job {*}$skips -deps {
		$varallfile
	} -targets {
		$resultfile
	} -vars {
		mincoverage mingenoqual
	} -code {
		set temp [filetemp $target]
		exec cg vcf2tsv $dep \
			| cg select -q [subst {
				\$genoqual >= $mingenoqual && \$coverage >= $mincoverage && \$type ne "ins"
			}] -f {chromosome begin end} -s - \
			| cg regjoin {*}[compresspipe $target] > $temp
		file rename -force -- $temp $target
		cg_zindex $target
	}
}

proc var_gatkh_job {args} {
	# putslog [list var_gatkh_job {*}$args]
	set cmdline [clean_cmdline cg var_gatkh {*}$args]
	global appdir
	upvar job_logdir job_logdir
	set pre ""
	set opts {}
	set split 1
	set deps {}
	set regionfile {}
	set threads 2
	set cleanup 1
	set regmincoverage 3
	set ERC GVCF
	set mincoverage 8
	set mingenoqual 25
	set resultfiles 0
	set rootname {}
	set skips {}
	set resultfile {}
	set mem 15G
	set time 3:00:00
	cg_options var_gatkh args {
		-L - -deps {
			lappend deps [file_absolute $value]
		}
		-regionfile {
			set regionfile $value
		}
		-regmincoverage {
			set regmincoverage $value
		}
		-pre {
			set pre $value
		}
		-split {
			set split $value
		}
		-threads - -t {
			putslog "-threads $value ignored because gatkh does not support threads for now"
			set threads $value
		}
		-cleanup {
			set cleanup $value
		}
		-ERC - -emitRefConfidence {
			if {$value ni {BP_RESOLUTION GVCF NONE}} {error "option $value not supported for -ERC, must be one of: BP_RESOLUTION, GVCF, NONE"}
			set ERC $value
		}
		-mincoverage {
			set mincoverage $value
		}
		-mingenoqual {
			set mingenoqual $value
		}
		-resultfiles {
			set resultfiles $value
		}
		-rootname {
			set rootname $value
		}
		-datatype {
			# not actually used
		}
		-skip {
			lappend skips -skip $value
		}
		-opts {
			set opts $value
		}
		-mem {
			set mem $value
		}
		-time {
			set time $value
		}
	} {bamfile refseq resultfile} 2 3
	set bamfile [file_absolute $bamfile]
	set refseq [file_absolute $refseq]
	if {$resultfile eq ""} {
		set resultfile [file dir $bamfile]/${pre}var-gatkh-[file_rootname $bamfile].tsv.zst
	} else {
		set resultfile [file_absolute $resultfile]
	}
	set resulttail [file tail $resultfile]
	set destdir [file dir $resultfile]
	if {$rootname eq ""} {
		set root [file_rootname $resultfile]
	} else {
		set root $rootname
	}
	# resultfiles
	set varfile $resultfile
	set vcffile [file root [gzroot $varfile]].vcf.gz
	set varallfile $destdir/${pre}varall-$root.gvcf.gz
	set uvarfile $destdir/${pre}uvar-$root.tsv.zst
	set sregfile $destdir/${pre}sreg-$root.tsv.zst
	set regclusterfile $destdir/reg_cluster-$root.tsv.zst
	set resultlist [list $varfile $sregfile $varallfile $vcffile $regclusterfile]
	if {$resultfiles} {
		return $resultlist
	}
	lappend skips -skip $resultlist
	if {$regionfile ne ""} {
		set regionfile [file_absolute $regionfile]
	} else {
		set regionfile [bam2reg_job {*}$skips -mincoverage $regmincoverage $bamfile]
	}
	lappend deps $regionfile
	# logfile
	job_logfile $destdir/var_gatkh_$resulttail $destdir $cmdline \
		{*}[versions bwa bowtie2 samtools gatk picard java gnusort8 zst os]
	# start
	## Produce gatkh SNP calls
	set gatkrefseq [gatk_refseq_job $refseq]
	set dep $bamfile
	set bamindex $bamfile.[indexext $bamfile]
	set deps [list $bamfile $gatkrefseq $bamindex {*}$deps]
	set cache [file dir $varallfile]/cache_varall_gatkh_[file tail $refseq].temp
	job_cleanup_add $cache
	job $varallfile {*}$skips -mem $mem -time $time -deps $deps -targets {
		$varallfile $varallfile.tbi
	} -vars {
		opts regionfile gatkrefseq refseq root ERC varallfile cache
	} -code {
		analysisinfo_write $dep $varallfile analysis $root sample $root varcaller gatkh varcaller_version [version gatk] varcaller_cg_version [version genomecomb] varcaller_region [filename $regionfile]
		set emptyreg [reg_isempty $regionfile]
		if {$emptyreg && [file exists $cache]} {
			copywithindex $cache $target $cache.tbi.temp $target.tbi
		} else {
			if {$regionfile ne ""} {
				set bedfile [tempbed $regionfile $refseq]
				lappend opts -L $bedfile
			}
			# -finishedpattern is a hack to catch an error that sometimes seems to happen, after fully processing the data
			# Do not use redirect to stdout/stderr, as the code needs the output to check if actual analysis was finished
			gatkexec -finishedpattern {HaplotypeCaller done\. Elapsed time} {-XX:ParallelGCThreads=1 -Xms512m -Xmx4g} HaplotypeCaller \
				{*}$opts -R $gatkrefseq \
				-I $dep \
				-O $varallfile.temp.gz \
				--annotate-with-num-discovered-alleles \
				-ERC $ERC \
				-G StandardAnnotation \
				-G StandardHCAnnotation \
				-G AS_StandardAnnotation \
				--max-reads-per-alignment-start 0
			file rename -force -- $varallfile.temp.gz $varallfile
			file rename -force -- $varallfile.temp.gz.tbi $varallfile.tbi
			# file delete $varallfile.temp
			if {$emptyreg && ![file exists $cache]} {
				file_copy -force $target $cache
				file_copy -force $target.tbi $cache.tbi.temp
			}
		}
	}
	job ${pre}gvcf2tsv-$root {*}$skips -deps {
		$varallfile $gatkrefseq
	} -targets {
		$uvarfile $vcffile
	} -vars {
		sample split pre root gatkrefseq varallfile mincoverage mingenoqual refseq vcffile uvarfile
	} -skip {
		$varfile
	} -code {
		set tempvcf [filetemp_ext $vcffile]
		analysisinfo_write $dep $target varcaller_mincoverage $mincoverage varcaller_mingenoqual $mingenoqual varcaller_cg_version [version genomecomb]
		gatkexec {-XX:ParallelGCThreads=1 -Xms512m -Xmx4g} GenotypeGVCFs \
			-R $gatkrefseq \
			-V $varallfile \
			-O $tempvcf \
			-G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation
		catch {file delete [gzroot $tempvcf].idx}
		result_rename $tempvcf $vcffile
		set fields {chromosome begin end type ref alt quality alleleSeq1 alleleSeq2}
		lappend fields [subst {sequenced=if(\$genoqual < $mingenoqual || \$coverage < $mincoverage,"u","v")}]
		lappend fields [subst {zyg=if(\$genoqual < $mingenoqual || \$coverage < $mincoverage,"u",\$zyg)}]
		lappend fields *
		exec cg vcf2tsv -split $split -meta [list refseq [file tail $refseq]] -removefields {
			name filter AN AC AF AA ExcessHet InbreedingCoeff MLEAC MLEAF NDA RPA RU STR
		} $vcffile | cg select -f $fields | cg zst > $uvarfile.temp
		file rename -force -- $uvarfile.temp $uvarfile
	}
	maketabix_job $vcffile
	# annotvar_clusters_job works using jobs
	annotvar_clusters_job {*}$skips -deletesrc $cleanup $uvarfile $varfile
	# make sreg
	sreg_gatkh_job ${pre}sreg-$root $varallfile $sregfile $mincoverage $mingenoqual $skips
	return $resultlist
}

proc cg_var_gatkh {args} {
	set args [job_init {*}$args]
	var_gatkh_job {*}$args
	job_wait
}
