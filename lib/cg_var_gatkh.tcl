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
	} {bamfile refseq}
	set bamfile [file_absolute $bamfile]
	set refseq [file_absolute $refseq]
	set destdir [file dir $bamfile]
	set bamtail [file tail $bamfile]
	if {$rootname eq ""} {
		set root gatkh-[file_rootname $bamtail]
	} else {
		set root $rootname
	}
	# resultfiles
	set varfile ${pre}var-$root.tsv
	set sregfile ${pre}sreg-$root.tsv
	set varallfile ${pre}varall-$root.gvcf
	set resultlist [list $destdir/$varfile.zst $destdir/$sregfile.zst $destdir/$varallfile.gz $destdir/reg_cluster-$root.tsv.zst]
	if {$resultfiles} {
		return $resultlist
	}
	if {$regionfile ne ""} {
		set regionfile [file_absolute $regionfile]
	} else {
		set regionfile [bam2reg_job {*}$skips -mincoverage $regmincoverage $bamfile]
	}
	lappend deps $regionfile
	# logfile
	set cmdline [list cg var_gatkh]
	foreach option {
		split deps bed pre regionfile regmincoverage
	} {
		if {[info exists $option]} {
			lappend cmdline -$option [get $option]
		}
	}
	lappend cmdline {*}$opts $bamfile $refseq
	job_logfile $destdir/var_gatkh_$bamtail $destdir $cmdline \
		{*}[versions bwa bowtie2 samtools gatk picard java gnusort8 zst os]
	# start
	## Produce gatkh SNP calls
	set keeppwd [pwd]
	cd $destdir
	set gatkrefseq [gatk_refseq_job $refseq]
	set dep $bamtail
	set bamtailindex $bamtail.[indexext $bamtail]
	set deps [list $bamtail $gatkrefseq $bamtailindex {*}$deps]
	job $varallfile {*}$skips -mem 15G -deps $deps -targets {
		$varallfile.gz $varallfile.gz.tbi
	} -vars {
		opts regionfile gatkrefseq refseq root ERC varallfile
	} -code {
		analysisinfo_write $dep $varallfile sample $root varcaller gatkh varcaller_version [version gatk] varcaller_cg_version [version genomecomb] varcaller_region [filename $regionfile]
		set emptyreg [reg_isempty $regionfile]
		set cache [file dir $target]/cache_var_gatkh_[file tail $refseq].temp
		if {$emptyreg && [file exists $cache]} {
			copywithindex $cache $target $cache.tbi.temp $target.tbi
		} else {
			if {$regionfile ne ""} {
				set bedfile [tempbed $regionfile $refseq]
				lappend opts -L $bedfile
			}
			# -finishedpattern is a hack to catch an error that sometimes seems to happen, after fully processing the data
			# Do not use redirect to stdout/stderr, as the code needs the output to check if actual analysis was finished
			gatkexec -finishedpattern {HaplotypeCaller done\. Elapsed time} {-XX:ParallelGCThreads=1 -d64 -Xms512m -Xmx4g} HaplotypeCaller \
				{*}$opts -R $gatkrefseq \
				-I $dep \
				-O $varallfile.temp.gz \
				--annotate-with-num-discovered-alleles \
				-ERC $ERC \
				-G StandardAnnotation \
				-G StandardHCAnnotation \
				-G AS_StandardAnnotation \
				--max-reads-per-alignment-start 0
			file rename -force -- $varallfile.temp.gz $varallfile.gz
			file rename -force -- $varallfile.temp.gz.tbi $varallfile.gz.tbi
			# file delete $varallfile.temp
			if {$emptyreg && ![file exists $cache]} {
				file copy -force $target $cache
				file copy -force $target.tbi $cache.tbi.temp
			}
		}
	}
	job ${pre}gvcf2tsv-$root {*}$skips -deps {
		$varallfile.gz $gatkrefseq
	} -targets {
		${pre}uvar-$root.tsv
	} -vars {
		sample split pre root gatkrefseq varallfile mincoverage mingenoqual refseq
	} -skip {
		$varfile.zst
	} -code {
		analysisinfo_write $dep $target varcaller_mincoverage $mincoverage varcaller_mingenoqual $mingenoqual varcaller_cg_version [version genomecomb]
		gatkexec {-XX:ParallelGCThreads=1 -d64 -Xms512m -Xmx4g} GenotypeGVCFs \
			-R $gatkrefseq \
			-V $varallfile.gz \
			-O ${pre}uvar-$root.temp.vcf \
			-G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation
		catch {file delete ${pre}uvar-$root.temp.vcf.idx}
		file rename -force -- ${pre}uvar-$root.temp.vcf ${pre}uvar-$root.vcf
		set fields {chromosome begin end type ref alt quality alleleSeq1 alleleSeq2}
		lappend fields [subst {sequenced=if(\$genoqual < $mingenoqual || \$coverage < $mincoverage,"u","v")}]
		lappend fields [subst {zyg=if(\$genoqual < $mingenoqual || \$coverage < $mincoverage,"u",\$zyg)}]
		lappend fields *
		exec cg vcf2tsv -split $split -meta [list refseq [file tail $refseq]] -removefields {
			name filter AN AC AF AA ExcessHet InbreedingCoeff MLEAC MLEAF NDA RPA RU STR
		} ${pre}uvar-$root.vcf | cg select -f $fields > ${pre}uvar-$root.tsv.temp
		file rename -force -- ${pre}uvar-$root.tsv.temp ${pre}uvar-$root.tsv
	}
	# annotvar_clusters_job works using jobs
	annotvar_clusters_job {*}$skips ${pre}uvar-$root.tsv $varfile.zst
	# make sreg
	sreg_gatkh_job ${pre}sreg-$root $varallfile $sregfile.zst $mincoverage $mingenoqual $skips
	# cleanup
	if {$cleanup} {
		set cleanupfiles [list \
			${pre}var-$root.vcf \
			${pre}uvar-$root.tsv ${pre}uvar-$root.vcf ${pre}uvar-$root.tsv.index
		]
		set cleanupdeps [list $varfile $varallfile.gz]
		cleanup_job clean_${pre}var-$root $cleanupfiles $cleanupdeps
	}
	cd $keeppwd
	return $resultlist
}

proc cg_var_gatkh {args} {
	set args [job_init {*}$args]
	var_gatkh_job {*}$args
	job_wait
}
