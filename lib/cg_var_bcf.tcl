proc var_bcf_tools {} {
	return {samtools bcftools}
}

# identical to var_sam, but
proc var_bcf_job {args} {
	upvar job_logdir job_logdir
	set cmdline [clean_cmdline cg var_bcf {*}$args]
	set pre ""
	set opts {}
	set split 1
	set deps {}
	set regionfile {}
	set BAQ 1
	set BQ 0
	set cleanup 1
	set regmincoverage 3
	set threads 2
	set mincoverage 5
	set minqual 30
	set resultfiles 0
	set rootname {}
	set skips {}
	set callmethod m
	set resultfile {}
	set mem 5G
	set time 2:00:00
	cg_options var_bcf args {
		-l - deps {
			lappend deps $value
		}
		-regionfile {
			set regionfile $value
			lappend deps $value
		}
		-regmincoverage {
			set regmincoverage $value
		}
		-callmethod {
			set callmethod [string index $value 0]
			if {$callmethod ni "c m"} {
				error "callmethod must be one of: c (consensus) or m (multiallelic)"
			}
		}
		-pre {
			set pre $value
		}
		-split {
			set split $value
		}
		-BQ {
			set BQ $value
		}
		-BAQ {
			set BAQ $value
		}
		-mincoverage {
			set mincoverage $value
		}
		-minqual {
			set minqual $value
		}
		-threads {
			set threads $value
			# not used yet
		}
		-cleanup {
			set cleanup $value
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
	set refseq [refseq $refseq]
	if {$resultfile eq ""} {
		set resultfile [file dir $bamfile]/${pre}var-bcf-[file_rootname $bamfile].tsv.zst
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
	set uvarfile $destdir/${pre}uvar-$root.tsv.zst
	set sregfile $destdir/${pre}sreg-$root.tsv.zst
	set varallfile $destdir/${pre}varall-$root.tsv.zst
	set regclusterfile $destdir/reg_cluster-$root.tsv
	set resultlist [list $varfile $sregfile $varallfile {} $regclusterfile]
	lappend skips -skip $resultlist
	if {$resultfiles} {
		return $resultlist
	}
	if {$regionfile ne ""} {
		set regionfile [file_absolute $regionfile]
	} else {
		set regionfile [bam2reg_job {*}$skips -mincoverage $regmincoverage $bamfile]
	}
	# logfile
	job_logfile $destdir/var_bcf_$resulttail $destdir $cmdline \
		{*}[versions bwa bowtie2 samtools bcftools picard java gnusort8 zst os]
	# start
	# make sure reference sequence is indexed
	job ${pre}var_bcf_faidx {*}$skips -deps {$refseq} -targets {$refseq.fai} -code {
		exec samtools faidx $dep
	}
	set deps [list $bamfile $refseq $refseq.fai {*}$deps]
	set cache [file dir $varallfile]/cache_varall_bcf_[file tail $refseq].temp
	job_cleanup_add $cache
	job bcfvarall_${pre}varall-$root {*}$skips -deps $deps -cores $threads -mem $mem -time $time -targets {
		$varallfile
	} -vars {
		refseq opts BQ BAQ regionfile root threads callmethod split cache
	} -code {
		analysisinfo_write $dep $target analysis $root sample $root varcaller samtools varcaller_version [version samtools] varcaller_cg_version [version genomecomb] varcaller_region [filename $regionfile]
		set emptyreg [reg_isempty $regionfile]
		if {$emptyreg && [file exists $cache]} {
			file_copy $cache $target
		} else {
			if {$regionfile ne ""} {
				set bedfile [tempbed $regionfile $refseq]
				lappend opts -T $bedfile
			}
			if {![true $BAQ]} {
				lappend opts --no-BAQ
			}
			if {[catch {version samtools 1}]} {
				error "bcftools calling needs samtools v > 1"
			} else {
				# bcftools -v for variant only
				# -t DP: Number of high-quality bases (per sample)
				# -t SP: Phred-scaled strand bias P-value
				# exec samtools mpileup --uncompressed -t DP,SP --min-BQ $BQ --fasta-ref $refseq {*}$opts $dep 2>@ stderr | bcftools call --threads $threads -$callmethod - > $target.temp 2>@ stderr
				exec bcftools mpileup -Ou --fasta-ref $refseq \
					--count-orphans --max-depth 1000 --min-BQ $BQ \
					--annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP \
					{*}$opts $dep 2>@ stderr \
					| bcftools call \
						--format-fields GQ,GP \
						-$callmethod \
						--threads $threads \
					| cg vcf2tsv -skiprefindels 1 \
						-split $split \
						-meta [list refseq [file tail $refseq]] \
						-removefields {name filter AN AC AF AA INDEL G3 HWE CLR UGT CGT PCHI2 QCHI2 PR} \
					{*}[compresspipe $target] > $target.temp 2>@ stderr
			}
			file rename -force -- $target.temp $target
			if {$emptyreg && ![file exists $cache]} {
				file_copy $target $cache
			}
		}
	}
	zstindex_job {*}$skips $varallfile
	job ${pre}var-$root {*}$skips -deps {
		$varallfile
	} -targets {
		$uvarfile
	} -skip {
		${pre}var-$root.tsv ${pre}var-$root.tsv.analysisinfo
	} -vars {
		root pre
	} -code {
		analysisinfo_write $dep $target analysis $root sample $root varcaller_mincoverage 5 varcaller_minquality 30 varcaller_cg_version [version genomecomb]
		set tempfile [filetemp_ext $target]
		cg select -overwrite 1 -q {
			$alt ne "." && $alleleSeq1 ne "." && $quality >= 10 && $totalcoverage > 4
			&& $zyg ni "r o"
		} -f {
			chromosome begin end type ref alt quality alleleSeq1 alleleSeq2
			{sequenced=if($quality < 30 || $totalcoverage < 5,"u","v")}
			{zyg=if($quality < 30 || $totalcoverage < 5,"u",$zyg)}
			*
		} $dep $tempfile
		result_rename $tempfile $target
 	}
	# annotvar_clusters_job works using jobs
	annotvar_clusters_job {*}$skips $uvarfile $resultfile
	# find regions
	sreg_sam_job ${pre}sreg-$root $varallfile $sregfile $mincoverage $minqual $skips
	# cleanup
	if {$cleanup} {
		set cleanupfiles [list \
			$uvarfile [gzroot $uvarfile].index [gzroot $uvarfile].temp \
			[file root $varallfile].vcf \
			[file root $varallfile].vcf.idx \
		]
		set cleanupdeps [list $resultfile $varallfile]
		cleanup_job clean_${pre}var-$root $cleanupfiles $cleanupdeps
	}
	return $resultlist
}

proc cg_var_bcf {args} {
	set args [job_init {*}$args]
	var_bcf_job {*}$args
	job_wait
}
