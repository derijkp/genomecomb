proc var_sam_tools {} {
	return {samtools bcftools}
}

proc sreg_sam_job {job varallfile resultfile {mincoverage 5} {minqual 30} {skips {}}} {
	upvar job_logdir job_logdir
	job $job {*}$skips -deps {
		$varallfile
	} -targets {
		$resultfile
	} -vars {
		mincoverage minqual 
	} -code {
		set temp [filetemp $target]
		set temp2 [filetemp $target]
		cg select -overwrite 1 -q [subst {
			\$quality >= $minqual && \$totalcoverage >= $mincoverage && \$type ne "ins"
		}] -f {chromosome begin end} $dep $temp
		file_write $temp2 "# regions selected from [gzroot $dep]: \$quality >= 30 && \$totalcoverage >= 5\n"
		cg regjoin $temp >> $temp2
		compress $temp2 $target 1 0
		file delete $temp
	}
}

proc var_sam_job {args} {
	upvar job_logdir job_logdir
	set pre ""
	set opts {}
	set split 1
	set deps {}
	set regionfile {}
	set BQ 0
	set cleanup 1
	set regmincoverage 3
	set threads 2
	set resultfiles 0
	set rootname {}
	set skips {}
	set callmethod c
	set resultfile {}
	cg_options var_sam args {
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
	} {bamfile refseq resultfile} 2 3
	set bamfile [file_absolute $bamfile]
	set refseq [refseq $refseq]
	if {$resultfile eq ""} {
		set resultfile [file dir $bamfile]/${pre}var-sam-[file_rootname $bamfile].tsv.zst
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
	set vcffile {}
	set varallfile $destdir/${pre}varall-$root.tsv.zst
	set varallvcf $destdir/${pre}varall-$root.vcf.zst
	set uvarfile $destdir/${pre}uvar-$root.tsv.zst
	set sregfile $destdir/${pre}sreg-$root.tsv.zst
	set regclusterfile $destdir/reg_cluster-$root.tsv
	set resultlist [list $varfile $sregfile $varallfile $vcffile $regclusterfile]
	if {$resultfiles} {
		return $resultlist
	}
	if {$regionfile ne ""} {
		set regionfile [file_absolute $regionfile]
	} else {
		set regionfile [bam2reg_job {*}$skips -mincoverage $regmincoverage $bamfile]
	}
	# logfile
	set cmdline [list cg var_sam]
	foreach option {
		split deps bed pre BQ regionfile regmincoverage
	} {
		if {[info exists $option]} {
			lappend cmdline -$option [get $option]
		}
	}
	lappend cmdline {*}$opts $bamfile $refseq
	job_logfile $destdir/var_sam_[file tail $resultfile] $destdir $cmdline \
		{*}[versions bwa bowtie2 samtools bcftools picard java gnusort8 zst os]
	# start
	# make sure reference sequence is indexed
	job ${pre}var_sam_faidx {*}$skips -deps {$refseq} -targets {$refseq.fai} -code {
		exec samtools faidx $dep
	}
	set deps [list $bamfile $refseq $refseq.fai {*}$deps]
	job ${pre}varall-$root {*}$skips -deps $deps -cores $threads -targets {
		$varallvcf
	} -skip {
		$varallfile
	} -vars {
		refseq opts BQ regionfile root threads callmethod
	} -code {
		analysisinfo_write $dep $target sample $root varcaller samtools varcaller_version [version samtools] varcaller_cg_version [version genomecomb] varcaller_region [filename $regionfile]
		set emptyreg [reg_isempty $regionfile]
		set cache [file dir $target]/cache_var_gatk_[file tail $refseq].temp.zst
		if {$emptyreg && [file exists $cache]} {
			copywithindex $cache $target $cache.tbi.temp $target.tbi
		} else {
			if {$regionfile ne ""} {
				set bedfile [tempbed $regionfile $refseq]
				lappend opts -l $bedfile
			}
			if {[catch {version samtools 1}]} {
				catch_exec samtools mpileup -uDS -Q $BQ -f $refseq {*}$opts $dep | bcftools view -cg - {*}[compresspipe $target] > $target.temp.zst
			} else {
				# bcftools -v for variant only
				# -t DP: Number of high-quality bases (per sample)
				# -t SP: Phred-scaled strand bias P-value
				catch_exec samtools mpileup --uncompressed -t DP,SP --min-BQ $BQ --fasta-ref $refseq {*}$opts $dep | bcftools call --threads $threads -$callmethod - {*}[compresspipe $target] > $target.temp.zst
			}
			file rename -force -- $target.temp.zst $target
			if {$emptyreg && ![file exists $cache]} {
				file copy -force $target $cache
			}
		}
	}
	job ${pre}varall-sam2tsv-$root {*}$skips -deps {
		$varallvcf
	} -targets {
		$varallfile
	} -vars {
		split refseq
	} -code {
		analysisinfo_write $dep $target
		cg vcf2tsv -skiprefindels 1 -split $split -meta [list refseq [file tail $refseq]] -removefields {name filter AN AC AF AA INDEL G3 HWE CLR UGT CGT PCHI2 QCHI2 PR} $dep $target.temp.zst
		file rename -force -- $target.temp.zst $target
	}
	# zst_job ${pre}varall-$root.tsv -i 1
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
		analysisinfo_write $dep $target varcaller_mincoverage 5 varcaller_minquality 30 varcaller_cg_version [version genomecomb]
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
		file rename -force -- $tempfile $target
	}
	# annotvar_clusters_job works using jobs
	annotvar_clusters_job {*}$skips $uvarfile $varfile
	# find regions
	sreg_sam_job ${pre}sreg-$root $varallfile $sregfile 5 30 $skips
	# cleanup
	if {$cleanup} {
		set cleanupfiles [list \
			$uvarfile [gzroot $uvarfile].index [gzroot $uvarfile].temp \
			$varallvcf \
			[gzroot $varallvcf].idx \
		]
		set cleanupdeps [list $varfile $varallfile]
		cleanup_job clean_${pre}var-$root $cleanupfiles $cleanupdeps
	}
	return $resultlist
}

proc cg_var_sam {args} {
	set args [job_init {*}$args]
	var_sam_job {*}$args
	job_wait
}
