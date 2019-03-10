proc sreg_sam_job {job varallfile resultfile {skips {}}} {
	upvar job_logdir job_logdir
	job $job {*}$skips -deps {$varallfile} -targets {$resultfile} -code {
		set temp [filetemp $target]
		set temp2 [filetemp $target]
		cg select -overwrite 1 -q {$quality >= 30 && $totalcoverage >= 5 && $type ne "ins"} -f {chromosome begin end} $dep $temp
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
		-skip {
			lappend skips -skip $value
		}
		default {
			lappend opts $key $value
		}
	} {bamfile refseq}
	set bamfile [file_absolute $bamfile]
	set refseq [refseq $refseq]
	set destdir [file dir $bamfile]
	if {$rootname eq ""} {
		set root sam-[file_rootname $bamfile]
	} else {
		set root $rootname
	}
	# resultfiles
	set varfile ${pre}var-$root.tsv
	set sregfile ${pre}sreg-$root.tsv
	set varallfile ${pre}varall-$root.tsv
	set resultlist [list $destdir/$varfile.zst $destdir/$sregfile.zst $destdir/$varallfile.zst $destdir/reg_cluster-$root.tsv.zst]
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
	job_logfile $destdir/var_sam_[file tail $bamfile] $destdir $cmdline \
		{*}[versions bwa bowtie2 samtools picard java gnusort8 zst os]
	set file [file tail $bamfile]
	# start
	set keeppwd [pwd]
	cd $destdir
	# make sure reference sequence is indexed
	job ${pre}var_sam_faidx {*}$skips -deps {$refseq} -targets {$refseq.fai} -code {
		exec samtools faidx $dep
	}
	set deps [list $file $refseq $refseq.fai {*}$deps]
	job ${pre}varall-$root {*}$skips -deps $deps \
	  -targets {${pre}varall-$root.vcf ${pre}varall-$root.vcf.analysisinfo} \
	  -skip [list ${pre}varall-$root.tsv ${pre}varall-$root.tsv.analysisinfo] \
	  -vars {refseq opts BQ regionfile root} \
	  -code {
		analysisinfo_write $dep $target sample $root varcaller samtools varcaller_version [version samtools] varcaller_cg_version [version genomecomb] varcaller_region [filename $regionfile]
		if {$regionfile ne ""} {
			set bedfile [tempbed $regionfile $refseq]
			lappend opts -l $bedfile
		}
		if {[catch {version samtools 1}]} {
			exec samtools mpileup -uDS -Q $BQ -f $refseq {*}$opts $dep 2>@ stderr | bcftools view -cg - > $target.temp 2>@ stderr
		} else {
			# bcftools -v for variant only
			# -t DP: Number of high-quality bases (per sample)
			# -t SP: Phred-scaled strand bias P-value
			exec samtools mpileup --uncompressed -t DP,SP --min-BQ $BQ --fasta-ref $refseq {*}$opts $dep 2>@ stderr | bcftools call -c - > $target.temp 2>@ stderr
		}
		file rename -force $target.temp $target
	}
	job ${pre}varall-sam2tsv-$root {*}$skips -deps {${pre}varall-$root.vcf} \
	-targets {${pre}varall-$root.tsv.zst ${pre}varall-$root.tsv.analysisinfo} -vars {split refseq} -code {
		analysisinfo_write $dep $target
		cg vcf2tsv -skiprefindels 1 -split $split -meta [list refseq [file tail $refseq]] -removefields {name filter AN AC AF AA INDEL G3 HWE CLR UGT CGT PCHI2 QCHI2 PR} $dep $target.temp.zst
		file rename -force $target.temp.zst $target
	}
	# zst_job ${pre}varall-$root.tsv -i 1
	zstindex_job {*}$skips ${pre}varall-$root.tsv.zst
	job ${pre}var-$root {*}$skips -deps {${pre}varall-$root.tsv} \
	-targets {${pre}uvar-$root.tsv ${pre}uvar-$root.tsv.analysisinfo} \
	-skip [list ${pre}var-$root.tsv ${pre}var-$root.tsv.analysisinfo] -vars {root pre} \
	-code {
		analysisinfo_write $dep $target varcaller_mincoverage 5 varcaller_minquality 30 varcaller_cg_version [version genomecomb]
		cg select -q {
			$alt ne "." && $alleleSeq1 ne "." && $quality >= 10 && $totalcoverage > 4
			&& $zyg ni "r o"
		} -f {
			chromosome begin end type ref alt quality alleleSeq1 alleleSeq2
			{sequenced=if($quality < 30 || $totalcoverage < 5,"u","v")}
			{zyg=if($quality < 30 || $totalcoverage < 5,"u",$zyg)}
			*
		} $dep $target.temp
		file rename -force $target.temp $target
	}
	# annotvar_clusters_job works using jobs
	annotvar_clusters_job {*}$skips ${pre}uvar-$root.tsv ${pre}var-$root.tsv.zst
	# find regions
	sreg_sam_job ${pre}sreg-$root ${pre}varall-$root.tsv ${pre}sreg-$root.tsv.zst $skips
	# cleanup
	if {$cleanup} {
		set cleanupfiles [list \
			${pre}uvar-$root.tsv ${pre}uvar-$root.tsv.index ${pre}uvar-$root.tsv.analysisinfo \
			${pre}varall-$root.vcf ${pre}varall-$root.vcf.analysisinfo \
			${pre}varall-$root.vcf.idx \
		]
		set cleanupdeps [list ${pre}var-$root.tsv ${pre}varall-$root.tsv]
		cleanup_job clean_${pre}var-$root $cleanupfiles $cleanupdeps
	}
	cd $keeppwd
	return $resultlist
}

proc cg_var_sam {args} {
	set args [job_init {*}$args]
	var_sam_job {*}$args
	job_wait
}
