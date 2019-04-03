proc gatk_refseq_job refseq {
	upvar job_logdir job_logdir
	set nrefseq [file root $refseq].fa
	if {![file exists $nrefseq] && $refseq ne $nrefseq} {
		mklink $refseq $nrefseq
	}
	job gatkrefseq_faidx-[file tail $nrefseq] -deps {$nrefseq} -targets {$nrefseq.fai} -code {
		exec samtools faidx $dep
	}
	set dict [file root $nrefseq].dict
	job gatkrefseq-[file tail $nrefseq] -deps {$nrefseq} -targets {$dict} -vars {nrefseq} -code {
		file delete $target.temp
		picard CreateSequenceDictionary R= $nrefseq O= $target.temp 2>@ stderr >@ stdout
		file rename -force $target.temp $target
	}
	return $nrefseq
}

proc annotvar_clusters_job {args} {
	upvar job_logdir job_logdir
	set skips {}
	cg_options annotvar_clusters args {
		-skip {
			lappend skips -skip $value
		}
	} {file resultfile} 2 2 {
		find clusters of variants
	}
	set root [file_rootname $file]
	job annotvar-clusters-$root {*}$skips -skip [list [gzroot $resultfile]] -deps {
		$file
	} -targets {
		reg_cluster-$root.tsv.zst reg_cluster-$root.tsv.zst.zsti
	} -code {
		analysisinfo_write $dep $target
		if {[file size $dep]} {
			cg clusterregions < $dep > $target.temp
			cg_zst $target.temp
			file rename -force $target.temp.zst $target
			cg_zstindex $target
		} else {
			file_write $target ""
			file_write $target.zsti ""
		}
	}
	job annotvar-annotclusters-$root {*}$skips -deps {
		$file reg_cluster-$root.tsv.zst
	} -targets {
		$resultfile
	} -code {
		analysisinfo_write $dep $target
		if {[file size $dep]} {
			cg annotate -stack 1 -v 2 -analysisinfo 0 $dep $target {*}[list_remove [lrange $deps 1 end] {}] 2>@ stderr >@ stdout
		} else {
			file_write $target ""
		}
	}
}

proc sreg_gatk_job {job varallfile resultfile {skips {}}} {
	upvar job_logdir job_logdir
	job $job {*}$skips -deps {
		$varallfile
	} -targets {
		$resultfile
	} -code {
		analysisinfo_write $dep $target
		set temp [filetemp $target]
		set temp2 [filetemp $target]
		cg select -overwrite 1 -q {$quality >= 30 && $totalcoverage >= 5 && $type ne "ins"} -f {chromosome begin end} $dep $temp
		file_write $temp2 "# regions selected from [gzroot $dep]: \$quality >= 30 && \$totalcoverage >= 5\n"
		cg regjoin $temp >> $temp2
		if {[file extension $target] eq ".zst"} {
			cg_zst -keep 0 -i 1 -o $target $temp2
		} else {
			file rename $temp2 $target
		}
		file delete $temp
	}
}

proc var_gatk_job {args} {
	upvar job_logdir job_logdir
	set pre ""
	set opts {}
	set split 1
	set deps {}
	set regionfile {}
	set threads 2
	set cleanup 1
	set regmincoverage 3
	set resultfiles 0
	set rootname {}
	set skips {}
	cg_options var_gatk args {
		-L - -deps {
			lappend deps $value
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
			set threads $value
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
	} {bamfile refseq}
	set bamfile [file_absolute $bamfile]
	set refseq [refseq $refseq]
	set destdir [file dir $bamfile]
	set bamtail [file tail $bamfile]
	if {$rootname eq ""} {
		set root gatk-[file_rootname $bamtail]
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
	lappend deps $regionfile
	# logfile
	set cmdline [list cg var_gatk]
	foreach option {
		split deps bed pre regionfile regmincoverage
	} {
		if {[info exists $option]} {
			lappend cmdline -$option [get $option]
		}
	}
	lappend cmdline {*}$opts $bamfile $refseq
	job_logfile $destdir/var_gatk_$bamtail $destdir $cmdline \
		{*}[versions bwa bowtie2 samtools gatk gatk3 picard java gnusort8 zst os]
	# start
	## Produce gatk SNP calls
	set keeppwd [pwd]
	cd $destdir
	set gatkrefseq [gatk_refseq_job $refseq]
	set deps [list $bamtail $gatkrefseq $bamtail.bai {*}$deps]
	job ${pre}varall-$root {*}$skips -mem [job_mempercore 5G $threads] -cores $threads \
	-deps $deps -targets {
		${pre}varall-$root.vcf
	} -skip {
		$varallfile
	} -vars {
		gatk opts regionfile gatkrefseq refseq threads root
	} -code {
		analysisinfo_write $dep $target sample $root varcaller gatk varcaller_version [version gatk3] varcaller_cg_version [version genomecomb] varcaller_region [filename $regionfile]
		set emptyreg [reg_isempty $regionfile]
		set cache [file dir $target]/cache_var_gatk_[file tail $refseq].temp]
		if {$emptyreg && [file exists $cache]} {
			file copy $cache $target
		} else {
			if {$regionfile ne ""} {
				set bedfile [tempbed $regionfile $refseq]
				lappend opts -L $bedfile
			}
			gatk3exec {-XX:ParallelGCThreads=1 -d64 -Xms512m -Xmx4g} UnifiedGenotyper \
				{*}$opts -nct $threads -R $dep2 -I $dep -o $target.temp \
				-stand_call_conf 10.0 -dcov 1000 \
				--annotateNDA \
				-glm SNP --output_mode EMIT_ALL_CONFIDENT_SITES 2>@ stderr >@ stdout
			file rename -force $target.temp $target
			catch {file delete $target.temp.idx}
			if {$emptyreg && ![file exists $cache]} {
				file copy $target $cache
			}
		}
	}
	job ${pre}varall-gatk2tsv-$root {*}$skips -deps {
		${pre}varall-$root.vcf
	} -targets {
		$varallfile.zst
	} -vars {
		sample split refseq
	} -code {
		analysisinfo_write $dep $target
		cg vcf2tsv -split $split -meta [list refseq [file tail $refseq]] -removefields {
			name filter AN AC AF AA ExcessHet InbreedingCoeff MLEAC MLEAF NDA RPA RU STR
		} $dep $target.temp.zst
		file rename -force $target.temp.zst $target
	}
	# zst_job $varallfile -i 1
	zstindex_job {*}$skips $varallfile.zst
	# predict deletions separately, because gatk will not predict snps in a region where a deletion
	# was predicted in the varall
	# not using threads, as these cause (sporadic) errors (https://gatkforums.broadinstitute.org/gatk/discussion/3141/unifiedgenotyper-error-somehow-the-requested-coordinate-is-not-covered-by-the-read)
	job ${pre}delvar-$root {*}$skips \
	-deps $deps \
	-targets {${pre}delvar-$root.vcf} -skip {
		${pre}delvar-$root.tsv
	} -skip {
		$varfile
	} -vars {
		gatk opts regionfile gatkrefseq refseq
	} -code {
		set emptyreg [reg_isempty $regionfile]
		set cache [file dir $target]/cache_var_gatk_[file tail $refseq].temp
		if {$emptyreg && [file exists $cache]} {
			file copy $cache $target
		} else {
			if {$regionfile ne ""} {
				set bedfile [tempbed $regionfile $refseq]
				lappend opts -L $bedfile
			}
			gatk3exec {-XX:ParallelGCThreads=1 -d64 -Xms512m -Xmx4g} UnifiedGenotyper \
				{*}$opts -R $dep2 -I $dep -o $target.temp \
				-stand_call_conf 10.0 -dcov 1000 \
				--annotateNDA \
				-glm INDEL 2>@ stderr >@ stdout
			file rename -force $target.temp $target
			catch {file delete $target.temp.idx}
			if {$emptyreg && ![file exists $cache]} {
				file copy $target $cache
			}
		}
	}
	job ${pre}delvar-gatk2tsv-$root {*}$skips -deps {
		${pre}delvar-$root.vcf
	} -targets {
		${pre}delvar-$root.tsv
	} -skip {
		$varfile
	} -vars {
		sample split refseq
	} -code {
		set temp [filetemp_ext $target]
		exec cg vcf2tsv -split $split -meta [list refseq [file tail $refseq]] -removefields {
			name filter AN AC AF AA ExcessHet InbreedingCoeff MLEAC MLEAF NDA RPA RU STR
		} $dep | cg select -overwrite 1 -q {
			$alt ne "." && $alleleSeq1 ne "." &&$quality >= 10 && $totalcoverage > 4
			&& $zyg ni "r o"
		} -f {
			chromosome begin end type ref alt quality alleleSeq1 alleleSeq2 
			{sequenced=if($quality < 30 || $totalcoverage < 5,"u","v")}
			{zyg=if($quality < 30 || $totalcoverage < 5,"u",$zyg)}
			*
		} - $temp
		file rename -force $temp $target
	}
	job ${pre}uvar-$root {*}$skips -deps {
		$varallfile.zst ${pre}delvar-$root.tsv
	} -targets {
		${pre}uvar-$root.tsv
	} -skip {
		$varfile
	} -vars {
		root pre
	} -code {
		analysisinfo_write $dep $target varcaller_mincoverage 5 varcaller_minquality 30 varcaller_cg_version [version genomecomb]
		cg select -overwrite 1 -q {$alt ne "." && $alleleSeq1 ne "." &&$quality >= 10 && $totalcoverage > 4} \
		-f {
			chromosome begin end type ref alt quality alleleSeq1 alleleSeq2 
			{sequenced=if($quality < 30 || $totalcoverage < 5,"u",if($zyg eq "r","r","v"))}
			{zyg=if($quality < 30 || $totalcoverage < 5,"u",$zyg)}
			*
		} $dep $target.temp
		cg cat -c f $target.temp $dep2 > $target.temp2
		cg select -overwrite 1 -s - $target.temp2 $target.temp3
		file rename -force $target.temp3 $target
		file delete $target.temp
		file delete $target.temp2
	}
	# annotvar_clusters_job works using jobs
	annotvar_clusters_job {*}$skips ${pre}uvar-$root.tsv $varfile.zst
	# make sreg
	sreg_gatk_job ${pre}sreg-$root $varallfile $sregfile.zst $skips
	## filter SNPs (according to seqanswers exome guide)
	# java -d64 -Xms512m -Xmx4g -jar $gatk -R $reference -T VariantFiltration -B:variant,VCF snp.vcf.recalibrated -o $outprefix.snp.filtered.vcf --clusterWindowSize 10 --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" --filterName "HARD_TO_VALIDATE" --filterExpression "DP < 5 " --filterName "LowCoverage" --filterExpression "QUAL < 30.0 " --filterName "VeryLowQual" --filterExpression "QUAL > 30.0 && QUAL < 50.0 " --filterName "LowQual" --filterExpression "QD < 1.5 " --filterName "LowQD" --filterExpression "SB > -10.0 " --filterName "StrandBias"
	# cleanup
	if {$cleanup} {
		set cleanupfiles [list \
			${pre}uvar-$root.tsv ${pre}uvar-$root.tsv.index \
			${pre}varall-$root.vcf ${pre}delvar-$root.vcf ${pre}delvar-$root.tsv \
		]
		set cleanupdeps [list $varfile $varallfile]
		cleanup_job clean_${pre}var-$root $cleanupfiles $cleanupdeps
	}
	cd $keeppwd
	return $resultlist
	# return [file join $destdir $varfile]
}

proc cg_var_gatk {args} {
	set args [job_init {*}$args]
	var_gatk_job {*}$args
	job_wait
}
