proc sreg_freebayes_job {job varallfile resultfile {skips {}}} {
	upvar job_logdir job_logdir
	job $job {*}$skips -deps {
		$varallfile
	} -targets {
		$resultfile
	} -code {
		set temp [filetemp $target]
		set temp2 [filetemp $target]
		cg select -overwrite 1 -q {$genoqual >= 30 && $totalcoverage >= 5 && $type ne "ins"} -f {chromosome begin end} $dep $temp
		file_write $temp2 "# regions selected from [gzroot $dep]: \$genoqual >= 30 && \$totalcoverage >= 5\n"
		cg regjoin $temp >> $temp2
		compress $temp2 $target 1 0
		file delete $temp
	}
}

proc var_freebayes_job {args} {
	upvar job_logdir job_logdir
	set pre ""
	set opts {}
	set split 0
	set deps {}
	set regionfile {}
	set threads 2
	set cleanup 1
	set regmincoverage 3
	set rootname {}
	set resultfiles 0
	set skips {}
	cg_options var_freebayes args {
		-L - -deps {
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
		-skip {
			lappend skips -skip $value
		}
		-opts {
			set opts $value
		}
	} {bamfile refseq} 2 2 {
		call variants using freebayes
	}
	set bamfile [file_absolute $bamfile]
	set refseq [refseq $refseq]
	set destdir [file dir $bamfile]
	set bamtail [file tail $bamfile]
	if {$rootname eq ""} {
		set root freebayes-[file_rootname $bamfile]
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
	set cmdline [list cg var_freebayes]
	foreach option {
		split deps bed pre regionfile regmincoverage
	} {
		if {[info exists $option]} {
			lappend cmdline -$option [get $option]
		}
	}
	lappend cmdline {*}$opts $bamfile $refseq
	job_logfile $destdir/var_freebayes_$bamtail $destdir $cmdline \
		{*}[versions bwa samtools freebayes picard java gnusort8 zst os]
	# start
	## Produce freebayes SNP calls
	set keeppwd [pwd]
	cd $destdir
	set deps [list $bamtail $refseq $bamtail.bai {*}$deps]
	job ${pre}varall-$root {*}$skips -mem 5G -deps $deps -targets {
		${pre}varall-$root.vcf
	} -skip {
		$varallfile
	} -vars {
		opts regionfile refseq root
	} -code {
		analysisinfo_write $dep $target sample $root varcaller freebayes varcaller_version [version freebayes] varcaller_cg_version [version genomecomb] varcaller_region [filename $regionfile]
		set emptyreg [reg_isempty $regionfile]
		set cache [file dir $target]/cache_var_freebayes_[file tail $refseq].temp
		if {$emptyreg && [file exists $cache]} {
			file copy $cache $target
		} else {
			if {$regionfile ne ""} {
				set regionfilesize [file size $regionfile]
				if {!$emptyreg} {
					set bedfile [tempbed $regionfile $refseq]
				} else {
					set bedfile [tempfile].bed
					set temp [exec samtools view -H $dep]
					regexp "\tSN:(\[^ \t\]+)\t" $temp temp chr
					file_write $bedfile $chr\t0\t1\n
				}
				lappend opts -t $bedfile
			}
			exec freebayes {*}$opts \
				--genotype-qualities --report-monomorphic --exclude-unobserved-genotypes \
				-f $refseq $dep > $target.temp 2>@ stderr
			file rename -force $target.temp $target
			catch {file delete $target.temp.idx}
			if {$emptyreg && ![file exists $cache]} {
				file copy $target $cache
			}
		}
	}
	job ${pre}varall-freebayes2tsv-$root {*}$skips -deps {
		${pre}varall-$root.vcf
	} -targets {
		$varallfile.zst
	} -vars {sample split refseq} \
	-code {
		analysisinfo_write $dep $target
		cg vcf2tsv -split $split -meta [list refseq [file tail $refseq]] -removefields {name filter AN AC AF AA ExcessHet InbreedingCoeff MLEAC MLEAF NDA RPA RU STR} $dep $target.temp.zst
		file rename -force $target.temp.zst $target
	}
	# zst_job $varallfile -i 1
	zstindex_job $varallfile.zst

	job ${pre}uvar-$root {*}$skips -deps {
		$varallfile
	} -targets {
		${pre}uvar-$root.tsv
	} -skip {
		$varfile
	} -vars {
		root pre
	} -code {
		analysisinfo_write $dep $target varcaller_mincoverage 5 varcaller_minquality 30 varcaller_cg_version [version genomecomb]
		cg select -q {$zyg ne "r" && $genoqual >= 10 && $totalcoverage > 4} \
		-f {
			chromosome begin end type ref alt quality alleleSeq1 alleleSeq2 
			{sequenced=if($genoqual < 30 || $totalcoverage < 5,"u",if($zyg eq "r","r","v"))}
			{zyg=if($genoqual < 30 || $totalcoverage < 5,"u",$zyg)}
			*
		} $dep $target.temp
		cg select -s - $target.temp $target.temp2
		file rename -force $target.temp2 $target
		file delete $target.temp
	}
	# annotvar_clusters_job works using jobs
	annotvar_clusters_job {*}$skips ${pre}uvar-$root.tsv $varfile.zst
	# make sreg
	sreg_freebayes_job ${pre}sreg-$root $varallfile $sregfile.zst $skips
	if {$cleanup} {
		set cleanupfiles [list \
			${pre}uvar-$root.tsv ${pre}uvar-$root.tsv.index \
			${pre}varall-$root.vcf \
		]
		set cleanupdeps [list $varfile $varallfile]
		cleanup_job clean_${pre}var-$root $cleanupfiles $cleanupdeps
	}
	cd $keeppwd
	return $resultlist
}

proc cg_var_freebayes {args} {
	set args [job_init {*}$args]
	set result [var_freebayes_job {*}$args]
	job_wait
	return $result
}
