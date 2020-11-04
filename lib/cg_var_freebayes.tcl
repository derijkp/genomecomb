proc var_sam_tools {} {
	return {freebayes}
}

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
	set resultfile {}
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
		-datatype {
			# not actually used
		}
		-skip {
			lappend skips -skip $value
		}
		-opts {
			set opts $value
		}
	} {bamfile refseq resultfile} 2 3 {
		call variants using freebayes
	}
	set bamfile [file_absolute $bamfile]
	set refseq [refseq $refseq]
	if {$resultfile eq ""} {
		set resultfile [file dir $bamfile]/${pre}var-freebayes-[file_rootname $bamfile].tsv.zst
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
	set resultlist [list $varfile $sregfile $varallfile $regclusterfile]
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
	job_logfile $destdir/var_freebayes_$resulttail $destdir $cmdline \
		{*}[versions bwa samtools freebayes picard java gnusort8 zst os]
	# start
	## Produce freebayes SNP calls
	set bamindex $bamfile.[indexext $bamfile]
	set deps [list $bamfile $refseq $bamindex {*}$deps]
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
					set temp [exec samtools view --no-PG -H $dep]
					regexp "\tSN:(\[^ \t\]+)\t" $temp temp chr
					file_write $bedfile $chr\t0\t1\n
				}
				lappend opts -t $bedfile
			}
			exec freebayes {*}$opts \
				--genotype-qualities --report-monomorphic --exclude-unobserved-genotypes \
				-f $refseq $dep > $target.temp 2>@ stderr
			file rename -force -- $target.temp $target
			catch {file delete $target.temp.idx}
			if {$emptyreg && ![file exists $cache]} {
				file copy $target $cache
			}
		}
	}
	job ${pre}varall-freebayes2tsv-$root {*}$skips -deps {
		${pre}varall-$root.vcf
	} -targets {
		$varallfile
	} -vars {
		sample split refseq
	} -code {
		analysisinfo_write $dep $target
		cg vcf2tsv -split $split -meta [list refseq [file tail $refseq]] -removefields {name filter AN AC AF AA ExcessHet InbreedingCoeff MLEAC MLEAF NDA RPA RU STR} $dep $target.temp.zst
		file rename -force -- $target.temp.zst $target
	}
	# zst_job $varallfile -i 1
	zstindex_job $varallfile
	job ${pre}uvar-$root {*}$skips -deps {
		$varallfile
	} -targets {
		$uvarfile
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
		cg select -overwrite 1 -s - $target.temp $target.temp2[gzext $target]
		result_rename $target.temp2[gzext $target] $target
		file delete $target.temp
	}
	# annotvar_clusters_job works using jobs
	annotvar_clusters_job {*}$skips $uvarfile $varfile
	# make sreg
	sreg_freebayes_job ${pre}sreg-$root $varallfile $sregfile $skips
	if {$cleanup} {
		set cleanupfiles [list \
			$uvarfile [gzroot $uvarfile].index [gzroot $uvarfile].temp \
			${pre}varall-$root.vcf \
		]
		set cleanupdeps [list $varfile $varallfile]
		cleanup_job clean_${pre}var-$root $cleanupfiles $cleanupdeps
	}
	return $resultlist
}

proc cg_var_freebayes {args} {
	set args [job_init {*}$args]
	set result [var_freebayes_job {*}$args]
	job_wait
	return $result
}
