proc sreg_freebayes_job {job varallfile resultfile} {
	upvar job_logdir job_logdir
	job $job -deps {$varallfile} -targets {$resultfile} -code {
		set temp [filetemp $target]
		set temp2 [filetemp $target]
		cg select -overwrite 1 -q {$genoqual >= 30 && $totalcoverage >= 5 && $type ne "ins"} -f {chromosome begin end} $dep $temp
		file_write $temp2 "# regions selected from [gzroot $dep]: \$genoqual >= 30 && \$totalcoverage >= 5\n"
		cg regjoin $temp >> $temp2
		if {[file extension $target] eq ".lz4"} {
			cg_lz4 -keep 0 -i 1 -o $target $temp2
		} else {
			file rename $temp2 $target
		}
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
	cg_options var_sam args {
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
		default {
			lappend opts $key $value
		}
	} {bamfile refseq} 2 2 {
		call variants using freebayes
	}
	set bamfile [file_absolute $bamfile]
	set refseq [file_absolute $refseq]
	set destdir [file dir $bamfile]
	if {$regionfile ne ""} {
		set regionfile [file_absolute $regionfile]
	} else {
		set regionfile [bam2reg_job -mincoverage $regmincoverage $bamfile]
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
	job_logfile $destdir/var_freebayes_[file tail $bamfile] $destdir $cmdline \
		{*}[versions bwa samtools freebayes picard java gnusort8 lz4 os]
	set file [file tail $bamfile]
	set root [file_rootname $file]
	if {$root eq ""} {set root [file root $file]}
	# start
	## Produce freebayes SNP calls
	set keeppwd [pwd]
	cd $destdir
	set deps [list $file $refseq $file.bai {*}$deps]
	job ${pre}varall-freebayes-$root -mem 5G -cores $threads -deps $deps -targets {
		${pre}varall-freebayes-$root.vcf
		${pre}varall-freebayes-$root.vcf.analysisinfo
	} -skip [list ${pre}varall-freebayes-$root.tsv ${pre}varall-freebayes-$root.tsv.analysisinfo] -vars {
		opts regionfile refseq threads root
	} -code {
		analysisinfo_write $dep $target sample freebayes-$root varcaller freebayes varcaller_version [version freebayes] varcaller_cg_version [version genomecomb] varcaller_region [filename $regionfile]
		if {$regionfile ne ""} {
			set bedfile [tempbed $regionfile $refseq]
			lappend opts -t $bedfile
		}
		exec freebayes {*}$opts \
			--genotype-qualities --report-monomorphic --exclude-unobserved-genotypes \
			-f $refseq $dep > $target.temp 2>@ stderr
		file rename -force $target.temp $target
	}
	job ${pre}varall-freebayes2tsv-$root -deps {
		${pre}varall-freebayes-$root.vcf
	} -targets {
		${pre}varall-freebayes-$root.tsv.lz4
		${pre}varall-freebayes-$root.tsv.analysisinfo
	} -vars {sample split} \
	-code {
		analysisinfo_write $dep $target
		cg vcf2tsv -split $split -removefields {name filter AN AC AF AA ExcessHet InbreedingCoeff MLEAC MLEAF NDA RPA RU STR} $dep $target.temp.lz4
		file rename -force $target.temp.lz4 $target
	}
	# lz4_job ${pre}varall-freebayes-$root.tsv -i 1
	lz4index_job ${pre}varall-freebayes-$root.tsv.lz4

	job ${pre}uvar-freebayes-$root -deps {
		${pre}varall-freebayes-$root.tsv
	} -targets {
		${pre}uvar-freebayes-$root.tsv
		${pre}uvar-freebayes-$root.tsv.analysisinfo
	} -skip [list ${pre}var-freebayes-$root.tsv ${pre}var-freebayes-$root.tsv.analysisinfo] -vars {
		root pre
	} -code {
		analysisinfo_write $dep $target varcaller_mincoverage 5 varcaller_minquality 30 varcaller_cg_version [version genomecomb]
		cg select -q {$alt ne "." && $alleleSeq1 ne "." &&$genoqual >= 10 && $totalcoverage > 4} \
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
	annotvar_clusters_job ${pre}uvar-freebayes-$root.tsv ${pre}var-freebayes-$root.tsv.lz4
	# make sreg
	sreg_freebayes_job ${pre}sreg-freebayes-$root ${pre}varall-freebayes-$root.tsv ${pre}sreg-freebayes-$root.tsv.lz4
	if {$cleanup} {
		set cleanupfiles [list \
			${pre}uvar-freebayes-$root.tsv ${pre}uvar-freebayes-$root.tsv.index ${pre}uvar-freebayes-$root.tsv.analysisinfo \
			${pre}varall-freebayes-$root.vcf ${pre}varall-freebayes-$root.vcf.analysisinfo \
		]
		set cleanupdeps [list ${pre}var-freebayes-$root.tsv ${pre}varall-freebayes-$root.tsv]
		cleanup_job clean_${pre}var-freebayes-$root $cleanupfiles $cleanupdeps
	}
	cd $keeppwd
	return [file join $destdir ${pre}var-freebayes-$root.tsv]
}

proc cg_var_freebayes {args} {
	set args [job_init {*}$args]
	set result [var_freebayes_job {*}$args]
	job_wait
	return $result
}