proc version_strelka {} {
	set version ?
	set strelkadir [searchpath STRELKAADIR strelka strelka*]
	set version [catch_exec $strelkadir/bin/configureStrelkaGermlineWorkflow.py --version]
	return $version
}

proc sreg_strelka_job {job varallfile resultfile {mincoverage 8} {mingenoqual 25} {skips {}}} {
	upvar job_logdir job_logdir
	job $job {*}$skips -deps {$varallfile} -targets {$resultfile} -vars {mincoverage mingenoqual} -code {
		set temp [filetemp $target]
		if {![file size $dep]} {
			file_write $target ""
		} else {
			exec cg vcf2tsv $dep \
				| cg select -q [subst {
					(def(\$genoqual,0) >= $mingenoqual || def(\$GQX,0) >= $mingenoqual) && def(\$coverage,0) >= $mincoverage && \$type ne "ins"
				}] -f {chromosome begin end} \
				| cg regjoin {*}[compresspipe $target] > $temp
			file rename -force $temp $target
			if {[file extension $target] eq ".lz4"} {
				exec lz4index $target
			}
		}
	}
}

proc var_strelka_job {args} {
	upvar job_logdir job_logdir
	set cmdline "[list cd [pwd]] \; [list cg sv_manta {*}$args]"
	set pre ""
	set opts {}
	set split 1
	set deps {}
	set regionfile {}
	set threads 2
	set cleanup 1
	set regmincoverage 3
	set mincoverage 8
	set mingenoqual 25
	set resultfiles 0
	set rootname {}
	set skips {}
	cg_options var_strelka args {
		-L - -deps {
			lappend deps $value
		}
		-exome {
			set exome $value
		}
		-regionfile {
			set regionfile $value
		}
		-regmincoverage {
			set regmincoverage $value
		}
		-mincoverage {
			set mincoverage $value
		}
		-mingenoqual {
			set mingenoqual $value
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
		default {
			lappend opts $key $value
		}
	} {bamfile refseq}
	set bamfile [file_absolute $bamfile]
	set refseq [refseq $refseq]
	set destdir [file dir $bamfile]
	set file [file tail $bamfile]
	if {$rootname eq ""} {
		set root strelka-[file_rootname $file]
	} else {
		set root $rootname
	}
	if {![info exists exome]} {
		if {[file exists [gzfile $destdir/reg_targets-*.tsv]]} {
			set exome 1
		} else {
			set exome 0
		}
	}
	if {$exome} {
		lappend opts --exome
	}
	# resultfiles
	set varfile ${pre}var-$root.tsv
	set sregfile ${pre}sreg-$root.tsv
	set varallfile ${pre}varall-$root.gvcf
	set resultlist [list $destdir/$varfile.lz4 $destdir/$sregfile.lz4 $destdir/$varallfile.gz $destdir/reg_cluster-$root.tsv.lz4]
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
	job_logfile $destdir/var_strelka_[file tail $bamfile] $destdir $cmdline \
		{*}[versions strelka java gnusort8 lz4 os]
	# start
	## Produce strelka SNP calls
	set keeppwd [pwd]
	cd $destdir
	set resultgvcf $varallfile
	set resultname $varallfile
	set resultvcf [file root [gzroot $varfile]].vcf
	job ${pre}varall-$root {*}$skips -mem [job_mempercore 5G $threads] -cores $threads -skip [list \
		$varallfile $varallfile.analysisinfo \
	] -deps [list \
		$bamfile $refseq $bamfile.bai {*}$deps \
	] -targets {
		$resultgvcf.gz $resultgvcf.gz.tbi $resultgvcf.analysisinfo
		$resultvcf.gz $resultvcf.gz.tbi $resultvcf.analysisinfo
	} -vars {
		bamfile resultgvcf resultvcf opts regionfile refseq threads root refseq varfile
	} -code {
		analysisinfo_write $bamfile $resultgvcf.gz sample $root varcaller strelka varcaller_version [version strelka] varcaller_cg_version [version genomecomb] varcaller_region [filename $regionfile]
		analysisinfo_write $bamfile $resultvcf.gz sample $root varcaller strelka varcaller_version [version strelka] varcaller_cg_version [version genomecomb] varcaller_region [filename $regionfile]
		set zerosize 0
		if {$regionfile ne ""} {
			set regionfile [gzfile $regionfile]
			set bedfile [tempbed $regionfile $refseq]
			if {[file size $bedfile] == 0} {
				set zerosize 1
			} else {
				if {[file extension $bedfile] ne ".gz"} {
					exec bgzip $bedfile
					set bedfile $bedfile.gz
				}
				exec tabix -f -p bed $bedfile
				lappend opts --callRegions=$bedfile
			}
		}
		if {!$zerosize} {
			file delete -force $varfile.strelkarun
			file mkdir $varfile.strelkarun
			set strelkadir [searchpath STRELKAADIR strelka strelka*]
			exec $strelkadir/bin/configureStrelkaGermlineWorkflow.py {*}$opts \
				--bam $dep \
				--referenceFasta $refseq \
				--runDir $varfile.strelkarun
			catch_exec $varfile.strelkarun/runWorkflow.py -m local -j $threads
			file rename -force $varfile.strelkarun/results/variants/genome.S1.vcf.gz $resultgvcf.gz
			file rename -force $varfile.strelkarun/results/variants/genome.S1.vcf.gz.tbi $resultgvcf.gz.tbi
			file rename -force $varfile.strelkarun/results/variants/variants.vcf.gz $resultvcf.gz
			file rename -force $varfile.strelkarun/results/variants/variants.vcf.gz.tbi $resultvcf.gz.tbi
			file delete -force $varfile.strelkarun
		} else {
			putslog "empty regionfile -> write empty $resultgvcf.gz"
			file_write $resultgvcf.gz ""
			file_write $resultgvcf.gz.tbi ""
			file_write $resultgvcf.analysisinfo ""
			file_write $resultvcf.gz ""
			file_write $resultvcf.gz.tbi ""
			file_write $resultvcf.analysisinfo ""
		}
	}
	job ${pre}vcf2tsv-$root {*}$skips -deps {
		$resultvcf.gz $resultvcf.analysisinfo
	} -targets {
		${pre}uvar-$root.tsv ${pre}uvar-$root.tsv.analysisinfo
	} -vars {
		pre root resultvcf sample split mincoverage mingenoqual type refseq
	} -skip {
		$varfile.lz4 $varfile.analysisinfo
	} -code {
		analysisinfo_write $dep $target varcaller_mincoverage $mincoverage varcaller_mingenoqual $mingenoqual varcaller_cg_version [version genomecomb]
		set fields {chromosome begin end type ref alt quality alleleSeq1 alleleSeq2}
		lappend fields [subst {sequenced=if(def(\$genoqual,0) < $mingenoqual || def(\$coverage,0) < $mincoverage,"u","v")}]
		lappend fields [subst {zyg=if(def(\$genoqual,0) < $mingenoqual || def(\$coverage,0) < $mincoverage,"u",\$zyg)}]
		lappend fields *
		if {[file size $resultvcf.gz]} {
			exec cg vcf2tsv -split $split -meta [list refseq [file tail $refseq]] -removefields {
				name filter AN AC AF AA ExcessHet InbreedingCoeff MLEAC MLEAF NDA RPA RU STR
			} $resultvcf.gz | cg select -f $fields > ${pre}uvar-$root.tsv.temp
		} else {
			file_write ${pre}uvar-$root.tsv.temp ""
		}
		file rename -force ${pre}uvar-$root.tsv.temp ${pre}uvar-$root.tsv
	}
	# annotvar_clusters_job works using jobs
	annotvar_clusters_job {*}$skips ${pre}uvar-$root.tsv $varfile.lz4
	# make sreg
	sreg_strelka_job ${pre}sreg-$root $varallfile $sregfile.lz4 $mincoverage $mingenoqual $skips
	# cleanup
	if {$cleanup} {
		set cleanupfiles [list \
			${pre}var-$root.vcf \
			${pre}uvar-$root.tsv ${pre}uvar-$root.tsv.index
		]
		set cleanupdeps [list $varfile $resultgvcf.gz]
		cleanup_job clean_${pre}var-$root $cleanupfiles $cleanupdeps
	}
	cd $keeppwd
	return $resultlist
	# return [file join $destdir $varfile]
}

proc cg_var_strelka {args} {
	set args [job_init {*}$args]
	var_strelka_job {*}$args
	job_wait
}
