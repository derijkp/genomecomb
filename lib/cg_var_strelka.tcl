proc var_strelka_tools {} {
	return {strelka}
}

proc version_strelka {} {
	set version ?
	set strelkadir [searchpath STRELKAADIR strelka strelka*]
	set version [catch_exec $strelkadir/bin/configureStrelkaGermlineWorkflow.py --version]
	return $version
}

proc sreg_strelka_job {job varallfile resultfile {mincoverage 8} {mingenoqual 25} {skips {}}} {
	upvar job_logdir job_logdir
	job $job {*}$skips -deps {
		$varallfile
	} -targets {
		$resultfile 
	} -vars {
		mincoverage mingenoqual
	} -code {
		set temp [filetemp $target]
		if {![file size $dep]} {
			file_write $target ""
		} else {
			exec cg vcf2tsv $dep \
				| cg select -q [subst {
					(def(\$genoqual,0) >= $mingenoqual || def(\$GQX,0) >= $mingenoqual) && def(\$coverage,0) >= $mincoverage && \$type ne "ins"
				}] -f {chromosome begin end} \
				| cg regjoin {*}[compresspipe $target] > $temp
			file rename -force -- $temp $target
			cg_zindex $target
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
	set datatype {}
	set skips {}
	cg_options var_strelka args {
		-L - -deps {
			lappend deps $value
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
		-datatype {
			set datatype $value
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
		set root strelka-[file_rootname $bamtail]
	} else {
		set root $rootname
	}
	if {$datatype eq ""} {
		if {[file exists [ampliconsfile $destdir]]} {
			set datatype targeted
		} elseif {[file exists [targetfile $destdir]]} {
			set datatype exome
		} else {
			set datatype genome
		}
	}
	if {$datatype eq "exome"} {
		lappend opts --exome
	} elseif {$datatype eq "amplicons"} {
		lappend opts --targeted
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
	job_logfile $destdir/var_strelka_[file tail $bamfile] $destdir $cmdline \
		{*}[versions strelka java gnusort8 zst os]
	# start
	## Produce strelka SNP calls
	set keeppwd [pwd]
	cd $destdir
	set resultgvcf $varallfile
	set resultname $varallfile
	set resultvcf [file root [gzroot $varfile]].vcf
	set bamfileindex $bamfile.[indexext $bamfile]
	job ${pre}varall-$root {*}$skips -mem [job_mempercore 5G $threads] -cores $threads -skip {
		$varallfile $varfile $varfile.analysisinfo
	} -deps [list \
		$bamfile $refseq $bamfileindex {*}$deps \
	] -targets {
		$resultgvcf.gz $resultgvcf.gz.tbi
		$resultvcf.gz $resultvcf.gz.tbi
	} -vars {
		bamfile resultgvcf resultvcf opts regionfile refseq threads root refseq varfile
	} -code {
		analysisinfo_write $bamfile $resultgvcf.gz sample $root varcaller strelka varcaller_version [version strelka] varcaller_cg_version [version genomecomb] varcaller_region [filename $regionfile]
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
			set strelkarun [scratchdir]/strelkarun
			file delete -force $strelkarun
			file mkdir $strelkarun
			set strelkadir [searchpath STRELKAADIR strelka strelka*]
			exec $strelkadir/bin/configureStrelkaGermlineWorkflow.py {*}$opts \
				--bam $dep \
				--referenceFasta $refseq \
				--runDir $strelkarun
			if {$threads == 1} {
				catch_exec $strelkarun/runWorkflow.py -m local -j 1
			} else {
				catch_exec $strelkarun/runWorkflow.py -m local -j $threads
			}
			file rename -force -- $strelkarun/results/variants/genome.S1.vcf.gz $resultgvcf.gz
			file rename -force -- $strelkarun/results/variants/genome.S1.vcf.gz.tbi $resultgvcf.gz.tbi
			file rename -force -- $strelkarun/results/variants/variants.vcf.gz $resultvcf.gz
			file rename -force -- $strelkarun/results/variants/variants.vcf.gz.tbi $resultvcf.gz.tbi
			file delete -force $strelkarun
		} else {
			putslog "empty regionfile -> write empty $resultgvcf.gz"
			file_write $resultgvcf.gz ""
			file_write $resultgvcf.gz.tbi ""
			file_write $resultgvcf.analysisinfo ""
			file_write $resultvcf.gz ""
			file_write $resultvcf.gz.tbi ""
			file_write $resultvcf.analysisinfo ""
		}
		analysisinfo_write $bamfile $resultvcf.gz sample $root varcaller strelka varcaller_version [version strelka] varcaller_cg_version [version genomecomb] varcaller_region [filename $regionfile]
	}
	job ${pre}vcf2tsv-$root {*}$skips -deps {
		$resultvcf.gz
	} -targets {
		${pre}uvar-$root.tsv
	} -vars {
		pre root resultvcf sample split mincoverage mingenoqual type refseq
	} -skip {
		$varfile
	} -code {
		analysisinfo_write $dep $target varcaller_mincoverage $mincoverage varcaller_mingenoqual $mingenoqual varcaller_cg_version [version genomecomb]
		set fields {chromosome begin end type ref alt quality alleleSeq1 alleleSeq2}
		lappend fields [subst {sequenced=if(def(\$genoqual,0) < $mingenoqual || def(\$coverage,0) < $mincoverage,"u","v")}]
		lappend fields [subst {zyg=if(def(\$genoqual,0) < $mingenoqual || def(\$coverage,0) < $mincoverage,"u",\$zyg)}]
		lappend fields *
		if {[file size $resultvcf.gz]} {
			exec cg vcf2tsv -split $split -meta [list refseq [file tail $refseq]] -removefields {
				name AN AC AF AA ExcessHet InbreedingCoeff MLEAC MLEAF NDA RPA RU STR
			} $resultvcf.gz | cg select -f $fields > ${pre}uvar-$root.tsv.temp
		} else {
			file_write ${pre}uvar-$root.tsv.temp ""
		}
		file rename -force -- ${pre}uvar-$root.tsv.temp ${pre}uvar-$root.tsv
	}
	# annotvar_clusters_job works using jobs
	annotvar_clusters_job {*}$skips ${pre}uvar-$root.tsv $varfile.zst
	# make sreg
	sreg_strelka_job ${pre}sreg-$root $varallfile $sregfile.zst $mincoverage $mingenoqual $skips
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
