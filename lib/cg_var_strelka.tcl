proc strelka_env {} {
	set strelkadir [findstrelka]
	set ::env(PATH) $strelkadir/bin:$::env(PATH)
	set ::env(LD_LIBRARY_PATH) $strelkadir/lib:$::env(LD_LIBRARY_PATH)
}

proc findstrelka {} {
	global strelkadir
	if {![info exists strelkadir]} {
		if {![catch {exec which configureStrelkaGermlineWorkflow.py} temp]} {
			set temp [file_resolve $temp]
			set strelkadir [file dir $temp]
		} else {
			set strelkadir [searchpath STRELKADIR strelka strelka*]
			set strelkadir [file_resolve $strelkadir]
			if {![file isdir $strelkadir]} {set strelkadir [file dir $strelkadir]}
		}
	}
	return $strelkadir
}

proc var_strelka_tools {} {
	return {strelka}
}

proc version_strelka {} {
	set version ?
	strelka_env
	set version [catch_exec configureStrelkaGermlineWorkflow.py --version]
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
	set cmdline [clean_cmdline cg var_strelka {*}$args]
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
	set resultfile {}
	set mem 5G
	set time 3:00:00
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
		set resultfile [file dir $bamfile]/${pre}var-strelka-[file_rootname $bamfile].tsv.zst
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
	set varfile $resultfile
	set vcffile [file root [gzroot $varfile]].vcf.gz
	set varallfile $destdir/${pre}varall-$root.gvcf.gz
	set uvarfile $destdir/${pre}uvar-$root.tsv.zst
	set sregfile $destdir/${pre}sreg-$root.tsv.zst
	set regclusterfile $destdir/reg_cluster-$root.tsv.zst
	set resultlist [list $varfile $sregfile $varallfile $vcffile $regclusterfile]
	if {$resultfiles} {
		return $resultlist
	}
	lappend skips -skip $resultlist
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
	set resultgvcf $varallfile
	set resultname $varallfile
	set resultvcf $vcffile
	set bamfileindex $bamfile.[indexext $bamfile]
	job ${pre}varall-$root {*}$skips -mem $mem -time $time -cores $threads -skip {
		$varallfile $varfile $varfile.analysisinfo
	} -deps [list \
		$bamfile $refseq $bamfileindex {*}$deps \
	] -targets {
		$resultgvcf $resultgvcf.tbi
		$resultvcf $resultvcf.tbi
	} -vars {
		bamfile resultgvcf resultvcf opts regionfile refseq threads root refseq varfile
	} -code {
		analysisinfo_write $bamfile $resultgvcf analysis $root sample $root varcaller strelka varcaller_version [version strelka] varcaller_cg_version [version genomecomb] varcaller_region [filename $regionfile]
		set zerosize 0
		if {$regionfile ne ""} {
			set regionfile [gzfile $regionfile]
			set bedfile [tempbed $regionfile $refseq]
			if {[file size $bedfile] == 0} {
				set zerosize 1
			} else {
				if {[file extension $bedfile] ne ".gz"} {
					if {[file exists $bedfile.gz]} {file delete $bedfile.gz}
					exec bgzip $bedfile
					set bedfile $bedfile.gz
				}
				exec tabix -f -p bed $bedfile
				lappend opts --callRegions=$bedfile
			}
		}
		if {!$zerosize} {
			strelka_env
			set strelkarun [scratchdir]/strelkarun
			file delete -force $strelkarun
			file mkdir $strelkarun
			exec configureStrelkaGermlineWorkflow.py {*}$opts \
				--bam $dep \
				--referenceFasta $refseq \
				--runDir $strelkarun
			if {$threads == 1} {
				catch_exec $strelkarun/runWorkflow.py -m local -j 1
			} else {
				catch_exec $strelkarun/runWorkflow.py -m local -j $threads
			}
			file rename -force -- $strelkarun/results/variants/genome.S1.vcf.gz $resultgvcf
			file rename -force -- $strelkarun/results/variants/genome.S1.vcf.gz.tbi $resultgvcf.tbi
			file rename -force -- $strelkarun/results/variants/variants.vcf.gz $resultvcf
			file rename -force -- $strelkarun/results/variants/variants.vcf.gz.tbi $resultvcf.tbi
			file delete -force $strelkarun
		} else {
			putslog "empty regionfile -> write empty $resultgvcf"
			file_write $resultgvcf ""
			file_write $resultgvcf.tbi ""
			file_write [analysisinfo_file $resultgvcf] ""
			file_write $resultvcf ""
			file_write $resultvcf.tbi ""
			file_write [analysisinfo_file $resultvcf] ""
		}
		analysisinfo_write $bamfile $resultvcf analysis $root sample $root varcaller strelka varcaller_version [version strelka] varcaller_cg_version [version genomecomb] varcaller_region [filename $regionfile]
	}
	job ${pre}vcf2tsv-$root {*}$skips -deps {
		$resultvcf
	} -targets {
		$uvarfile
	} -vars {
		pre root resultvcf sample split mincoverage mingenoqual type refseq uvarfile
	} -skip {
		$varfile
	} -code {
		analysisinfo_write $dep $target varcaller_mincoverage $mincoverage varcaller_mingenoqual $mingenoqual varcaller_cg_version [version genomecomb]
		set fields {chromosome begin end type ref alt quality alleleSeq1 alleleSeq2}
		lappend fields [subst {sequenced=if(def(\$genoqual,0) < $mingenoqual || def(\$coverage,0) < $mincoverage,"u","v")}]
		lappend fields [subst {zyg=if(def(\$genoqual,0) < $mingenoqual || def(\$coverage,0) < $mincoverage,"u",\$zyg)}]
		lappend fields *
		if {[file size $resultvcf]} {
			exec cg vcf2tsv -split $split -meta [list refseq [file tail $refseq]] -removefields {
				name AN AC AF AA ExcessHet InbreedingCoeff MLEAC MLEAF NDA RPA RU STR
			} $resultvcf | cg select -f $fields | cg zst > $uvarfile.temp
		} else {
			file_write $uvarfile.temp ""
		}
		file rename -force -- $uvarfile.temp $uvarfile
	}
	# annotvar_clusters_job works using jobs
	annotvar_clusters_job {*}$skips $uvarfile $varfile
	# make sreg
	sreg_strelka_job ${pre}sreg-$root $varallfile $sregfile $mincoverage $mingenoqual $skips
	# cleanup
	if {$cleanup} {
		set cleanupfiles [list \
			$uvarfile [gzroot $uvarfile].index [gzroot $uvarfile].temp \
		]
		set cleanupdeps [list $varfile $resultgvcf]
		cleanup_job clean_${pre}var-$root $cleanupfiles $cleanupdeps
	}
	return $resultlist
	# return [file join $destdir $varfile]
}

proc cg_var_strelka {args} {
	set args [job_init {*}$args]
	var_strelka_job {*}$args
	job_wait
}
