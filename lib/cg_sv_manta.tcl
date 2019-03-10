proc version_manta {} {
	set version ?
	set mantadir [searchpath MANTADIR manta manta*]
	set help [catch_exec $mantadir/bin/configManta.py -h]
	regexp {Version: ([^\n]*)} $help temp version
	return $version
}

proc sv_manta_job {args} {
	upvar job_logdir job_logdir
	set cmdline "[list cd [pwd]] \; [list cg sv_manta {*}$args]"
	set refseq {}
	set opts {}
	set split 1
	set threads 2
	set cleanup 0
	set regmincoverage 3
	set resultfiles 0
	set rootname {}
	set skips {}
	set resultfile {}
	cg_options sv_manta args {
		-refseq {
			set refseq $value
		}
		-split {
			set split $value
		}
		-exome {
			set exome $value
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
		-skip {
			lappend skips -skip $value
		}
		default {
			if {[regexp {^-..} $key]} {set key -$key}
			lappend opts $key $value
		}
	} {bamfile resultfile} 1 2 {
		run SV calls using manta
	}
	set bamfile [file_absolute $bamfile]
	set refseq [refseq $refseq]
	if {$resultfile eq ""} {
		set root manta-[file_rootname $bamfile]
		set resultfile [file dir $bamfile]/sv-$root.tsv.zst
	} else {
		set root [file_rootname $resultfile]
	}
	set resultanalysisinfo [gzroot $resultfile].analysisinfo
	set destdir [file dir $resultfile]
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
	set resultlist [list $resultfile $resultanalysisinfo]
	if {$resultfiles} {
		return $resultlist
	}
	# logfile
	job_logfile $destdir/sv_manta_[file tail $resultfile] $destdir $cmdline \
		{*}[versions manta gnusort8 zst os]
	# start
	set gatkrefseq [gatk_refseq_job $refseq]
	## Produce manta sv calls
	if {[file extension $bamfile] eq ".cram"} {
 		set bamfileindex $bamfile.crai
	} else {
		set bamfileindex $bamfile.bai
	}
	job sv_manta-$root.vcf {*}$skips -mem [expr {1*$threads}]G -cores $threads \
	-skip [list $resultfile $resultfile.analysisinfo] \
	-deps {
		$bamfile $refseq $bamfileindex $refseq.fai
	} -targets {
		$resultfile.mantarun/results/variants/diploidSV.vcf.gz $resultfile.mantarun.analysisinfo
	} -vars {
		resultfile manta opts gatkrefseq threads root
	} -code {
		set mantadir [searchpath MANTADIR manta manta*]
		analysisinfo_write $dep $resultfile.mantarun sample $root varcaller manta varcaller_version [version manta] varcaller_cg_version [version genomecomb]
		exec $mantadir/bin/configManta.py {*}$opts \
			--bam $dep \
			--referenceFasta $gatkrefseq \
			--runDir $resultfile.mantarun
		catch_exec $resultfile.mantarun/runWorkflow.py -m local -j $threads
		file delete -force $resultfile.mantarun/workspace
	}
	# 
	job sv-manta-vcf2tsv-$root {*}$skips -deps {
		$resultfile.mantarun/results/variants/diploidSV.vcf.gz
	} -targets {
		$resultfile $resultanalysisinfo
	} -vars {
		sample split resultfile
	} -code {
		analysisinfo_write $resultfile.mantarun $target
		cg vcf2tsv -split $split -removefields {name filter AN AC AF AA ExcessHet InbreedingCoeff MLEAC MLEAF NDA RPA RU STR} $resultfile.mantarun/results/variants/diploidSV.vcf.gz $target.temp[gzext $target]
		file rename -force $target.temp[gzext $target] $target
		hardlink $resultfile.mantarun/results/variants/diploidSV.vcf.gz [file root [gzroot $target]].vcf.gz
		hardlink $resultfile.mantarun/results/variants/diploidSV.vcf.gz.tbi [file root [gzroot $target]].vcf.gz.tbi
	}
	# cleanup
	return $resultlist
}

proc cg_sv_manta {args} {
	set args [job_init {*}$args]
	sv_manta_job {*}$args
	job_wait
}

