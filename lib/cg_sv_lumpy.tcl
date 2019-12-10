proc version_lumpy {} {
	set version ?
	catch {
		exec lumpy
	} c
	regexp {\(v ([^)]*)\)} $c temp version
	return $version
}

proc sv_lumpy_job {args} {
	upvar job_logdir job_logdir
	set cmdline "[list cd [pwd]] \; [list cg sv_lumpy {*}$args]"
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
	cg_options sv_lumpy args {
		-refseq {
			set refseq $value
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
		-exome {
			# notused
		}
		-skip {
			lappend skips -skip $value
		}
		default {
			if {[regexp {^-..} $key]} {set key -$key}
			lappend opts $key $value
		}
	} {bamfile resultfile} 1 2
	set bamfile [file_absolute $bamfile]
	set refseq [refseq $refseq]
	if {$resultfile eq ""} {
		set root lumpy-[file_rootname $bamfile]
		set resultfile [file dir $bamfile]/sv-$root.tsv.zst
	} else {
		set root [file_rootname $resultfile]
	}
	set vcffile [file root [gzroot $resultfile]].vcf
	set resultanalysisinfo [gzroot $resultfile].analysisinfo
	set destdir [file dir $resultfile]
	# check dependencies (and adapt config if needed)
	
	# resultfiles
	set resultlist [list $resultfile $resultanalysisinfo]
	if {$resultfiles} {
		return $resultlist
	}
	# logfile
	job_logfile $destdir/sv_lumpy_[file tail $resultfile] $destdir $cmdline \
		{*}[versions lumpy gnusort8 zst os]
	# start
	set gatkrefseq [gatk_refseq_job $refseq]
	## Produce lumpy sv calls
	set bamfileindex $bamfile.[indexext $bamfile]
	job sv_lumpy-$root.vcf {*}$skips -mem [expr {1*$threads}]G -cores $threads \
	-skip [list $resultfile $resultfile.analysisinfo] \
	-deps {
		$bamfile $refseq $bamfileindex $refseq.fai
	} -targets {
		$vcffile $vcffile.analysisinfo
	} -vars {
		resultfile vcffile lumpy opts refseq threads root
	} -code {
		analysisinfo_write $dep $vcffile sample $root varcaller lumpy varcaller_version [version lumpy] varcaller_cg_version [version genomecomb]
		catch_exec lumpyexpress {*}$opts \
			-B $dep \
			-T [scratchdir] \
			-R $refseq \
			-o $vcffile.temp.vcf 2>@ stderr >@ stdout
		file rename -force -- $vcffile.temp.vcf $vcffile
	}
	# 
	job sv-lumpy-vcf2tsv-$root {*}$skips -deps {
		$vcffile
	} -targets {
		$resultfile $resultanalysisinfo
	} -vars {
		vcffile sample split resultfile
	} -code {
		analysisinfo_write $dep $target
		cg vcf2tsv -split $split -removefields {
			name filter AN AC AF AA ExcessHet InbreedingCoeff MLEAC MLEAF NDA RPA RU STR
		} $vcffile $target.temp[gzext $target]
		file rename -force -- $target.temp[gzext $target] $target
	}
	# cleanup
	return $resultlist
}

proc cg_sv_lumpy {args} {
	set args [job_init {*}$args]
	sv_lumpy_job {*}$args
	job_wait
}

