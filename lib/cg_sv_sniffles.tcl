proc version_sniffles {} {
	set version ?
	catch {exec sniffles} msg
	if {![regexp {Version: *([^\n]+)} $msg temp version]} {
		set exe [exec which sniffles]
		if {![catch {file link $exe} exe]} {
			regsub {^sniffles-} [file tail $exe] {} version
		}
	}
	return $version
}

proc sv_sniffles_job {args} {
	upvar job_logdir job_logdir
	set cmdline "[list cd [pwd]] \; [list cg sv_sniffles {*}$args]"
	set refseq {}
	set opts {}
	set split 1
	set threads 2
	set cleanup 1
	set regmincoverage 3
	set resultfiles 0
	set skips {}
	set min_support 2
	set min_seq_size 300
	set resultfile {}
	cg_options sv_sniffles args {
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
		-skip {
			lappend skips -skip $value
		}
		-maxdist {
			lappend opts -d $value
		}
		-min_support {
			set min_support $value
		}
		-min_seq_size {
			set min_seq_size $value
		}
		-snifflesopts {
			lappend opts {*}$value
		}
	} {bamfile resultfile} 1 2
	set bamfile [file_absolute $bamfile]
	set refseq [refseq $refseq]
	if {$resultfile eq ""} {
		set root sniffles-[file_rootname $bamfile]
		set resultfile [file dir $bamfile]/sv-$root.tsv.zst
	} else {
		set root [file_rootname $resultfile]
	}
	set destdir [file dir $resultfile]
	set resultanalysisinfo [gzroot $resultfile].analysisinfo
	set vcffile [file root [gzroot $resultfile]].vcf
	# resultfiles
	set resultlist [list $resultfile $resultanalysisinfo $vcffile]
	if {$resultfiles} {
		return $resultlist
	}
	# logfile
	job_logfile $destdir/sv_sniffles_[file tail $resultfile] $destdir $cmdline \
		{*}[versions sniffles gnusort8 zst os]
	# start
	## Produce sniffles SNP calls
	set keeppwd [pwd]
	cd $destdir
	set bamfileindex $bamfile.[indexext $bamfile]
	job sv_$root.vcf {*}$skips -mem [job_mempercore 1G $threads] -cores $threads \
	-skip [list $resultfile $resultanalysisinfo] \
	-deps {
		$bamfile $refseq $bamfileindex
	} -targets {
		$vcffile $vcffile.analysisinfo
	} -vars {
		sniffles opts refseq threads root min_support min_seq_size
	} -code {
		analysisinfo_write $dep $target sample $root varcaller sniffles varcaller_version [version sniffles] varcaller_cg_version [version genomecomb]
		exec sniffles {*}$opts --threads $threads --genotype --skip_parameter_estimation \
			--min_support $min_support --min_seq_size $min_seq_size \
			-m $dep -v $target.temp 2>@ stderr >@ stdout
		file rename -force -- $target.temp $target
	}
	# 
	job sv_vcf2tsv-$root {*}$skips -deps {
		$vcffile
	} -targets {
		$resultfile $resultanalysisinfo
	} -vars {
		sample split
	} -code {
		analysisinfo_write $dep $target
		cg vcf2tsv -locerror correct -split $split -removefields {name filter AN AC AF AA ExcessHet InbreedingCoeff MLEAC MLEAF NDA RPA RU STR} $dep $target.temp[gzext $target]
		file rename -force -- $target.temp[gzext $target] $target
	}
	# cleanup
	cd $keeppwd
	return $resultlist
}

proc cg_sv_sniffles {args} {
	set args [job_init {*}$args]
	sv_sniffles_job {*}$args
	job_wait
}
