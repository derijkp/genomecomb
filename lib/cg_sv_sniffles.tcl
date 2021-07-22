proc sv_sniffles_sortdistrreg {} {
	return 1
}

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
	set snifflesopts {}
	set resultfile {}
	set region {}
	set sample {}
	cg_options sv_sniffles args {
		-refseq {
			set refseq $value
		}
		-split {
			set split $value
		}
		-threads - -t {
			if {$value > 2} {
				puts stderr "warning: reduced -threads to 2 for running sniffles (use distrreg for more efficient distribution, you can overrule this with -sniffles-threads)"
				set value 2
			}
			set threads $value
		}
		-cleanup {
			set cleanup $value
		}
		-resultfiles {
			set resultfiles $value
		}
		-preset {
			# not used
		}		
		-skip {
			lappend skips -skip $value
		}
		-maxdist {
			lappend opts -d $value
		}
		-region {
			set region $value
		}
		-sample {
			set sample $value
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
	foreach {key value} [specialopts -sniffles] {
		switch $key {
			-regmincoverage {set regmincoverage $value}
			-min_support {set min_support $value}
			-min_seq_size {set min_seq_size $value}
			-threads {set threads $value}
			default {
				lappend opts $key $value
			}
		}
	}
	set bamfile [file_absolute $bamfile]
	set refseq [refseq $refseq]
	if {$resultfile eq ""} {
		set root sniffles-[file_rootname $bamfile]
		set resultfile [file dir $bamfile]/sv-$root.tsv.zst
	} else {
		set root [file_rootname $resultfile]
	}
	if {$sample eq ""} {set sample $root}
	set destdir [file dir $resultfile]
	set resultanalysisinfo [analysisinfo_file $resultfile]
	set vcffile [file root [gzroot $resultfile]].vcf
	# resultfiles
	set resultlist [list $resultfile $resultanalysisinfo $vcffile]
	if {$resultfiles} {
		return $resultlist
	}
	lappend skips -skip $resultlist
	# logfile
	job_logfile $destdir/sv_sniffles_[file tail $resultfile] $destdir $cmdline \
		{*}[versions sniffles gnusort8 zst os]
	# start
	## Produce sniffles SNP calls
	set keeppwd [pwd]
	cd $destdir
	set bamfileindex $bamfile.[indexext $bamfile]
	job sv_sniffles_$root.vcf {*}$skips -mem 1G -cores $threads \
	-deps {
		$bamfile $refseq $bamfileindex
	} -targets {
		$vcffile $vcffile.analysisinfo
	} -vars {
		bamfile sniffles opts refseq threads root sample min_support min_seq_size region threads
	} -code {
		analysisinfo_write $dep $target sample $sample varcaller sniffles varcaller_version [version sniffles] varcaller_cg_version [version genomecomb]
		if {$region ne ""} {
			set usebam [scratchfile].bam
			exec samtools view -h -b -1 --threads $threads $bamfile {*}[samregions $region $refseq] > $usebam
			exec samtools index $usebam
		} else {
			set usebam $bamfile
		}
		if {[catch {
			exec sniffles {*}$opts --threads $threads --skip_parameter_estimation \
				--genotype --cluster \
				--min_support $min_support --min_seq_size $min_seq_size \
				-m $usebam -v $target.temp 2>@ stderr >@ stdout
			file rename -force -- $target.temp $target
		} msg]} {
			# sniffles sometinmes (allways?) crashes on empty or small bam
			# only give error on larger bam, otherwise write empty result
			set temp [exec samtools view $usebam | head -100 | wc -l]
			if {$temp >= 100} {
				error $msg
			}
			file_write $target {}
		}
	}
	# 
	job sv_sniffles_vcf2tsv-$root {*}$skips -deps {
		$vcffile
	} -targets {
		$resultfile $resultanalysisinfo
	} -vars {
		sample split
	} -code {
		analysisinfo_write $dep $target
		if {[file size $dep] == 0} {
			file_write $target ""
		} else {
			cg vcf2tsv -locerror correct -split $split -removefields {name filter AN AC AF AA ExcessHet InbreedingCoeff MLEAC MLEAF NDA RPA RU STR} $dep $target.temp[gzext $target]
			file rename -force -- $target.temp[gzext $target] $target
		}
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
