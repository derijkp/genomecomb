proc version_npinv {} {
	set version ?
	set npinv [findjar npInv]
	catch {execjar npInv -h} help
	if {![regexp {Version:[ \t]*([^ ][^\n]*)} $help temp version]} {
		regexp {npInv(.*)} [file root [file tail $npinv]] temp version
	}
	return $version
}

proc sv_npinv_job {args} {
	upvar job_logdir job_logdir
	set cmdline "[list cd [pwd]] \; [list cg sv_npinv {*}$args]"
	set refseq {}
	set opts {}
	set min 50
	set max 10000000
	set split 1
	set threads 2
	set cleanup 0
	set resultfiles 0
	set rootname {}
	set skips {}
	set resultfile {}
	set region {}
	set sample {}
	cg_options sv_npinv args {
		-refseq {
			set refseq $value
		}
		-split {
			set split $value
		}
		-exome {
			set exome $value
		}
		-min - -npinv-min {
			set min $value
		}
		-max - -npinv-max {
			set max $value
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
		-preset {
			# not used
		}		
		-region {
			set region $value
		}
		-sample {
			set sample $value
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
	if {$resultfile eq ""} {
		set root npinv-[file_rootname $bamfile]
		set resultfile [file dir $bamfile]/sv-$root.tsv.zst
	} else {
		set root [file_rootname $resultfile]
	}
	if {$sample eq ""} {set sample $root}
	set resultanalysisinfo [analysisinfo_file $resultfile]
	set destdir [file dir $resultfile]
	set vcffile [file root [gzroot $resultfile]].vcf
	# resultfiles
	set resultlist [list $resultfile {} {} $vcffile $resultanalysisinfo]
	if {$resultfiles} {
		return $resultlist
	}
	lappend skips -skip $resultlist
	# logfile
	job_logfile $destdir/sv_npinv_[file tail $resultfile] $destdir $cmdline \
		{*}[versions npinv gnusort8 zst os]
	# start
	## Produce npinv sv calls
	set bamfileindex $bamfile.[indexext $bamfile]
	job sv_npinv-$root.vcf {*}$skips -mem 5G \
	-deps {
		$bamfile $bamfileindex
	} -targets {
		$vcffile $vcffile.analysisinfo
	} -vars {
		bamfile vcffile opts root min max region sample refseq
	} -code {
		analysisinfo_write $dep $vcffile sample $sample varcaller npinv varcaller_version [version npinv] varcaller_cg_version [version genomecomb]
		if {$region ne ""} {
			set usebam [scratchfile].bam
			catch_exec samtools view -h -b -1 $bamfile {*}[samregions $region $refseq] > $usebam
			exec samtools index $usebam
		} else {
			set usebam $bamfile
		}
		execjar -javaversion 1.8 npInv {*}$opts \
			--input $usebam \
			--output $vcffile.temp \
			--min $min -max $max
		file rename -force -- $vcffile.temp $vcffile
	}
	# 
	job sv_npinv-vcf2tsv-$root {*}$skips -deps {
		$vcffile $vcffile.analysisinfo
	} -targets {
		$resultfile $resultanalysisinfo
	} -vars {
		sample split resultfile
	} -code {
		analysisinfo_write $dep $target
		cg vcf2tsv -split $split -removefields {
			name filter AN AC AF AA ExcessHet InbreedingCoeff MLEAC MLEAF NDA RPA RU STR
		} $dep $target.temp[gzext $target]
		file rename -force -- $target.temp[gzext $target] $target
	}
	# cleanup
	return $resultlist
}

proc cg_sv_npinv {args} {
	set args [job_init {*}$args]
	sv_npinv_job {*}$args
	job_wait
}

