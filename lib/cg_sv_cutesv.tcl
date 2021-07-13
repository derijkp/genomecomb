proc version_cuteSV {} {
	set version ?
	lindex [exec cuteSV --version] end
}

proc sv_cuteSV_job {args} {
	upvar job_logdir job_logdir
	set cmdline "[list cd [pwd]] \; [list cg sv_cuteSV {*}$args]"
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
	set max_cluster_bias_DEL 100
	set diff_ratio_merging_DEL 0.3
	set preset {}
	set resultfile {}
	set region {}
	set sample {}
	cg_options sv_cuteSV args {
		-refseq {
			set refseq $value
		}
		-split {
			set split $value
		}
		-preset {
			set preset $value
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
		-cuteSVopts {
			lappend opts {*}$value
		}
	} {bamfile resultfile} 1 2
	foreach {key value} [specialopts -cuteSV] {
		switch $key {
			-regmincoverage {set regmincoverage $value}
			-min_support {set min_support $value}
			-min_seq_size {set min_seq_size $value}
			-max_cluster_bias_DEL {set max_cluster_bias_DEL $value}
			-diff_ratio_merging_DEL {set diff_ratio_merging_DEL $value}
			-threads {set threads $value}
			default {
				lappend opts $key $value
			}
		}
	}
	set bamfile [file_absolute $bamfile]
	set refseq [refseq $refseq]
	if {$resultfile eq ""} {
		if {$preset eq ""} {
			set root cuteSV-[file_rootname $bamfile]
		} else {
			set root cuteSV_${preset}-[file_rootname $bamfile]
		}
		set resultfile [file dir $bamfile]/sv-$root.tsv.zst
	} else {
		set root [file_rootname $resultfile]
	}
	if {$sample eq ""} {set sample $root}
	set destdir [file dir $resultfile]
	set resultanalysisinfo [analysisinfo_file $resultfile]
	set vcffile [file root [gzroot $resultfile]].vcf
	
	if {$preset eq "pacbioclr" || $preset eq "pacbio"} {
		lappend opts --max_cluster_bias_INS 100 --diff_ratio_merging_INS	0.3 \
			--max_cluster_bias_DEL 200 --diff_ratio_merging_DEL	0.5
	} elseif {$preset eq "pacbioccs"} {
		lappend opts --max_cluster_bias_INS 1000 --diff_ratio_merging_INS	0.9 \
			--max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL	0.5
	} elseif {$preset eq "ont"} {
		lappend opts --max_cluster_bias_INS 100 --diff_ratio_merging_INS	0.3 \
			--max_cluster_bias_DEL 100 --diff_ratio_merging_DEL	0.3
	} else {
		lappend opts --max_cluster_bias_DEL $max_cluster_bias_DEL --diff_ratio_merging_DEL	$max_cluster_bias_DEL
	}
	# resultfiles
	set resultlist [list $resultfile $resultanalysisinfo $vcffile]
	if {$resultfiles} {
		return $resultlist
	}
	# logfile
	job_logfile $destdir/sv_cuteSV_[file tail $resultfile] $destdir $cmdline \
		{*}[versions cuteSV gnusort8 zst os]
	# start
	## Produce cuteSV SNP calls
	set keeppwd [pwd]
	cd $destdir
	set bamfileindex $bamfile.[indexext $bamfile]
	job sv_cutesv_$root.vcf {*}$skips -mem 1G -cores $threads \
	-skip [list $resultfile $resultanalysisinfo] \
	-deps {
		$bamfile $refseq $bamfileindex
	} -targets {
		$vcffile $vcffile.analysisinfo
	} -vars {
		bamfile cuteSV opts refseq threads root sample min_support min_seq_size region
	} -code {
		set workdir [scratchdir]
		analysisinfo_write $dep $target sample $sample varcaller cuteSV varcaller_version [version cuteSV] varcaller_cg_version [version genomecomb]
		if {$region ne ""} {
			set usebam [scratchfile].bam
			exec samtools view -h -b -1 --threads $threads $bamfile {*}[samregions $region $refseq] > $usebam
			exec samtools index $usebam
		} else {
			set usebam $bamfile
		}
		exec cuteSV {*}$opts --threads $threads --genotype \
			--min_support $min_support \
			$usebam $refseq $target.temp $workdir 2>@ stderr >@ stdout
		file rename -force -- $target.temp $target
	}
	# 
	job sv_cutesv_vcf2tsv-$root {*}$skips -deps {
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

proc cg_sv_cuteSV {args} {
	set args [job_init {*}$args]
	sv_cuteSV_job {*}$args
	job_wait
}
