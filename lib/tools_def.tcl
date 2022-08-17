proc adapterfile {{adapterfile {}}} {
	if {$adapterfile eq ""} {
		return $::externdir/adaptors.fa
	} elseif {![file exists $adapterfile]} {
		error "adapterfile $adapterfile does not exists"
	} else {
		return [file_absolute $adapterfile]
	}
}

proc ampliconsfile {sampledir {ref {}}} {
	if {$ref ne {}} {
		append ref _
	}
	lindex [jobglob $sampledir/reg_${ref}amplicons*.tsv $sampledir/reg_amplicons*.tsv $sampledir/reg_*_amplicons*.tsv] 0
}

proc targetfile {sampledir {ref {}}} {
	if {$ref ne {}} {
		append ref _
	}
	lindex [jobglob $sampledir/reg_${ref}targets*.tsv $sampledir/reg_targets*.tsv $sampledir/reg_*_targets*.tsv] 0
}

proc targetfile_job {sampledir {dbdir {}}} {
	upvar job_logdir job_logdir
	set dbdir [dbdir $dbdir]
	set ref [dbdir_ref $dbdir]
	set targetfile [targetfile $sampledir]
	if {[jobfileexists $targetfile]} {
		if {[job_file_or_link_exists $targetfile]} {
			set link [gzlink $targetfile]
			if {[file exists $link]} {
				# correct if linking to uncompressed where actual regfile is compressed
				file delete $targetfile
				gzmklink $link $targetfile
				return [gzfile $targetfile]
			}
		}
		return $targetfile
	}
	# get targetfile (if possible) based on info_capture.tsv file (old system)
	set capturefile $sampledir/info_capture.txt
	if {![file exists $capturefile]} {
		if {[file tail [file dir $sampledir]] eq "samples"} {
			set take2 [gzfile [file dir [file dir $sampledir]]/reg_${ref}_targets.tsv]
			set capturefile [gzfile [file dir [file dir $sampledir]]/info_capture.tsv]
		} else {
			set take2 [gzfile [file dir $sampledir]/reg_${ref}_targets.tsv]
			set capturefile [gzfile [file dir $sampledir]/info_capture.tsv]
		}
		if {[jobfileexists $take2]} {
			set targetfile $sampledir/reg_${ref}_targets.tsv[gzext $take2]
			set take2 [find_link $take2 1]
			mklink_asjob $take2 $targetfile
			return $targetfile
		}
	}
	if {[file exists $capturefile]} {
		set targetfile $sampledir/reg_${ref}_targets.tsv
		set capture [string trim [file_read $capturefile]]
		set oritargetfile [gzfile $dbdir/extra/reg_${ref}_exome_$capture.tsv]
		if {![file exists $oritargetfile]} {
			array set transa {seqcapv3 SeqCap_EZ_v3 sure4 SureSelectV4 sure5 SureSelectV5 sure5utr SureSelectV5UTR}
			set capture [get transa($capture) $capture]
			set oritargetfile [gzfile $dbdir/extra/reg_${ref}_exome_$capture.tsv]
		}
		if {[file exists $oritargetfile]} {
			mklink_asjob $oritargetfile $targetfile[gzext $oritargetfile] 1
			return $targetfile[gzext $oritargetfile]
		}
	}
	set ampliconsfile [ampliconsfile $sampledir]
	if {[jobfileexists $ampliconsfile]} {
		job reports_amplicons2targetfile -deps {$ampliconsfile} -targets {$targetfile.zst} -vars {sample dbdir ref} -code {
			cg regcollapse $dep | cg zst > $target
		}
		return $targetfile
	}
	jobglob $sampledir/reg_${ref}_targets.tsv
}
