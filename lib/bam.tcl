proc bam_chrs {bamfile} {
	set bamheader [exec samtools view --no-PG -H $bamfile]
	list_unmerge [regexp -all -inline {SN:([^\t]+)} $bamheader] 1 result
	return $result
}

proc bam_index {args} {
	foreach file $args {
		set ext [file extension $file]
		if {$ext in ".bam .cram" || ($ext eq ".gz" && [file extension [file root $file]] eq ".sam")} {
			exec samtools index $file >@ stdout 2>@ stderr
		}
	}
}

proc bam_index_job {args} {
	upvar job_logdir job_logdir
	set skips {}
	set optional 1
	cg_options bam_index args {
		-skip {
			lappend skips -skip $value
		}
		-optional {
			set optional $value
		}
	} {bam} 1 1
	if {[file extension [gzroot $bam]] eq ".sam"} return
	set pre [lindex [split $bam -] 0]
	set root [file_rootname $bam]
	set bamindex $bam.[indexext $bam]
	job bamindex-[file_part $pre-$root end] {*}$skips -optional $optional -deps [list $bam] -targets [list $bamindex] -code {
		putslog "making $target"
		exec samtools index $dep >@ stdout 2>@ stderr
	}
}

proc sam_empty file {
	if {[catch {exec samtools view --no-PG $file | head -1} test]} {
		return 0
	} else {
		return 1
	}
}
