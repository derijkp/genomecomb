proc bam_chrs {bamfile} {
	set bamheader [exec samtools view -H $bamfile]
	list_unmerge [regexp -all -inline {SN:([^\t]+)} $bamheader] 1 result
	return $result
}

proc bam_index_job {args} {
	upvar job_logdir job_logdir
	set skips {}
	cg_options bam_index args {
		-skip {
			lappend skips -skip $value
		}
	} {bam} 1 1
	if {[file extension $bam] eq ".sam"} return
	set pre [lindex [split $bam -] 0]
	set root [file_rootname $bam]
	set dir [file dir $bam]
	set ext [file extension $bam]
	if {$ext eq ".cram"} {
		job cramindex-$pre-$root -deps [list $bam] -targets [list $bam.crai] -code {
			exec samtools index $dep >@ stdout 2>@ stderr
			puts "making $target"
		}
	} else {
		job bamindex-$pre-$root -deps [list $bam] -targets [list $bam.bai] -code {
			exec samtools index $dep >@ stdout 2>@ stderr
			puts "making $target"
		}
	}
}
