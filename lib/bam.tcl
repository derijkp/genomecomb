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
	set threads 1
	cg_options bam_index args {
		-skip {
			lappend skips -skip $value
		}
		-optional {
			set optional $value
		}
		-threads {
			set threads $value
		}
	} {bam} 1 1
	if {[file extension [gzroot $bam]] eq ".sam"} return
	set pre [lindex [split $bam -] 0]
	set root [file_rootname $bam]
	set bamindex $bam.[indexext $bam]
	job [job_relfile2name bamindex- [file tail $bam]] {*}$skips -optional $optional \
	-deps [list $bam] \
	-targets [list $bamindex] \
	-vars {
		threads
	} -code {
		putslog "making $target"
		set tempdir [dirtemp $dep]
		mklink $dep $tempdir/[file tail $dep]
		exec samtools index -@ $threads $tempdir/[file tail $dep] >@ stdout 2>@ stderr
		file rename -force $tempdir/[file tail $dep].[indexext $dep] $dep.[indexext $dep]
		file delete -force $tempdir
	}
}

proc sam_empty file {
	set f [open [list | samtools view --no-PG $file]]
	set read [gets $f line]
	catch {close $f}
	if {$read == -1} {
		return 1
	} else {
		return 0
	}
}

proc sam_filter {list} {
	array set a {
		PAIRED        1
		PROPER_PAIR   2
		UNMAP         4
		MUNMAP        8
		REVERSE      16
		MREVERSE     32
		READ1        64
		READ2       128
		SECONDARY   256
		QCFAIL      512
		DUP        1024
		SUPPL      2048
	}
	set filter 0
	set els {}
	foreach el [list_remove [split $list ",; "] {}] {
		lappend els $a($el)
	}
	format %.0f [lmath_sum $els]
}
