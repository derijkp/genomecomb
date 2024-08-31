proc bam_index {args} {
	set threads 1
	cg_options bam_index args {
		-threads {
			set threads $value
		}
	} {bam} 1 ...
	set bams [list $bam {*}$args]
	foreach dep $bams {
		set ext [file extension $dep]
		if {$ext in ".bam .cram" || ($ext eq ".gz" && [file extension [file root $dep]] eq ".sam")} {
			set ext [indexext [gzroot $dep]]
			set target $dep.$ext
			if {[file exists $target]} continue
			putslog "making $target"
			set tempdir [dirtemp $dep]
			mklink $dep $tempdir/[file tail $dep]
			exec samtools index -@ $threads $tempdir/[file tail $dep] >@ stdout 2>@ stderr
			exec touch -r $dep $tempdir/[file tail $dep].$ext
			file rename -force $tempdir/[file tail $dep].$ext $dep.$ext
			file delete -force $tempdir
		}
	}
}

proc bam_index_job {args} {
	upvar job_logdir job_logdir
	set skips {}
	set optional 1
	set threads 1
	cg_options bam_index args {
		-threads {
			set threads $value
		}
		-skip {
			lappend skips -skip $value
		}
		-optional {
			set optional $value
		}
	} {bam} 1 ...
	set bams [list $bam {*}$args]
	foreach bam $bams {
		set ext [file extension $bam]
		if {$ext in ".bam .cram" || ($ext eq ".gz" && [file extension [file root $dep]] eq ".sam")} {
			# can only make index for bam, cram, or gz compressed sam
		} else continue
		set bamindex $bam.[indexext [gzroot $bam]]
		if {[file exists $bamindex] && [file size $bamindex] == 0} {
			file rename $bamindex $bamindex.old
		}
		job [job_relfile2name bamindex- [file tail $bam]] {*}$skips -optional $optional \
		-deps [list $bam] \
		-targets [list $bamindex] \
		-vars {
			threads
		} -code {
			bam_index -threads $threads $dep
		}
	}
}

proc cg_bam_index {args} {
	set args [job_init {*}$args]
	bam_index_job {*}$args
	job_wait
}
