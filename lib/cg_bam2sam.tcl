proc cg_bam2sam {args} {
	upvar job_logdir job_logdir
	set sortchr 1
	set samfile -
	set dbdir {}
	cg_options bam2sam args {
		-sortchr {
			set sortchr $value
		}
		-dbdir {
			set dbdir $value
		}
	} {bamfile samfile} 1 2
	set bamfile [file_absolute $bamfile]
	set ext [file extension $bamfile]
	if {$ext eq ".cram"} {
		refcram $dbdir
	}
	if {!$sortchr} {
		if {$samfile eq "-"} {
			exec samtools view -h $bamfile >@ stdout
		} else {
			exec samtools view -h $bamfile {*}[compresspipe $samfile] > $samfile
		}
	} else {
		set o [wgzopen $samfile]
		set header [split [exec samtools view -H $bamfile] \n]
		set pos 0
		foreach line $header {
			if {[regexp ^@SQ $line]} break
			puts $o $line
			incr pos
		}
		set list {}
		while 1 {
			lappend list $line
			incr pos
			set line [lindex $header $pos]
			if {![regexp ^@SQ $line]} break
		}
		set list [ssort -n $list]
		puts $o [join $list \n]
		set len [llength $header]
		while {$pos < $len} {
			incr pos
			puts $o $line
			set line [lindex $header $pos]
			incr pos
		}
		bam_index_job $bamfile
		foreach line $list {
			flush $o
			if {![regexp {\tSN:([^\t]+)\t} $line temp chr]} {
				error "format error in SQ line: $line"
			}
			exec samtools view $bamfile $chr >@ $o
		}
		gzclose $o
	}
}
