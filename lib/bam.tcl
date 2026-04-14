proc bam_chrs {bamfile} {
	set bamheader [catch_exec samtools view --no-PG -H $bamfile]
	list_unmerge [regexp -all -inline {SN:([^\t]+)} $bamheader] 1 result
	return $result
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

proc usebam {bam {max cramv3} {index 1}} {
	set ext [file extension $bam]
	if {$ext eq ".bam"} {
		set tempbam [tempdir]/[file root [file tail $bam]].bam
		mklink $bam $tempbam
		mklink $bam.bai $tempbam.bai
	} elseif {$max eq "cramv3"} {
		set fh [open $bam rb]
		set hdr [read $fh 6]
		close $fh
		binary scan $hdr a4cc magic major minor
		if {$magic ne "CRAM"} {
			error "$bam is not a cram file"
		} elseif {$major == 3 && $minor == 1} {
			set tempbam [tempdir]/[file root [file tail $bam]].bam
			exec samtools view -h -b $bam > $tempbam
			exec samtools index $tempbam
		} else {
			set tempbam [tempdir]/[file root [file tail $bam]]$ext
			mklink $bam $tempbam
			mklink $bam.crai $tempbam.crai
		}
	} elseif {$max eq "bam"} {
		set tempbam [tempdir]/[file root [file tail $bam]].bam
		exec samtools view -h -b $bam > $tempbam
		exec samtools index $tempbam
	}
	return $tempbam
}
