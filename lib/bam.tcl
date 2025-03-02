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
