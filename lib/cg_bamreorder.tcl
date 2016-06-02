proc cg_bamreorder {src dest ref} {
	if {![file exists $ref.fai]} {
		exec samtools faidx $ref
	}
	set index [csv_parse [file_read $ref.fai] \t]
	set refchrs [list_subindex $index 0]
	set reflens [list_subindex $index 3]
	# check bam file to see if we need to rename some contigs
	set bamheader [exec samtools view -H $src]
	set bamchrs {}
	set bamlens {}
	foreach line [split $bamheader \n] {
		if {[regexp {SN:([^ ]+)\tLN:([0-9]+)} $line temp chr len]} {
			lappend bamchrs $chr
			lappend bamlens $len
		}
	}
	# putsvars bamchrs bamlens
	set torename {}
	foreach chr $bamchrs len $bamlens {
		set pos [lsearch $refchrs $chr]
		set newchr $chr
		if {$pos == -1} {
			if {[regexp ^chr $chr]} {
				set newchr [string range $chr 3 end]
				set pos [lsearch $refchrs $newchr]
			} else {
				set newchr chr$chr
				set pos [lsearch $refchrs $newchr]
			}
		}
		if {$pos == -1} {
			puts stderr "WARNING: contig $chr was not found in ref"
		} else {
			set reflen [lindex $reflens $pos]
			if {$len != $reflen} {
				puts stderr "WARNING: dropping contig $chr: has different size from $newchr in ref"
				set newchr _rem_$chr
			}
			if {$newchr ne $chr} {
				lappend torename SN:$chr SN:$newchr
			}
		}
	}

	if {[llength $torename]} {
		set bamheader [string_change $bamheader $torename]
		set tempfile [tempfile]
		file_write $tempfile $bamheader
		set tempsrc $dest.temp
		exec samtools reheader $tempfile $src > $tempsrc
		exec samtools index $tempsrc
	} else {
		set tempsrc $src
	}
	picard ReorderSam R=$ref I=$tempsrc O=$dest.temp2 ALLOW_INCOMPLETE_DICT_CONCORDANCE=true VALIDATION_STRINGENCY=LENIENT 2>@ stderr >@ stdout
	exec samtools index $dest.temp2
	file delete $dest.temp  $dest.temp.bai
	file rename $dest.temp2 $dest
	file rename $dest.temp2.bai $dest.bai
}
