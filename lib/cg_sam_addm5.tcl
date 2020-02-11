proc sam_header_addm5 {header {refseq {}}} {
	if {![regexp {M5:} $header]} {
		dbdir $refseq
		if {![file exists [get ::env(REF_PATH) /complgen/refseq]/mapping.tsv]} {
			error "Could not find reference md5 mapping file (for cram): specify dbdir, reference or use REF_PATH"
		}
		set temp [file_read $::env(REF_PATH)/mapping.tsv]
		foreach line [split $temp \n] {
			set line [split $line \t]
			if {[llength $line] < 2} continue
			set a([lindex $line 1]) [lindex $line 2]
		}
		set newheader {}
		foreach line [split $header \n] {
			if {[regexp ^@SQ $line]} {
				if {![regexp {SN:([^ \t]+)} $line temp name]} {
					error "@SQ field without SN: for line $line"
				}
				set name [chr_clip $name]
				if {![info exists a($name)]} {
					set name chr$name
					if {![info exists a($name)]} {
						error "no md5 mapping found for sequence $name in $::env(REF_PATH)/mapping.tsv"
					}
				}
				append line "\tM5:$a($name)"
				lappend newheader $line
			} else {
				lappend newheader $line
			}
		}
		set header [join $newheader \n]
	}
	return $header
}

proc cg_sam_addm5 {{refseq {}}} {
	set header {}
	while {[gets stdin line] != -1} {
		if {[string index $line 0] ne "@"} break
		append header $line\n
	}
	puts [sam_header_addm5 $header $refseq]
	puts $line
	fcopy stdin stdout
}
