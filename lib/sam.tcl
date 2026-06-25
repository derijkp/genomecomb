#1	QNAME	Query (pair) NAME
#2	FLAG	bitwise FLAG
#3	RNAME	Reference sequence NAME
#4	POS	1-based leftmost POSition/coordinate of clipped sequence
#5	MAPQ	MAPping Quality (Phred-scaled)
#6	CIAGR	extended CIGAR string
#7	MRNM	Mate Reference sequence NaMe (\u2018=\u2019 if same as RNAME)
#8	MPOS	1-based Mate POSistion
#9	ISIZE	Inferred insert SIZE
#10	SEQ	query SEQuence on the same strand as the reference
#11	QUAL	query QUALity (ASCII-33 gives the Phred base quality)
#12	OPT	variable OPTional fields in the format TAG:VTYPE:VALUE

#0x0001	the read is paired in sequencing
#0x0002	the read is mapped in a proper pair
#0x0004	the query sequence itself is unmapped
#0x0008	the mate is unmapped
#0x0010	strand of the query (1 for reverse)
#0x0020	strand of the mate
#0x0040	the read is the first read in a pair
#0x0080	the read is the second read in a pair
#0x0100	the alignment is not primary
#0x0200	the read fails platform/vendor quality checks
#0x0400	the read is either a PCR or an optical duplicate

proc sam_paired flags {
	expr {($flags & 0x001) != 0}
}

proc sam_mapped_pair flags {
	expr {($flags & 0x002) != 0}
}

proc sam_unmapped flags {
	expr {($flags & 0x004) != 0}
}

proc sam_unmapped_mate flags {
	expr {($flags & 0x008) != 0}
}

proc sam_reverse flags {
	expr {($flags & 0x010) != 0}
}

proc sam_reverse_mate flags {
	expr {($flags &	0x020) != 0}
}

proc sam_pair_first flags {
	expr {($flags & 0x040) != 0}
}

proc sam_pair_second flags {
	expr {($flags & 0x080) != 0}
}

proc sam_notprimary flags {
	expr {($flags & 0x100) != 0}
}

proc sam_failed flags {
	expr {($flags & 0x200) != 0}
}

proc sam_duplicate flags {
	expr {($flags & 0x400) != 0}
}

if 0 {

proc sam_flag flag {
	foreach {mask var} {
		0x001	paired
		0x002	mapped_pair
		0x004	unmapped
		0x008	unmapped_mate
		0x010	reverse
		0x020	reverse_mate
		0x040	pair_first
		0x080	pair_second
		0x100	notprimary
		0x200	failed
		0x400	duplicate
	} {
		if {[expr {$flags & $mask}]} {set a($var) 1} else {set a($var) 0}
	}
}
foreach {QNAME FLAG RNAME POS MAPQ CIGAR MRNM MPOS ISIZE SEQ QUAL OPT} $line break
}

proc sam_readgroupdata_fix {readgroupdata} {
	set result {}
	foreach {key value} $readgroupdata {
		if {[string length $key] != 2} {
			set value $key=$value
			set key CO
		}
		lappend result $key $value
	}
	return $result
}

proc sam_readgroup {readgroupdata sample {format RG} args} {
	if {$readgroupdata eq "-"} {return ""}
	unset -nocomplain a
	array set a [list PL illumina ID $sample PU $sample SM $sample]
	foreach {key value} $args {
		set a($key) $value
	}
	if {$a(PL) eq "illumina" && ![info exists a(LB)]} {
		set a(LB) solexa-123
	}
	foreach {key value} $readgroupdata {
		if {[string length $key] != 2} {
			set value $key=$value
			set key CO
		}
		set a($key) $value
	}
	if {$format eq "RG"} {
		set rg @RG\\tID:$a(ID)
		unset a(ID)
		foreach {key value} [array get a] {
			if {[string length $key] != 2} {
				set value $key=$value
				set key CO
			}
			append rg "\\t$key:$value"
		}
	} else {
		set rg [list ID:$a(ID)]
		unset a(ID)
		foreach {key value} [array get a] {
			if {[string length $key] != 2} {
				set value $key=$value
				set key CO
			}
			lappend rg "$key:$value"
		}
		if {$format eq "list"} {
			return $rg
		}
	}
	return $rg
}

proc sam_header_get {samfile} {
	if {![file size $samfile]} {return ""}
	if {[gziscompressed $samfile]} {
		set header [exec cg zcat $samfile | samtools view --no-PG -H]
	} else {
		set header [catch_exec samtools view --no-PG -H $samfile]
	}
	return $header
}

proc sam_header_extract {samheader key} {
	set result {}
	foreach line [split $samheader \n] {
		if {[regexp $key\\t $line]} {
			lappend result $line
		}
	}
	return $result
}

proc sam_header_add {header args} {
	unset -nocomplain a
	set newkeys {}
	set keys {}
	foreach line [split $header \n] {
		set key [lindex [split $line \t] 0]
		lappend a($key) $line
		list_addnew keys $key
	}
	foreach line $args {
		set key [lindex [split $line \t] 0]
		if {![info exists a($key)]} {list_addnew newkeys $key}
		list_addnew a($key) $line
	}
	if {[llength $newkeys]} {
		if {[lindex $keys end] eq @PG} {
			set keys [linsert $keys end-1 {*}$newkeys]
		} else {
			lappend keys {*}$newkeys
		}
	}
	set newheader {}
	foreach key $keys {
		foreach line $a($key) {
			lappend newheader $line
		}
	}
	return [join $newheader \n]
}

