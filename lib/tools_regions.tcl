proc region2bins {start end {level 0}} {
	if {$level == 0} {
		set bins 0
		set stop 0
	} else {
		set bins {}
		array set levels {1 1 2 9 3 73 4 585 5 4681}
		set stop $levels($level)
	}
	set pre 4681
	set start [expr {$start >> 14}]
	set end [expr {$end >> 14}]
	while {$pre > $stop} {
		for {set i $start} {$i <= $end} {incr i} {
			lappend bins [expr {$pre+$i}]
		}
		set pre [expr {$pre >> 3}]
		set start [expr {$start >> 3}]
		set end [expr {$end >> 3}]
	}
	return [lsort -integer $bins]
}

proc overlap {start1 end1 start2 end2} {
	if {$start2 >= $end1} {return [expr {$end1-$start2}]}
	if {$end2 < $start1} {return [expr {$end2-$start1}]}
	if {$start2 > $start1} {set start1 $start2}
	if {$end2 < $end1} {set end1 $end2}
	expr {$end1-$start1}
}

proc reg_compare {loc1 loc2} {
	if {![llength $loc1]} {return 1}
	if {![llength $loc2]} {return -1}
	foreach {chr1 start1 end1} $loc1 break
	foreach {chr2 start2 end2} $loc2 break
	set chrcomp [loc_compare $chr1 $chr2]
	if {$chrcomp != 0} {return $chrcomp}
	if {$start1 == $start2} {return 0}
	if {$start2 >= $end1} {return [expr {$end1 - $start2 -1}]}
	if {$end2 <= $start1} {return [expr {$start1 - $end2 + 1}]}
	return 0
}

proc tempbed {regionfile {reffile {}}} {
	if {[file extension [gzroot $regionfile]] eq ".bed"} {
		set bedfile $regionfile
	} elseif {$reffile eq ""} {
		set bedfile [tempfile].bed
		tsv2bed $regionfile $bedfile
	} else {
		set bedfile [tempfile].bed
		gatkworkaround_tsv2bed $regionfile $reffile $bedfile
	}
	return $bedfile
}

proc reg_isempty {file} {
	if {$file eq ""} {return 0}
	if {[file extension $file] eq ".bed"} {
		if {[file size $file] == 0} {return 1} else {return 0}
	} else {
		set f [gzopen $file]
		set header [tsv_open $f]
		set isempty [eof $f]
		gzclose $f
		return $isempty
	}
}
