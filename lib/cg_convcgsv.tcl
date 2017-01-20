#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc convcgsv {srcfile dstfile} {
	catch {close $o} ; catch {gzclose $f}
	set f [gzopen $srcfile]
	set header [tsv_open $f comment]
	regexp {#FORMAT_VERSION[\t ]+([0-9.]+)} $comment temp version
	if {$version > 1.3} {
		set new 1
		set poss [list_cor $header {Id LeftStrand LeftChr LeftPosition LeftLength TransitionLength RightStrand RightChr RightPosition RightLength Distance score scoreType}]
	} else {
		set new 0
		set poss [list_cor $header {id	leftStrand	leftChr	leftPosition	leftLength transitionLength	rightStrand	rightChr	rightPosition	rightLength distance score scoreType}]
	}
	set nheader {chromosome begin end type id strand1 start1 end1 size strand2 chr2 start2 end2 score scoretype transitionLength}
	lappend nheader {*}[list_sub $header -exclude $poss]
	set tempfile [filetemp $dstfile]
	set o [open $tempfile w]
	puts $o [join $nheader \t]
	set last [expr {[llength $header]-1}]
	while {![eof $f]} {
		set line [split [gets $f] \t]
		set line [lrange $line 0 $last]
		if {![llength $line]} continue
		foreach {id strand1 chr1 end1 len1 transitionLength strand2 chr2 start2 len2 size score scoretype} [list_sub $line $poss] break
		if {$id eq "W-2"} continue
		set start1 [expr {$end1-$len1}]
		set end2 [expr {$start2+$len2}]
		if {![isint $size]} {
			set size -1
		}
		if {$chr1 ne $chr2} {
			set type trans
		} elseif {$strand1 ne $strand2} {
			set type inv
			set size [expr {max($size,0)}]
		} else {
			set type del
		}
		if {$type eq "trans"} {
			set s [expr {350-($end1-$start1)}]
			if {$s < 10} {set s 10}
			set end [expr {$end1+$s}]
		} elseif {$type eq "inv"} {
			set end [expr {$end1+$size}]
		} else {
			set end $start2
		}
		if {$end1 < $end} {
			set begin $end1
		} else {
			set begin $end
			set end $end1
		}
		set result [list $chr1 $begin $end $type $id $strand1 $start1 $end1 $size $strand2 $chr2 $start2 $end2 $score $scoretype $transitionLength]
		set temp [list_sub $line -exclude $poss]
		if {[llength $temp]} {
			lappend result {*}$temp
		}
		puts $o [join $result \t]
	}
	close $o
	gzclose $f
	cg select -s - $tempfile $dstfile
	file delete $tempfile
}

proc cg_convcgsv {args} {
	if {[llength $args] != 2} {
		errorformat convcgsv
	}
	foreach {srcfile dstfile} $args break
	set tempfile [filetemp $dstfile]
	convcgsv $srcfile $tempfile
	file rename -force $tempfile $dstfile
}
