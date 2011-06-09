#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

# < 0 if comp1 < comp2
# > 0 if comp1 > comp2
proc comparepos {comp1 comp2} {
	if {![isint [lindex $comp1 1]]} {return 1}
	if {![isint [lindex $comp2 1]]} {return -1}
	foreach {chr1 pos1 end1 type1} $comp1 break
	foreach {chr2 pos2 end2 type2} $comp2 break
	set chr1 [chr2num $chr1]
	set chr2 [chr2num $chr2]
	if {$chr1 ne $chr2} {
		return [expr {$chr1-$chr2}]
	} elseif {$pos1 != $pos2} {
		return [expr {$pos1-$pos2}]
	} elseif {$end1 != $end2} {
		return [expr {$end1-$end2}]
	} elseif {$type1 ne $type2} {
		if {$type1 < $type2} {return -1} else {return 1}
	} else {
		return 0
	}
}

proc compare_annot_getline {f} {
	set cur1 {}
	while {![eof $f]} {
		set cur1 [split [gets $f] \t]
		if {[llength $cur1]} break
	}
	return $cur1
}

proc lset_always {varName pos value} {
	upvar $varName var
	if {[catch {lset var $pos $value}]} {
		set var [list_concat $var [list_fill [expr {$pos - [llength $var]}] {}]]
		lappend var $value
	}
}

array set chrtrans {M 95 MT 95 X 96 Y 97}
proc chr2num {chr} {
	if {[isint $chr]} {return $chr}
	regsub ^chr $chr {} chr
	set nchr [get ::chrtrans($chr) $chr]
	return $nchr
}

proc chr2clean {chr} {
	if {[isint $chr]} {return $chr}
	regsub ^chr $chr {} chr
	if {$chr eq "MT"} {set chr M}
	return $chr
}

