#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

# comparepos replaced by loc_compare in extension
# < 0 if comp1 < comp2
# > 0 if comp1 > comp2
# chr_compare is currently also implemented by loc_compare in extension

proc chr_compare {chr1 chr2} {
	loc_compare $chr1 $chr2
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

proc chr_clip chr {
	regsub {^[Cc][Hh][Rr]-?} $chr {} chr
	if {$chr eq "MT"} {set chr M}
	return $chr
}

