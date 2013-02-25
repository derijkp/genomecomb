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

proc compare {args} {
	if {[llength $args]  != 6} {error "compare only with 2 ids"}
	foreach {seq1 a11 a21 seq2 a12 a22} $args break
	if {$seq1 eq "r"} {
		if {$seq2 eq "r"} {
			return r
		} elseif {$seq2 eq "v"} {
			return df
		}
	} elseif {$seq1 eq "u" || $seq2 eq "u"} {
		return un
	} elseif {$seq1 eq "v"} {
		if {$seq2 eq "r"} {
			return df
		} elseif {$seq2 eq "v"} {
			if {($a11 eq $a12 && $a21 eq $a22) || ($a11 eq $a22 && $a21 eq $a12)} {
				return sm
			} else {
				return mm
			}
		}
	}
	return ?
}

if 0 {
compare u A A r A A
compare r A A r A A
compare v A C r A A
compare r A A v A C
compare v A C v A C
compare v A C v C A
compare v A C v A A
}

proc zyg {args} {
	if {[llength $args]  != 5} {error "wrong # of args for zyg: only one allowed"}
	foreach {seq1 a1 a2 ref alt} $args break
	if {$seq1 eq "r"} {
		return r
	} elseif {$seq1 eq "u"} {
		return u
	} elseif {$seq1 eq "v"} {
		if {$a1 eq $alt && $a2 eq $alt} {
			return m
		} elseif {($a1 eq $ref && $a2 eq $alt) || ($a1 eq $alt && $a2 eq $ref)} {
			return t
		} elseif {$a1 eq $alt || $a2 eq $alt} {
			return c
		} else {
			return o
		}
	}
	return ?
}

if 0 {
# t
zyg v A C A C 
zyg v C A A C 
# m
zyg v C C A C
# c
zyg v C G A C
zyg v G C A C
# o
zyg v G A A C
}

