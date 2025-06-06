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

# natural less then
proc loc_lt {chr1 chr2} {
	if {[loc_compare $chr1 $chr2] < 0} {return 1} else {return 0}
}

# natural less then or equal
proc loc_lte {chr1 chr2} {
	if {[loc_compare $chr1 $chr2] <= 0} {return 1} else {return 0}
}

# natural greater then
proc loc_gt {chr1 chr2} {
	if {[loc_compare $chr1 $chr2] > 0} {return 1} else {return 0}
}

# natural greater then or equal
proc loc_gte {chr1 chr2} {
	if {[loc_compare $chr1 $chr2] >= 0} {return 1} else {return 0}
}

proc distr_chr_compare {chr1 start limitchr limitstart limitend} {
	set comp [loc_compare $chr1 $limitchr]
	if {
		$comp == 0 
		|| ([string index $limitchr end] != "_") && [string match [chr_clip $limitchr]* [chr_clip $chr1]]
	} {
		if {$limitstart eq ""} {return 0}
		if {$start < $limitstart} {return -1}
		if {$start < $limitend} {return 0}
		return 1
	}
	return $comp
}

proc lloc_compare {loc1 loc2} {
	loc_compare [join $loc1 " "] [join $loc2 " "]
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
	if {[llength $args]  == 4} {
		foreach {a1 a2 ref alt} $args break
		set seq1 v
	} elseif {[llength $args]  == 5} {
		foreach {seq1 a1 a2 ref alt} $args break
	} else {error "wrong # of args for zyg, must be: ?sequenced? a1 a2 ref alt"}
	if {$seq1 eq "r"} {
		return r
	} elseif {$seq1 eq "u"} {
		return u
	} elseif {$alt eq "?"} {
		if {$a1 == "?" || $a2 == "?"} {return ?}
		if {$a1 ne $a2} {
			return t
		} elseif {$a1 eq $ref} {
			return r
		} else {
			return m
		}
	} elseif {$seq1 eq "v"} {
		if {$a1 == "?" || $a2 == "?"} {return v}
		if {$alt eq ""} {
			set alt {{}}
		} elseif {[llength $alt] <= 1} {
			set alt [split $alt ,]
		}
		if {$a1 in $alt} {
			if {$a1 eq $a2} {
				return m
			} elseif {$a2 eq $ref} {
				return t
			} else {
				return c
			}
		} elseif {$a2 in $alt} {
			if {$a1 eq $ref} {
				return t
			} else {
				return c
			}
		} elseif {$a1 eq $ref && $a2 eq $ref} {
			return r
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

