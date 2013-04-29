#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc annot_region_init {regfile} {
	global annot
	catch {annot_region_close $regfile}
	set fr [gzopen [gzfile $regfile]]
	set header [tsv_open $fr]
	set poss [tsv_basicfields $header 3]
	if {[lsearch $poss -1] != -1} {
		set poss [list_cor $header {chr patchstart pos}]
	}
	if {[lsearch $poss -1] != -1} {
		error "Wrong header"
	}
	set line [split [gets $fr] \t]
	foreach {chr2 start2 end2} [list_sub $line $poss] break
	set annot(reg,$regfile) [list $fr $chr2 $start2 $end2 $poss]
}

proc annot_region_get {regfile chr1 start1 end1} {
	global annot
	foreach {fr chr2 start2 end2 poss2} $annot(reg,$regfile) break
	set cur 0
	while 1 {
		# putsvars chr1 chr2 start1 end1 start2 end2
		set chrcompar [chr_compare $chr2 $chr1]
		if {$chrcompar > 0} break
		if {$chrcompar == 0} {
			if {$start2 > $end1} break
			if {$start2 == $end1} {
				if {($start2 == $end2) || ($start1 == $end1)} {set cur 1}
				break
			}
			if {$start2 > $start1} {
				set cur 1
				break
			} elseif {$end2 >= $start1} {
				set cur 1
				break
			}
		}
		set line2 [getnotempty $fr]
		if {![llength $line2]} break
		foreach {chr2 start2 end2} [list_sub $line2 $poss2] break
	}
	set annot(reg,$regfile) [list $fr $chr2 $start2 $end2 $poss2]
	return $cur
}

proc annot_region_in {regfile chr1 start1 end1} {
	global annot
	foreach {fr chr2 start2 end2 poss2} $annot(reg,$regfile) break
	set cur 0
	while 1 {
		# putsvars chr1 chr2 start1 end1 start2 end2
		set chrcompar [chr_compare $chr2 $chr1]
		if {$chrcompar > 0} break
		if {$chrcompar == 0} {
			if {$start2 > $end1} break
			if {$start2 == $end1} {
				if {($start2 == $end2) || ($start1 == $end1)} {set cur 1}
				break
			}
			if {($start1 >= $start2) && ($end1 <= $end2)} {
				set cur 1
				break
			}
			if {($start2 > $start1) || ($end2 >= $start1)} {
				break
			}
		}
		set line2 [getnotempty $fr]
		if {![llength $line2]} break
		foreach {chr2 start2 end2} [list_sub $line2 $poss2] break
	}
	set annot(reg,$regfile) [list $fr $chr2 $start2 $end2 $poss2]
	return $cur
}

proc annot_region_close {regfile} {
	global annot
	foreach {fr} $annot(reg,$regfile) break
	catch {close $fr}
	unset annot(reg,$regfile)
}

proc annot_coverage_init {dir sample} {
	global annot
	catch {annot_coverage_close $dir}
	set annot(cov,$dir,$sample) [list -1 {} 0 {} {} {}]
}

proc annot_coverage_get {dir sample chr begin {force 0}} {
	global annot
	set present 1
	foreach {curchr chrfile present poss type obj} [get annot(cov,$dir,$sample) {{} {} 0}] break
	if {$chr ne $curchr || $force} {
		if {$present} {
			switch $type {
				{} {
					tsv_index_close $chrfile offset
				}
				bcol {
					foreach {refbcol covbcol} $obj break
					bcol_close $refbcol ; bcol_close $covbcol
				}
			}
		}
		set nchr [chr_clip $chr]
		set chrfile [multicompar_reannot_find $dir $sample coverage/coverage-*-$nchr.bcol coverage/coverage-*-chr$nchr.bcol coverage/coverage-$nchr-*.bcol coverage/coverage-chr$nchr-*.bcol]
		if {[file exists $chrfile]} {
			set reffile [multicompar_reannot_find $dir $sample coverage/refScore-*-$nchr.bcol coverage/refScore-*-chr$nchr.bcol coverage/refScore-$nchr-*.bcol coverage/refScore-chr$nchr-*.bcol]
			set type bcol
		} else {
			set chrfile [multicompar_reannot_find $dir $sample coverage/coverageRefScore-$nchr-*.tsv coverage/coverage-$nchr-*.tsv coverage/coverageRefScore-chr$nchr-*.tsv coverage/coverage-chr$nchr-*.tsv]
			set type tsv
		}
		if {[file exists $chrfile]} {
			set present 1
		} elseif {[inlist {Y chrY} $chr]} {
			set present 0
		} else {
			puts stderr "coverage(RefScore) file not found ($dir/coverage/coverage(RefScore)-*-$nchr.(tsv|bcol))"
			set present 0
		}
		if {!$present} {
			set annot(cov,$dir,$sample) [list $chr {} $present $poss $type {}]
			return {u u}
		}
		if {$type eq "bcol"} {
			if {[file exists $reffile]} {
				set refbcol [bcol_open $reffile]
			} else {
				set refbcol {}
			}
			set covbcol [bcol_open $chrfile]
			set obj [list $refbcol $covbcol]
			set annot(cov,$dir,$sample) [list $chr $chrfile $present $poss $type $obj]
		} else {
			tsv_index_open $chrfile offset 1
			set header [tsv_index_header $chrfile]
			set rpos [lsearch $header refScore]
			set cpos [lsearch $header uniqueSequenceCoverage]
			if {$cpos == -1} {set cpos [lsearch $header coverage]}
			set poss [list $rpos $cpos]
			set annot(cov,$dir,$sample) [list $chr $chrfile $present $poss $type {}]
		}
	}
	if {!$present} {return {u u}}
	if {$type eq "bcol"} {
		foreach {refbcol covbcol} $obj break
		if {$refbcol eq ""} {
			set refscore ?
		} else {
			set refscore [bcol_get $refbcol $begin $begin]
		}
		set coverage [bcol_get $covbcol $begin $begin]
		return [list $refscore $coverage]
	} else {
		ifcatch {tsv_index_get $chrfile offset $begin} cline -regexp {
			"not found in" {set cline {u u u}}
		}
		return [list_sub $cline $poss]
	}
}

proc annot_coverage_close {dir sample} {
	global annot
	foreach {curchr chrfile present poss type obj} [get annot(cov,$dir,$sample) {{} {} 0 {} {} {}}] break
	switch $type {
		tsv {
			tsv_index_close $chrfile offset
		}
		bcol {
			foreach {refbcol covbcol} $obj break
			if {$refbcol ne ""} {bcol_close $refbcol}
			bcol_close $covbcol
		}
	}
	unset annot(cov,$dir,$sample)
}

proc annot_varall_init {varallfile sample header} {
	global annot
	set af [gzopen $varallfile]
	set aheader [tsv_open $af]
	set loc2poss [tsv_basicfields $aheader 6 0]
	set poss1 {}
	set poss2 {}
	set i 0
	foreach field $aheader {
		set pos [lsearch $header $field-$sample]
		if {[inlist $loc2poss $pos]} continue
		if {$pos != -1} {
			lappend poss1 $pos
			lappend poss2 $i
		}
		incr i
	}
	set loc2poss [lrange $loc2poss 0 2]
	set line2 [split [gets $af] \t]
	set loc2 [list_sub $line2 $loc2poss]
	set annot(varallinfo,$varallfile,$sample) [dict create af $af poss1 $poss1 poss2 $poss2 line2 $line2 loc2 $loc2 loc2poss $loc2poss]
}

proc annot_varall_annot {varallfile sample loc force lineVar} {
	global annot
	upvar $lineVar line
	set info $annot(varallinfo,$varallfile,$sample)
	set af [dict get $info af]
	set line2 [dict get $info line2]
	set loc2 [dict get $info loc2]
	set loc2poss [dict get $info loc2poss]
	while {![eof $af]} {
		set d [loc_compare $loc $loc2]
		if {$d == 0} {
			set poss1 [dict get $info poss1]
			set cur [list_sub $line $poss1]	
			set newvalues [list_sub $line2 [dict get $info poss2]]
			foreach c $cur pos $poss1 newvalue $newvalues {
				if {(!$force) && $c ne "?"} continue
				lset line $pos $newvalue
			}
			break
		} elseif {$d < 0} {
			break
		}
		set line2 [split [gets $af] \t]
		set loc2 [list_sub $line2 $loc2poss]
	}
	dict set annot(varallinfo,$varallfile,$sample) line2 $line2
	dict set annot(varallinfo,$varallfile,$sample) loc2 $loc2
}

proc annot_varall_close {varallfile sample} {
	global annot
	set info $annot(varallinfo,$varallfile,$sample)
	close [dict get $info af]
	unset annot(varallinfo,$varallfile,$sample)
}
