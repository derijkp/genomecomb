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
	set poss [tsv_basicfields $header 3 $regfile]
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
			putslog "coverage(RefScore) file not found ($dir/coverage/coverage(RefScore)-*-$nchr.(tsv|bcol))"
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

proc annot_varall_getline {infoVar} {
	upvar $infoVar info
	set af [dict get $info af]
	set dchr [dict get $info dchr]
	set dend [dict get $info dend]
	set loc2poss [dict get $info loc2poss]
	set type2pos [dict get $info type2pos]
	if {[dict exists $info cache]} {
		set line2 [dict get $info cache]
	} else {
		set line2 [split [gets $af] \t]
	}
	set loc2 [list_sub $line2 $loc2poss]
	foreach {chr2 begin2} $loc2 break
	foreach {pchr pbegin} $loc2 break
	set result {}
	while 1 {
		set type2 [lindex $line2 $type2pos]
		if {$type2 in "del sub"} {
			foreach {dchr dbegin dend} $loc2 break
			dict set info dchr $dchr
			dict set info dend $dend
		} elseif {$type2 eq "snp"} {
			set result $line2
		}
		set line2 [split [gets $af] \t]
		if {[eof $af]} {
			if {![llength $line2]} break
		} elseif {![llength $line2]} continue
		set loc2 [list_sub $line2 $loc2poss]
		foreach {chr2 begin2} $loc2 break
		if {[loc_compare $chr2 $pchr] != 0 || $begin2 != $pbegin} {
			if {[llength $result]} break
		}
	}
	dict set info cache $line2
	return $result
}

proc annot_varall_init {varallfile sample header} {
	global annot

	catch {close $af} ; set af [gzopen $varallfile]
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
	set type2pos [lindex $loc2poss 3]
	set loc2poss [lrange $loc2poss 0 2]
	set info [dict create af $af poss1 $poss1 poss2 $poss2 loc2poss $loc2poss type2pos $type2pos dend -1 dchr ""]
	set line2 [annot_varall_getline info]
	dict set info line2 $line2
	set annot(varallinfo,$varallfile,$sample) $info

#	annot_varall_getline info
	return $info
}

proc annot_varall_annot {sampleaVar sample loc force lineVar} {
	global annot
	upvar $lineVar line
	upvar $sampleaVar samplea
	foreach {chr begin end type} $loc break
	set varallfile $samplea(varall,$sample)
	set info $annot(varallinfo,$varallfile,$sample)
	if {$type eq "snp"} {
		set loc [lrange $loc 0 2]
	} else {
		set loc [list $chr $begin [expr {$begin+1}]]
	}
	set dend [dict get $info dend]
	if {$dend != -1} {
		set dchr [dict get $info dchr]
		if {[loc_compare $dchr $chr] != 0} {
			set dend -1
			dict set info dend -1
		}
	}
	set af [dict get $info af]
	set line2 [dict get $info line2]
	set loc2poss [dict get $info loc2poss]
	set seq [lindex $line $samplea(seq,$sample)]
	set zyg [lindex $line $samplea(zyg,$sample)]
	while 1 {
		if {[llength $line2]} {
			set loc2 [list_sub $line2 $loc2poss]
			set d [lloc_compare $loc $loc2]
		} else {
			# no more data in varall
			set d -1
		}
		if {$d == 0} {
#puts ------
#putsvars loc sample line line2 info
			set a1 [lindex $line $samplea(a1,$sample)]
			set a2 [lindex $line $samplea(a2,$sample)]
			set poss1 [dict get $info poss1]
			set cur [list_sub $line $poss1]	
			set newvalues [list_sub $line2 [dict get $info poss2]]
			foreach c $cur pos $poss1 newvalue $newvalues {
				if {(!$force) && $c ne "?"} continue
				lset line $pos $newvalue
			}
			if {($seq eq "?" || $zyg eq "?")} {
				if {$samplea(annot_regionfile,$sample)} {
					set r [annot_region_in $samplea(regionfile,$sample) $chr $begin $end]
					if {!$r} {set seq u}
				} else {
					set seq ?
				}
				if {$seq eq "u"} {
					set zyg u
				} elseif {$type ne "snp"} {
					# for simplicity, we indicate for indels reference
					# to be completely correct, we should find out potential 
					# differences in the deletion to annotate o
					# but this is (currently) not done
					set ref [lindex $line $samplea(ref,$sample)]
					if {$a1 == "?" || $a2 == "?"} {
						set a1 $ref
						lset line $samplea(a1,$sample) $ref
						lset line $samplea(a2,$sample) $ref
						set seq r ; set zyg r
					} else {
						set alt [lindex $line $samplea(alt,$sample)]
						set zyg [zyg $a1 $a2 $ref $alt]
						if {$zyg in "r o"} {set seq r} else {set seq v}
					}
				} elseif {$samplea(ref,$sample) != -1} {
					set ref [lindex $line $samplea(ref,$sample)]
					set alt [lindex $line $samplea(alt,$sample)]
					set a1 [lindex $line $samplea(a1,$sample)]
					set a2 [lindex $line $samplea(a2,$sample)]
					if {$end <= $dend} {
						# overlapping a deletion
						if {$a1 eq $ref} {
							set a2 @
							lset line $samplea(a2,$sample) @
						} else {
							set a1 @
							lset line $samplea(a1,$sample) @
						}
					}
					set zyg [zyg $a1 $a2 $ref $alt]
					if {$zyg in "r o"} {set seq r} else {set seq v}
				} else {
					set zyg ?
				}
				lset line $samplea(seq,$sample) $seq
				if {$samplea(zyg,$sample) != -1} {lset line $samplea(zyg,$sample) $zyg}
			}
#putsvars line
			break
		} elseif {$d < 0} {
			if {$seq eq "?" || $zyg eq "?"} {
#puts ------d<0
#putsvars loc sample line line2 info
				if {$end < $dend} {
					# in a deletion
					lset line $samplea(a1,$sample) @
					lset line $samplea(a2,$sample) @
					lset line $samplea(seq,$sample) r
					if {$samplea(zyg,$sample) != -1} {lset line $samplea(zyg,$sample) o}
				} else {
					lset line $samplea(seq,$sample) u
					if {$samplea(zyg,$sample) != -1} {lset line $samplea(zyg,$sample) u}
				}
#putsvars line
			}
			break
		}
		if {![llength $line2]} break
		set line2 [annot_varall_getline info]
		set dend [dict get $info dend]
		if {$dend != -1} {
			set dchr [dict get $info dchr]
			if {[loc_compare $dchr $chr] != 0} {
				set dend -1
				dict set info dend -1
			}
		}
	}
	dict set info line2 $line2
	set annot(varallinfo,$varallfile,$sample) $info
}

proc annot_varall_close {varallfile sample} {
	global annot
	set info $annot(varallinfo,$varallfile,$sample)
	catch {close [dict get $info af]}
	unset annot(varallinfo,$varallfile,$sample)
}

