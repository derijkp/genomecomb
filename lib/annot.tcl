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

proc annot_coverage_init {dir} {
	global annot
	catch {annot_coverage_close $dir}
	set annot(cov,$dir) {-1 {} 0 {} {} {}}
}

proc annot_coverage_get {dir chr begin} {
	global annot
	set present 1
	foreach {curchr chrfile present poss type obj} [get annot(cov,$dir) {{} {} 0}] break
	if {$chr ne $curchr} {
		switch $type {
			{} {
				tsv_index_close $chrfile offset
			}
			bcol {
				foreach {refbcol covbcol} $obj break
				bcol_close $refbcol ; bcol_close $covbcol
			}
		}
		set nchr [chr_clip $chr]
		set chrfile [gzfile $dir/coverage/coverage-*-$nchr.bcol $dir/coverage/coverage-*-chr$nchr.bcol $dir/coverage/coverage-$nchr-*.bcol $dir/coverage/coverage-chr$nchr-*.bcol]
		if {[file exists $chrfile]} {
			set reffile [gzfile $dir/coverage/refScore-*-$nchr.bcol $dir/coverage/refScore-*-chr$nchr.bcol $dir/coverage/refScore-$nchr-*.bcol $dir/coverage/refScore-chr$nchr-*.bcol]
			set type bcol
		} else {
			set chrfile [gzfile $dir/coverage/coverageRefScore-$nchr-*.tsv $dir/coverage/coverage-$nchr-*.tsv $dir/coverage/coverageRefScore-chr$nchr-*.tsv $dir/coverage/coverage-chr$nchr-*.tsv]
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
		if {!$present} {return {u u}}
		if {$type eq "bcol"} {
			if {[file exists $reffile]} {
				set refbcol [bcol_open $reffile]
			} else {
				set refbcol {}
			}
			set covbcol [bcol_open $chrfile]
			set obj [list $refbcol $covbcol]
			set annot(cov,$dir) [list $chr $chrfile $present $poss $type $obj]
		} else {
			tsv_index_open $chrfile offset 1
			set header [tsv_index_header $chrfile]
			set rpos [lsearch $header refScore]
			set cpos [lsearch $header uniqueSequenceCoverage]
			if {$cpos == -1} {set cpos [lsearch $header coverage]}
			set poss [list $rpos $cpos]
			set annot(cov,$dir) [list $chr $chrfile $present $poss $type {}]
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

proc annot_coverage_close {dir} {
	global annot
	foreach {curchr chrfile present poss type obj} [get annot(cov,$dir) {{} {} 0 {} {} {}}] break
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
	unset annot(cov,$dir)
}
