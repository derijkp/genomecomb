proc annot_region_init {regfile} {
	global annot
	catch {annot_region_close $regfile}
	set fr [open $regfile]
	set header [split [gets $fr] \t]
	set poss [list_cor $header {chromosome begin end}]
	if {[lsearch $poss -1] != -1} {
		set poss [list_cor $header {chr patchstart pos}]
	}
	if {[lsearch $poss -1] != -1} {
		error "Wrong header"
	}
	set line [split [gets $fr] \t]
	foreach {chr2 start2 end2} [list_sub $line $poss] break
	set nchr2 [chr2num $chr2]
	set annot(reg,$regfile) [list $fr $nchr2 $start2 $end2 $poss]
}

proc annot_region_get {regfile chr1 start1 end1} {
	global annot
	set nchr1 [chr2num $chr1]
	foreach {fr nchr2 start2 end2 poss2} $annot(reg,$regfile) break
	set cur 0
	while 1 {
		# putsvars chr1 chr2 start1 end1 start2 end2
		if {$nchr2 > $nchr1} break
		if {$nchr2 == $nchr1} {
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
		set nchr2 [chr2num $chr2]
	}
	set annot(reg,$regfile) [list $fr $nchr2 $start2 $end2 $poss2]
	return $cur
}

proc annot_region_in {regfile chr1 start1 end1} {
	global annot
	set nchr1 [chr2num $chr1]
	foreach {fr nchr2 start2 end2 poss2} $annot(reg,$regfile) break
	set cur 0
	while 1 {
		# putsvars chr1 chr2 start1 end1 start2 end2
		if {$nchr2 > $nchr1} break
		if {$nchr2 == $nchr1} {
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
		set nchr2 [chr2num $chr2]
	}
	set annot(reg,$regfile) [list $fr $nchr2 $start2 $end2 $poss2]
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
	set annot(cov,$dir) {-1 {} 0}
}

proc annot_coverage_get {dir chr begin} {
	global annot
	foreach {curchr chrfile present} [get annot(cov,$dir) {{} {} 0}] break
	if {$chr ne $curchr} {
		if {[llength $chrfile]} {
			tsv_index_close $chrfile offset
		}
		set chrfile [lindex [glob -nocomplain $dir/coverage/coverageRefScore-$chr-*-ASM.tsv $dir/coverage/coverageRefScore-$chr-*-ASM.tsv.rz $dir/coverage/coverageRefScore-$chr-*-ASM.tsv.gz] 0]
		if {[llength $chrfile]} {
			tsv_index_open $chrfile offset 1
			set present 1
		} else {
			set present 0
		}
		set annot(cov,$dir) [list $chr $chrfile $present]
	}
	if {!$present} {return {? ?}}
	ifcatch {tsv_index_get $chrfile offset $begin} cline -regexp {
		"not found in" {set cline {? ? ?}}
	}
	foreach {refscore coverage} {{} {}} break
	foreach {offset refscore coverage} $cline break
	return [list $refscore $coverage]
}

proc annot_coverage_close {dir} {
	global annot
	foreach {curchr chrfile present} [get annot(cov,$dir) {{} {} 0}] break
	if {[llength $chrfile]} {
		tsv_index_close $chrfile offset
	}
	unset annot(cov,$dir)
}
