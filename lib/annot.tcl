proc annot_region_init {regfile} {
	global annot
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

proc annot_region_close {regfile} {
	global annot
	foreach {fr} $annot(reg,$regfile) break
	catch {close $fr}
	unset annot(reg,$regfile)
}

proc annot_coverage_init {dir} {
	global annot
	set annot(cov,$dir) {-1 {} 0 {}}
}

proc annot_coverage_get {dir chr begin} {
	global annot
	foreach {curchr fc rem chrfile} $annot(cov,$dir) break
	if {$chr ne $curchr} {
		catch {close $fc}
		if {$rem && [file extension $chrfile] ne ".gz"} {file delete $chrfile}
		set chrfile [glob -nocomplain $dir/ASM/REF/coverageRefScore-$chr-*-ASM.tsv]
		if {![llength $chrfile]} {
			set chrfile [glob -nocomplain $dir/ASM/REF/coverageRefScore-$chr-*-ASM.tsv.gz]
			puts "temporarily uncompressing $chrfile"
			exec gunzip -c $chrfile > [file root $chrfile]
			set chrfile [glob -nocomplain $dir/ASM/REF/coverageRefScore-$chr-*-ASM.tsv]
			lset annot(cov,$dir) 2 1
		} else {
			lset annot(cov,$dir) 2 0
		}
		set chrfile [lindex $chrfile 0]
		lset annot(cov,$dir) 3 $chrfile
		set fc [open $chrfile]
		lset annot(cov,$dir) 1 $fc
		while {![eof $fc]} {
			set cline [split [gets $fc] \t]
			if {[string index $cline 0] eq ">"} break
		}
		lset annot(cov,$dir) 0 $chr
	}
	set cline [tsv_nextline $fc 0 $begin]
	foreach {refscore coverage} {{} {}} break
	foreach {offset refscore coverage} $cline break
	return [list $refscore $coverage]
}

proc annot_coverage_close {dir} {
	global annot
	foreach {curchr fc rem chrfile} $annot(cov,$dir) break
	catch {close $fc}
	unset annot(cov,$dir)
}
