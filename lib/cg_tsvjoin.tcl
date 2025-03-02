proc cg_tsvjoin_fullout {o sdefline1 poss1 id2 line2 poss2} {
	set temp $sdefline1
	foreach v $id2 pos $poss1 {
		lset temp $pos $v
	}
	puts $o [join $temp \t]\t[join [list_sub $line2 -exclude $poss2] \t]
}

proc cg_tsvjoin_next {file f idfields poss previdVar prevlineVar idVar lineVar} {
	upvar $previdVar previd
	upvar $prevlineVar prevline
	upvar $idVar id
	upvar $lineVar line
	set previd $id
	set prevline $line
	set r [gets $f line]
	if {$r < -1} {return $r}
	set line [split $line \t]
	set id [list_sub $line $poss]
	if {[nat_compare $previd $id] > 0} {
		error "file $file not properly sorted on $idfields ($previd comes before $id)"
	}
	return $r
}

proc cg_tsvjoin {args} {
	set pre1 {}
	set pre2 {}
	set idfields1 {}
	set idfields2 {}
	set type f
	set usecomments 1
	set sorted 0
	set default {}
	cg_options tsvjoin args {
		-pre1 {set pre1 $value}
		-pre2 {set pre2 $value}
		-idfields - -idfields1 {set idfields1 $value}
		-idfields2 {set idfields2 $value}
		-comments {
			if {$value ni "0 1 2"} {
				error "-comments should be one of:0,1,2"
			}
			set usecomments $value
		}
		-type {
			switch $value {
				full {set type f}
				inner {set type i}
				left {set type l}
				right {set type r}
				default {
					error "unknown join type, must be one of: full inner left right"
				}
			}
		}
		-default {
			set default $value
		}
		-sorted {
			set sorted $value
		}
	} {file1 file2 resultfile} 2 3
	if {![llength $idfields2]} {set idfields2 $idfields1}
	catch {close $f1} ; catch {close $f2}
	set f1 [gzopen $file1]
	set header1 [tsv_open $f1 comments1]
	set f2 [gzopen $file2]
	set header2 [tsv_open $f2 comments2]
	if {![llength $idfields1]} {
		set idfields1 [list_common $header1 $header2]
		set idfields2 $idfields1
	}
	set poss1 [list_cor $header1 $idfields1]
	if {[inlist $poss1 -1]} {error "missing idfields in $file1: [list_sub $idfields1 [list_find $poss1 -1]]"}
	set poss2 [list_cor $header2 $idfields2]
	if {[inlist $poss2 -1]} {error "missing idfields in $file2: [list_sub $idfields2 [list_find $poss2 -1]]"}
	set remain2 [list_lremove $header2 $idfields2]
	set sdefline1 [list_fill [llength $header1] $default]
	set defline1 [join $sdefline1 \t]
	set sdefline2 [list_fill [llength $remain2] $default]
	set defline2 [join $sdefline2 \t]
	if {[llength [list_common $remain2 [list_lremove $header1 $idfields1]]]} {
		if {$pre1 eq $pre2} {
			error "The same (non-id) field is present in both filesnd pre1 and pre2 are the same"
		}
	}
	if {$pre1 eq ""} {
		set newheader $header1
	} else {
		set newheader {}
		foreach field $header1 {
			if {$field in $idfields1} {
				lappend newheader $field
			} else {
				lappend newheader ${pre1}$field
			}
		}
	}
	if {$pre2 eq ""} {
		lappend newheader {*}$remain2
	} else {
		foreach field $remain2 {
			if {$field in $idfields2} {
				lappend newheader $field
			} else {
				lappend newheader ${pre2}$field
			}
		}
	}
	if {!$sorted} {
		gzclose $f1 ; gzclose $f2
		set tempfile1 [tempfile].zst
		set tempfile2 [tempfile].zst
		cg select -s $idfields1 $file1 $tempfile1
		cg select -s $idfields2 $file2 $tempfile2
		set f1 [gzopen $tempfile1]
		set header1 [tsv_open $f1 comments1]
		set f2 [gzopen $tempfile2]
		set header2 [tsv_open $f2 comments2]
	}
	if {[info exists resultfile]} {
		set o [wgzopen $resultfile]
	} else {
		set o stdout
	}
	if {$usecomments eq "1"} {
		puts -nonewline $o $comments1
	} elseif {$usecomments eq "2"} {
		puts -nonewline $o $comments1
	}
	puts $o [join $newheader \t]
	set previd1 {}
	set line1 [split [gets $f1] \t]
	set id1 [list_sub $line1 $poss1]
	set previd2 {}
	set line2 [split [gets $f2] \t]
	set id2 [list_sub $line2 $poss2]
	while 1 {
		# putsvars previd1 prevline1 id1 line1
		# putsvars previd2 prevline2 id2 line2
		set c [loc_compare $id1 $id2]
		if {$c == 0} {
			puts $o [join $line1 \t]\t[join [list_sub $line2 -exclude $poss2] \t]
			set r2 [cg_tsvjoin_next $file2 $f2 $idfields2 $poss2 previd2 prevline2 id2 line2]
			set r1 [cg_tsvjoin_next $file1 $f1 $idfields2 $poss1 previd1 prevline1 id1 line1]
			if {$r1 == -1 || $r2 == -1} break
			while 1 {
				if {[loc_compare $previd1 $id2] == 0} {
					puts $o [join $prevline1 \t]\t[join [list_sub $line2 -exclude $poss2] \t]
					set r2 [cg_tsvjoin_next $file2 $f2 $idfields2 $poss2 previd2 prevline2 id2 line2]
					if {$r2 < -1} break
				} elseif {[loc_compare $id1 $previd2] == 0} {
					puts $o [join $line1 \t]\t[join [list_sub $prevline2 -exclude $poss2] \t]
					set r1 [cg_tsvjoin_next $file1 $f1 $idfields2 $poss1 previd1 prevline1 id1 line1]
					if {$r1 < -1} break
				} else {
					break
				}
			}
		} elseif {$c > 0} {
			if {$type in "f r"} {
				cg_tsvjoin_fullout $o $sdefline1 $poss1 $id2 $line2 $poss2
			}
			set r2 [cg_tsvjoin_next $file2 $f2 $idfields2 $poss2 previd2 prevline2 id2 line2]
			if {$r2 == -1} break
		} else {
			if {$type in "f l"} {
				puts $o [join $line1 \t]\t$defline2
			}
			set r1 [cg_tsvjoin_next $file1 $f1 $idfields2 $poss1 previd1 prevline1 id1 line1]
			if {$r1 == -1} break
		}
	}

	while {$r1 > 0} {
		if {$type in "f l"} {
			puts $o [join $line1 \t]\t$defline2
		}
		set r1 [cg_tsvjoin_next $file1 $f1 $idfields2 $poss1 previd1 prevline1 id1 line1]
	}
	while {$r2 > 0} {
		if {$type in "f r"} {
			cg_tsvjoin_fullout $o $sdefline1 $poss1 $id2 $line2 $poss2
		}
		set r2 [cg_tsvjoin_next $file2 $f2 $idfields2 $poss2 previd2 prevline2 id2 line2]
	}
	gzclose $f1
	gzclose $f2

	if {$o ne "stdout"} {close $o}
}
