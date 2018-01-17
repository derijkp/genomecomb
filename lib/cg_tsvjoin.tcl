proc cg_tsvjoin {args} {
	set pre1 {}
	set pre2 {}
	set idfields1 {}
	set idfields2 {}
	set type f
	cg_options tsvjoin args {
		-pre1 {set pre1 $value}
		-pre2 {set pre2 $value}
		-idfields - -idfields1 {set idfields1 $value}
		-idfields2 {set idfields2 $value}
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
	} {file1 file2 resultfile} 2 3
	if {![llength $idfields2]} {set idfields2 $idfields1}
	if {[info exists resultfile]} {
		set o [open $resultfile w]
	} else {
		set o stdout
	}

	catch {close $f1} ; catch {close $f2}
	set f1 [gzopen $file1]
	set header1 [tsv_open $f1]
	set f2 [gzopen $file2]
	set header2 [tsv_open $f2]
	if {![llength $idfields1]} {
		set idfields1 [list_common $header1 $header2]
		set idfields2 $idfields1
	}
	set poss1 [list_cor $header1 $idfields1]
	if {[inlist $poss1 -1]} {error "missing idfields in $file1: [list_sub $idfields1 [list_find $poss1 -1]]"}
	set poss2 [list_cor $header2 $idfields2]
	if {[inlist $poss2 -1]} {error "missing idfields in $file2: [list_sub $idfields2 [list_find $poss2 -1]]"}
	set remain2 [list_lremove $header2 $idfields2]
	set defline1 [join [list_fill [llength $header1] {}] \t]
	set defline2 [join [list_fill [llength $remain2] {}] \t]
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
			lappend newheader ${pre1}$field
		}
	}
	if {$pre2 eq ""} {
		lappend newheader {*}$remain2
	} else {
		set newheader {}
		foreach field $remain2 {
			lappend newheader ${pre2}$field
		}
	}
	puts $o [join $newheader \t]
	set previd1 {}
	set line2 [split [gets $f2] \t]
	set id2 [list_sub $line2 $poss2]
	set previd2 $id2
	while 1 {
		if {[gets $f1 line1] == -1} {
			if {$type in "f r"} {
				while {[gets $f2 line2] != -1} {
					set line2 [split $line2 \t]
					puts $o $defline1\t[join [list_sub $line2 -exclude $poss2] \t]
				}
			}
			break
		}
		set sline1 [split $line1 \t]
		if {![llength $sline1]} continue
		set id1 [list_sub $sline1 $poss1]
		if {[loc_compare $previd1 $id1] > 0} {
			error "file $file1 not properly sorted on $idfields1 ($previd1 comes before $id1)"
		}
		set previd1 $id1
		set c [loc_compare $id1 $id2]
		while {$c > 0} {
			if {$type in "f r"} {
				puts $o $defline1\t[join [list_sub $line2 -exclude $poss2] \t]
			}
			set line2 [split [gets $f2] \t]
			set id2 [list_sub $line2 $poss2]
			if {[loc_compare $previd2 $id2] > 0} {
				error "file $file2 not properly sorted on $idfields2 ($previd2 comes before $id2)"
			}
			set previd2 $id2
			set c [loc_compare $id1 $id2]
		}
		if {$c == 0} {
			puts $o [join $line1 \t]\t[join [list_sub $line2 -exclude $poss2] \t]
			set line2 [split [gets $f2] \t]
			set id2 [list_sub $line2 $poss2]
			if {[loc_compare $previd2 $id2] > 0} {
				error "file $file2 not properly sorted on $idfields2 ($previd2 comes before $id2)"
			}
			set previd2 $id2
		} else {
			if {$type in "f l"} {
				puts $o [join $line1 \t]\t$defline2
			}
		}
	}

	gzclose $f1
	gzclose $f2
	if {$o ne "stdout"} {close $o}
}
