if 0 {

package require Tclx
signal -restart error SIGINT
lappend auto_path ~/dev/completegenomics/lib
package require Extral
cd /complgen/sv
	cd /complgen/sv
	set outfile svcompar_test.tsv
	set outfile svcompar_GS102_GS103/svcompar_GS102_GS103-20.tsv
	set svfile1 sv79-20-pairs-sv.tsv
	set svfile2 sv78-20-pairs-sv.tsv
	set outfile svcompar_sv78_sv79-20.tsv
	set svfile1 GS103/GS103-9-paired-sv.tsv
	set svfile2 oldGS103-9-paired-sv.tsv
	set svfile1 GS103/GS103-9-paired-sv.tsv
	set svfile2 GS102/GS102-9-paired-sv.tsv
	set svfile1 GS103/GS103-20-paired-sv.tsv
	set svfile2 GS102/GS102-20-paired-sv.tsv

}

proc svmulticompar_groupdists {dlist} {
	set curdist [lindex $dlist 0 5]
	set result {}
	set newlist {}
	foreach line $dlist {
		set d [lindex $line 5]
		set diff [expr {abs($d-$curdist)}]
		set sizediff [min [max [expr {round(0.1*$curdist)}] 35] 2000]
		if {$diff > $sizediff} {
			lappend result $newlist
			set curdist [lindex $line 5]
			set newlist [list $line]
		} else {
			lappend newlist $line
		}
	}
	if {[llength $newlist]} {
		lappend result $newlist
	}
	return $result
}

proc svmulticompar_groupcompatible {plist} {
	set plist [lsort -integer -index 3 $plist]
	set curend1 [lindex $plist 0 3]
	set curstart2 [lindex $plist 0 9]
	set result {}
	set newlist {}
	foreach line $plist {
		set end1 [lindex $line 3]
		set start2 [lindex $line 9]
		if {[overlap $end1 $start2 $curend1 $curstart2] > 1} {
			lappend newlist $line
		} else {
			lappend result $newlist
			set newlist {}
			set curend1 $end1
			set curstart2 $start2
		}
	}
	if {[llength $newlist]} {
		lappend result $newlist
	}
	return $result
}

proc svmulticompar_compar {line1 line2} {
	if {![llength $line1]} {return -1}
	if {![llength $line2]} {return 1}
	foreach {src chr1 start1} $line1 break
	foreach {src chr2 start2} $line2 break
	if {$chr1 ne $chr2} {
		set nchr1 [chr2num $chr1]
		set nchr2 [chr2num $chr2]
		return [expr {$nchr2 - $nchr1}]
	}
	set end1 [lindex $line1 9]
	set end2 [lindex $line2 9]
	if {$end2 < $start1} {
		return -1
	} elseif {$end1 < $start2} {
		return 1
	} else {
		return 0
	}
}

proc svmulticompar_getline {f poss {type 1}} {
	while 1 {
		set line [split [gets $f] \t]
		if {[llength $line]} break
		if {[eof $f]} {return {}}
	}
	set cur [list_sub $line $poss]
	if {[lindex $cur 3] eq "trans"} {
		set endpos [expr {[lindex $cur 2]+200}]
		lset cur 4 0
		lset cur 7 $endpos
		lset cur 8 [expr {$endpos+400}]
	}
	list_concat $type $cur $line
}

proc svmulticompar_write {o id group poss2 dummy1 dummy2 {ddummy1 {}} {ddummy2 {}}} {
#		set start1 [lmath_max [list_subindex $group 2]]
#		set end1 [lmath_max [list_subindex $group 3]]
#		set size [expr {round([lmath_average [list_subindex $group 5]])}]
#		set start2 [lmath_max [list_subindex $group 8]]
#		set end2 [lmath_max [list_subindex $group 9]]
	unset -nocomplain todo
	foreach line $group {
		lappend todo([lindex $line 0]) $line
	}
	if {![info exists todo(1)] || ![info exists todo(2)]} {
		set udummy1 $dummy1; set udummy2 $dummy2
	} else {
		set udummy1 $ddummy1; set udummy2 $ddummy2
	}
	set max [max [llength [get todo(1) ""]] [llength [get todo(2) ""]]]
	if {$max > 1} {append id -$max}
	foreach l1 [get todo(1) ""] l2 [get todo(2) ""] {
		if {[llength $l2]} {
			set oline2 [lrange $l2 10 end]
			set cur2 [list_sub $oline2 $poss2]
		} else {
			set oline2 $udummy2
		}
		if {[llength $l1]} {
			set oline1 [lrange $l1 10 end]
		} else {
			set oline1 $udummy1
			set oline1 [lreplace $oline1 1 9 {*}$cur2]
		}
		lset oline1 0 $id
		puts $o [join [list_concat $oline1 $oline2] \t]
	}
}

proc svmulticompar_getlist {f1 poss1 len1 line1Var f2 poss2 len2 line2Var} {
	upvar $line1Var line1
	upvar $line2Var line2
	set list {}
	set compar [svmulticompar_compar $line1 $line2]
	if {$compar >= 0} {
		set listchr [lindex $line1 1]	
		set liststart [lindex $line1 2]
		set listend [lindex $line1 9]
	} elseif {$compar < 0} {
		set listchr [lindex $line2 1]	
		set liststart [lindex $line2 2]
		set listend [lindex $line2 9]
	}
	while {![eof $f1] || ![eof $f2]} {
		set match 0
		if {[llength $line1]} {
			set listchr1 [lindex $line1 1]
			set liststart1 [lindex $line1 2]
			if {($listchr1 == $listchr) && ($liststart1 < $listend)} {
				lappend list $line1
				set listend [lindex $line1 9]
				set line1 [svmulticompar_getline $f1 $poss1 1]
				set match 1
			}
		}
		if {[llength $line2]} {
			set listchr2 [lindex $line2 1]
			set liststart2 [lindex $line2 2]
			if {($listchr2 == $listchr) && ($liststart2 < $listend)} {
				lappend list $line2
				set listend [lindex $line2 9]
				set line2 [svmulticompar_getline $f2 $poss2 2]
				set match 2
			}
		}
		if {!$match} break
	}
	return $list
}

proc svmulticompar {svfile1 svfile2} {

	set locfields {chr1 start1 end1 type size zyg chr2 start2 end2}
	if {![file exists $svfile1]} {
		set o [open $svfile1 w]
		puts $o id\t[join $locfields \t]
		close $o
	}

	catch {close $f1}; catch {close $f2}; catch {close $o}; catch {file delete $tempfile2}
	set o [open $svfile1.temp w]
	#
	# open compar file
	set f1 [open $svfile1]
	set header1 [split [gets $f1] \t]
	set len1 [llength $header1]
	set poss1 [list_cor $header1 $locfields]
	set dummy1 [list_fill [llength $header1] {}]
	set ddummy1 [list_fill [llength $header1] d]
	#
	# open add file
	set name [file root [file tail $svfile2]]
	set name [lindex [split $name -] end]
	set tempfile2 [tempfile]
	cg select -s "chr1 end1" < $svfile2 > $tempfile2
	set f2 [open $tempfile2]
	set header2 [split [gets $f2] \t]
	set len2 [llength $header2]
	set poss2 [list_cor $header2 $locfields]
	set dummy2 [list_fill [llength $header2] {}]
	set ddummy2 [list_fill [llength $header2] d]
	#
	# make new header
	set header $header1
	foreach field $header2 {append header \t${field}-$name}
	puts $o [join $header \t]
	#
	# go over files
	set line1 [svmulticompar_getline $f1 $poss1 1]
	set line2 [svmulticompar_getline $f2 $poss2 2]

	set did 1
	while {![eof $f1] || ![eof $f2]} {
		set list [svmulticompar_getlist $f1 $poss1 $len1 line1 $f2 $poss2 $len2 line2]
#if {[lindex $list 0 4] eq "trans"} {error STOP}
#puts [join $list \n]\n\n
# join $list \n\n
		if {[llength $list] == 1} {
			svmulticompar_write $o $did $list $poss2 $dummy1 $dummy2
			incr did
			continue
		}
		# split on type
		unset -nocomplain todo
		foreach line [lsort -integer -index 5 $list] {
			set type [lindex $line 4]
			set size [lindex $line 5]
			lappend todo($type) $line
		}
		set list {}
		foreach type [array names todo] {
			if {[llength $todo($type)] < 2} {
				lappend list $todo($type)
				continue
			}
			set plists [svmulticompar_groupdists $todo($type)]
			foreach plist $plists {
				set clists [svmulticompar_groupdists $plist]
				foreach clist $clists {
					lappend list $clist
				}
			}
		}
		foreach group $list {
			svmulticompar_write $o $did $group $poss2 $dummy1 $dummy2 $ddummy1 $ddummy2
			incr did
		}		


	}


	flush $o
	close $o
	close $f1
	close $f2
	file delete $tempfile2
	file rename -force $svfile1 $svfile1.old
	cg select -s {chr1 end1} $svfile1.temp $svfile1
	putslog "finished adding $name to $svfile1"

}
