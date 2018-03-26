#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc svopen {file side {name {}}} {
	if {$name eq ""} {
		set name [file root [file tail $file]]
		regsub {^[^-]*-} $name {} name
	}
	set fh [gzopen $file]
	set header [tsv_open $fh]
	set poss [tsv_basicfields $header 6 0]
	set pos [lsearch $poss -1]
	if {$pos != -1 && $pos < 4} {error "One of chromosome, begin, end, type (or equivalent) missing in file $file"}
	set temp [list_fill [llength $header] 0 1]
	set newheader {chromosome begin end type ref alt}
	foreach field [list_sub $header [list_lremove $temp $poss]] {
		if {[string first - $field] == -1} {append field -$name}
		lappend newheader $field
	}
	set poss [list -1 {*}$poss {*}[list_lremove $temp $poss]]
	if {[lindex $poss 5] == -1} {
		set makeref 1
	} else {
		set makeref -1
	}
	if {[lindex $poss 6] == -1} {
		set makealt [lsearch $header size]
		# if {$makealt == -1} {error "no alt and no size column for file $file"}
		set makealt [lsearch $poss $makealt]
		lappend makealt [lsearch $poss [lsearch $header chr2]]
		lappend makealt [lsearch $poss [lsearch $header start2]]
	} else {
		set makealt ""
	}
	return [list $fh $poss $side $makeref $makealt $newheader $header]
}

proc svclose {f} {
	foreach {fh} $f break
	gzclose $fh
}

proc svgetline {f} {
	foreach {fh poss side makeref makealt} $f break
	while 1 {
		if {[gets $fh line] == -1} {
			return {}
		}
		set line [split $line \t]
		if {[llength $line]} break
	}
	set line [list_sub $line $poss]
	lset line 0 $side
	foreach {side chr begin end type} $line break
	if {$makeref != -1 || [llength $makealt]} {
		lset line 1 [chr_clip $chr]
		if {$makeref != -1} {
			lset line 5 [expr {$end - $begin}]
		}
		if {[llength $makealt]} {
			if {$type eq "inv"} {
				lset line 6 i
			} elseif {$type eq "del"} {
				lset line 6 ""
			} elseif {$type eq "ins"} {
				lset line 6 [lindex $line [lindex $makealt 0]]
			} elseif {$type eq "trans"} {
				foreach {size chr2 start2} [list_sub $line $makealt] break
				lset line 6 \[[chr_clip $chr2]:$start2\[
			} else {
				lset line 6 ?
			}
		}
	} else {
		lset line 1 [chr_clip $chr]
		if {$type eq "trans"} {
			set temp [lindex $line 6]
			regsub {([\[\]])chr} $temp {\1} temp
			lset line 6 $temp
		}
	}
	return $line
}

proc svmulticompar_dist {sline1 sline2 {margin 30} {lmargin 300} {overlap 75}} {
	foreach {side1 chr1 begin1 end1 type1 ref1 alt1} $sline1 break
	foreach {side2 chr2 begin2 end2 type2 ref2 alt2} $sline2 break
	if {$type1 in "del inv"} {
		set overlap1 [max $begin1 $begin2]
		set overlap2 [min $end1 $end2]
		set psize [expr {100*($overlap2 - $overlap1)}]
		if {$end1 == $begin1 || [expr {$psize/($end1-$begin1)}] < $overlap} {return 2147483648}
		if {$end2 == $begin2 || [expr {$psize/($end2-$begin2)}] < $overlap} {return 2147483648}
		set margin $lmargin
	} elseif {$type1 eq "trans"
		&& [regexp {\[?([^:]+):([0-9]+)\]?} $alt1 temp tchr1 tbegin1]
		&& [regexp {\[?([^:]+):([0-9]+)\]?} $alt2 temp tchr2 tbegin2]
	} {
		if {$tchr1 ne $tchr2} {return 2147483648}
		set diff [expr {abs($tbegin2 - $tbegin1)}]
		if {$diff > $lmargin} {return 2147483648}
		return [expr {abs($begin2 - $begin1) + $diff}]
	}
	set enddiff [expr {abs($end2 - $end1)}]
	if {$enddiff > $margin} {return 2147483648}
	if {[isint $alt1]} {
		if {[isint $alt2]} {
			set altdiff [expr {abs($alt2 - $alt1)}]
			if {$altdiff > $margin} {return 2147483648}
		} else {
			set alt2 [string length $alt2]
			set altdiff [expr {abs($alt2 - $alt1)}]
			if {$altdiff > $margin} {return 2147483648}
		}
	} elseif {[isint $alt2]} {
		set alt1 [string length $alt1]
		set altdiff [expr {abs($alt2 - $alt1)}]
		if {$altdiff > $margin} {return 2147483648}
	} else {
		if {$alt1 ne $alt2} {return 2147483648}
		set altdiff 0
	}
	set diff [expr {abs($begin2 - $begin1) + $enddiff + $altdiff}]
	return $diff
}

proc svmulticompar_out {line1 line2 dummy1 dummy2} {
	if {[llength $line1]} {
		if {[llength $line2]} {
			list {*}[lrange $line1 1 end] {*}[list_sub $line2 {2 3 6}] {*}[lrange $line2 7 end]
		} else {
			list {*}[lrange $line1 1 end] {*}$dummy2
		}
	} elseif {[llength $line2]} {
		list {*}[lrange $line2 1 6] {*}$dummy1 {*}[list_sub $line2 {2 3 6}] {*}[lrange $line2 7 end]
	} else {
		error "both line1 and line2 empty"
	}
}

proc svmulticompar_groupdists {list1 list2 dummy1 dummy2 {margin 30} {lmargin 300} {overlap 75}} {
# putsvars list1 list2
# if {[llength $list1] > 1 && [llength $list2] > 1} {error STOP}
	set matches {}
	set p1 0
	foreach sline1 $list1 {
		set p2 0
		foreach sline2 $list2 {
			set diff [svmulticompar_dist $sline1 $sline2 $margin $lmargin $overlap]
			if {$diff < 2147483640} {
				lappend matches [list $diff $p1 $p2]
			}
			incr p2
		}
		incr p1
	}
	set matches [lsort -index 0 -integer $matches]
	unset -nocomplain done1
	unset -nocomplain done2
	set pos 0
	list_foreach {diff p1 p2} $matches {
		incr pos
		if {[info exists done2($p2)]} continue
		set done2($p2) 1
		if {[info exists done1($p1)]} {
			if {[lsearch [list_subindex [lrange $matches $pos end] 1] $p1] != -1} continue
		} else {
			lappend done1($p1) $p2
		}
	}
	set result {}
	unset -nocomplain done2
	set p1s [array names done1]
	foreach p1 $p1s {
		set line1 [lindex $list1 $p1]
		foreach p2 $done1($p1) {
			lappend result [svmulticompar_out $line1 [lindex $list2 $p2] $dummy1 $dummy2]
			set done2($p2) 1
		}
	}
	foreach line1 [list_sub $list1 -exclude $p1s] {
		lappend result [svmulticompar_out $line1 {} $dummy1 $dummy2]
	}
	foreach line2 [list_sub $list2 -exclude [array names done2]] {
		lappend result [svmulticompar_out {} $line2 $dummy1 $dummy2]
	}
	return $result
}

proc svmulticompar_getline {f poss {type 1}} {
	global cchr cpos locpos
	# -1 because count has not been prepended yet
	set chrpos [expr {$locpos(chromosome)-1}]
	set typepos [expr {$locpos(type)-1}]
	set sizepos [expr {$locpos(size)-1}]
	set end2pos [expr {$locpos(end2)-1}]
	while 1 {
		set line [split [gets $f] \t]
		if {[llength $line]} break
		if {[eof $f]} {return {}}
	}
	set cur [list_sub $line $poss]
	if {[lindex $cur $typepos] eq "trans"} {
		set end1pos [expr {$locpos(end1)-1}]
		set start2pos [expr {$locpos(start2)-1}]
		set endpos [expr {[lindex $cur $end1pos]+200}]
		lset cur $sizepos 0
		lset cur $start2pos $endpos
		lset cur $end2pos [expr {$endpos+400}]
	}
	if {![isint [lindex $cur $sizepos]]} {
		set beginpos [expr {$locpos(begin)-1}]
		set endpos [expr {$locpos(end)-1}]
		lset cur $sizepos [expr {[lindex $cur $endpos]-[lindex $cur $beginpos]}]
	}
	set temp [lindex $cur $chrpos]
	if {$temp ne $cchr} {
		set cchr $temp
		putslog "Starting chromosome $cchr"
		set cpos 1000000
	}
	if {[lindex $cur $end2pos] > $cpos} {
		putslog $cpos
		incr cpos 1000000
	}
	list_concat $type $cur $line
}

proc svmulticompar_getlist {f1 line1Var f2 line2Var maxmargin} {
	upvar $line1Var line1
	upvar $line2Var line2
	global locpos
	# lines have the following format: chromosome begin end type ref alt ...
# puts [join [list_subindex [list $line1 $line2] {0 1 2 3 4 5}] \n]
	set return 0
	if {![llength $line1]} {
		if {![llength $line2]} {return {}}
		set return 2
	} elseif {![llength $line2]} {
		set return 1
	} else {
		foreach {side1 chr1 begin1 end1} $line1 break
		foreach {side2 chr2 begin2 end2} $line2 break
		set chrcomp [loc_compare $chr1 $chr2]
		if {$chrcomp < 0} {
			set return 1
		} elseif {$chrcomp > 0} {
			set return 2
		} elseif {[expr {$begin1+$maxmargin}] <= $begin2} {
			set return 1
		} elseif {[expr {$begin2+$maxmargin}] <= $begin1} {
			set return 2
		}
	}
	if {$return == 1} {
		lappend linea 1
		set list [list $line1]
		set line1 [svgetline $f1]
		return $list
	} elseif {$return == 2} {
		lappend line2 2
		set list [list $line2]
		set line2 [svgetline $f2]
		return $list
	}
	set curchr $chr1
	if {$begin1 < $begin2} {
		set curstart [expr {$begin1 - $maxmargin}]
		set curend [expr {$begin2 + $maxmargin}]
	} else {
		set curstart [expr {$begin2 - $maxmargin}]
		set curend [expr {$begin1 + $maxmargin}]
	}
	set list [list $line1 $line2]
	set line1 [svgetline $f1]
	foreach {side1 chr1 begin1} $line1 break
	set line2 [svgetline $f2]
	foreach {side2 chr2 begin2} $line2 break
	while {[llength $line1] || [llength  $line2]} {
		set match 0
		if {[llength $line1]} {
			if {$chr1 == $curchr && $begin1 < $curend} {
				lappend list $line1
				set temp [expr {$begin1+$maxmargin}]
				if {$temp > $curend} {
					set curend $temp
				}
				set line1 [svgetline $f1]
				foreach {side1 chr1 begin1} $line1 break
				set match 1
			}
		}
		if {[llength $line2]} {
			if {$chr2 == $curchr && $begin2 < $curend} {
				lappend list $line2
				set temp [expr {$begin2+$maxmargin}]
				if {$temp > $curend} {
					set curend $temp
				}
				set line2 [svgetline $f2]
				foreach {side2 chr2 begin2} $line2 break
				set match 2
			}
		}
		if {!$match} break
	}
	return $list
}

proc svmulticompar {args} {
	global locpos
	set margin 30
	set lmargin 300
	set overlap 75
	cg_options svmulticompar args {
		-margin {set margin $value}
		-lmargin {set lmargin $value}
		-overlap {set overlap $value}
	} {svfile1 svfile2}
	set maxmargin [max $margin $lmargin]

	set locfields {chromosome begin end type ref alt}
	if {![file exists $svfile1]} {
		set name [file root [file tail $svfile2]]
		regsub {^[^-]*-} $name {} name
		catch {svclose $f} ; catch {close $o}
		set f [svopen $svfile2 1]
		foreach {fh side poss makeref makealt newheader} $f break
		set o [open $svfile1 w]
		puts $o [join [linsert $newheader 6 lbegin-$name lend-$name lalt-$name] \t]
		while 1 {
			set line [svgetline $f]
			if {![llength $line]} break
			puts $o [join [list {*}[lrange $line 1 6] {*}[list_sub $line {2 3 6}] {*}[lrange $line 7 end]] \t]
		}
		close $o
		svclose $f
		return
	}

	#
	# open compar file
	catch {svclose $f1}
	set f1 [svopen $svfile1 1]
	foreach {fh1 side1 poss1 makeref1 makealt1 header1} $f1 break
	if {[lrange $header1 0 5] ne {chromosome begin end type ref alt}} {
		error "error in compar file $svfile1, must start with fields: chromosome begin end type ref alt"
	}
	set len1 [llength $header1]
	set dummy1 [list_fill [expr {[llength $header1]-6}] {}]
	set ddummy1 [list_fill [expr {[llength $header1]-6}] d]
	set line1 [svgetline $f1]
	#
	# open add file
	set name [file root [file tail $svfile2]]
	regsub {^[^-]*-} $name {} name
	catch {svclose $f2}
	set f2 [svopen $svfile2 2 $name]
	foreach {fh2 side2 poss2 makeref2 makealt2 header2} $f2 break
	set len2 [llength $header2]
	set dummy2 [list_fill [expr {[llength $header2]-3}] {}]
	set ddummy2 [list_fill [expr {[llength $header2]-3}] d]
	set line2 [svgetline $f2]

	catch {close $o} ; set o [open $svfile1.temp w]
	# make new header
	set header [list {*}$header1 lbegin-$name lend-$name lalt-$name {*}[lrange $header2 6 end]]
	puts $o [join $header \t]
	#
	# go over files
	set ::cchr {}
	set ::cpos 0
	while {[llength $line1] || [llength $line2]} {
		# get overlapping lines from both files in the following format (locfields):
		# {src chr1 begin end type start1 end1 size zyg chr2 start2 end2 (complete line)}
		set list [svmulticompar_getlist $f1 line1 $f2 line2 $maxmargin]
#puts ------------
#putsvars list line1 line2
#puts [join [list_subindex $list {0 1 2 3 4 5}] \n]
#if {[lsearch [list_subindex $list 3] 2897835] != -1} {error STOPPED}
#if {[lsearch [list_subindex $list 3] 67759423] != -1} {error STOPPED}
#if {[lindex $list 0 4] eq "trans"} {error STOP}
#puts [join $list \n]\n\n
# join $list \n\n
		if {[llength $list] == 1 || [llength [list_remdup [list_subindex $list 0]]] == 1} {
			# from single source
			foreach line $list {
				set side [lindex $line 0]
				if {$side == 1} {
					puts $o [join [svmulticompar_out $line {} $dummy1 $dummy2] \t]
				} else {
					puts $o [join [svmulticompar_out {} $line $dummy1 $dummy2] \t]
				}
			}
			continue
		}
		#
		# multi source, match first
		# split on type
		unset -nocomplain todo
		unset -nocomplain typesa
		foreach line $list {
			set side [lindex $line 0]
			set type [lindex $line 4]
			set typesa($type) 1
			lappend todo($type,$side) $line
		}
		set list {}
		foreach type [lsort [array names typesa]] {
			if {![info exists todo($type,1)]} {
				foreach l2 $todo($type,2) {
					lappend list [svmulticompar_out {} $l2 $dummy1 $dummy2]
				}
			} elseif {![info exists todo($type,2)]} {
				foreach l1 $todo($type,1) {
					lappend list [svmulticompar_out $l1 {} $dummy1 $dummy2]
				}
			} else {
				# foreach {list1 list2} [list $todo($type,1) $todo($type,2)] break
				set plist [svmulticompar_groupdists $todo($type,1) $todo($type,2) $dummy1 $dummy2 $margin $lmargin $overlap]
				lappend list {*}$plist
			}
		}
		if {[llength $list] > 1} {
			set list [ssort -natural $list]
		}
		foreach temp $list {
			puts $o [join $temp \t]
		}
	}

	flush $o
	close $o
	svclose $f1
	svclose $f2
	# file delete $tempfile2
	# cg select -s {chr1 start1} $svfile1.temp $svfile1.temp2
	# file delete $svfile1.temp
	file rename -force $svfile1.temp $svfile1
	putslog "finished adding $name to $svfile1"

}

proc cg_svmulticompar {args} {
	set margin 30
	set lmargin 300
	set overlap 75
	cg_options svmulticompar args {
		-margin {set margin $value}
		-lmargin {set lmargin $value}
		-overlap {set overlap $value}
	} {compar_file svfile} 2 ...
	set files [list $svfile {*}$args]
	set done {}
	if {[file exists $compar_file]} {
		set list [cg select -h $compar_file]
		set poss [list_find -glob $list start1-*]
		set done [list_sub $list $poss]
		set done [list_regsub -all {^start1-} $done {}]
	}
	foreach file $files {
		set name [lindex [split [file tail [file root $file]] -] end]
		if {[inlist $done $name]} {
			putslog "Skipping $file: $name already present"
		} else {
			putslog "Adding $file"
			svmulticompar -margin $margin -lmargin $lmargin -overlap $overlap $compar_file $file
		}
	}
}
