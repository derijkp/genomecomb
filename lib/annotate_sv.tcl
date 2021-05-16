#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc annotate_sv_groupdists {list1 list2 emptyout {margin 30} {lmargin 300} {tmargin 300} {overlap 75}} {
# putsvars list1 list2 dummy1 dummy2 margin lmargin tmargin overlap
# if {[llength $list1] > 1 && [llength $list2] > 1} {error STOP}
	set matches {}
	set p1 0
	foreach sline1 $list1 {
		set p2 0
		foreach sline2 $list2 {
			# returns distance between Svs (begindiff + enddiff + altdiff)
			set diff [svmulticompar_dist $sline1 $sline2 $margin $lmargin $tmargin $overlap]
			# diff = 2147483640 if it goes out of the given margins (= no match)
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
		if {[info exists done1($p1)]} continue
		set done1($p1) $p2
	}
	set result {}
	set p1s [array names done1]
	foreach p1 $p1s {
		set line1 [lindex $list1 $p1]
		set line2 [lindex $list2 $done1($p1)]
		lappend result [list [join [lrange $line1 1 6] \t] [join [lrange $line2 7 end] \t]]
	}
	foreach line1 [list_sub $list1 -exclude $p1s] {
		lappend result [list [join [lrange $line1 1 6] \t] $emptyout]
	}
	return $result
}

proc annotatesv {file dbfile name annotfile dbinfo {margin 30} {lmargin 300} {tmargin 300} {overlap 75}} {
# putsvars file dbfile name annotfile dbinfo
	set maxmargin [max $margin $lmargin $tmargin]
	#
	# open varfile
	catch {svclose $f1}
	set f1 [svopen $file 1 $name -]
	foreach {fh1 side1 poss1 makeref1 makealt1 header1} $f1 break
	if {[lrange $header1 0 5] ne {chromosome begin end type ref alt}} {
		error "error in compar file $file, must have (or can deduce) fields: chromosome begin end type ref alt"
	}
	set line1 [svgetline $f1]
	#
	# open dbfile
	# set dataposs [dict get $dbinfo dataposs]
	set outfields [dict get $dbinfo outfields]
	if {![llength $outfields]} {
		error "no outfields found (or specified) for $dbfile: cannot annotate sv"
	}
	catch {svclose $f2}
	set f2 [svopen $dbfile 2 $name $outfields]
	foreach {fh2 side2 poss2 makeref2 makealt2 header2} $f2 break
	if {[lrange $header2 0 5] ne {chromosome begin end type ref alt}} {
		error "error in compar file $dbfile, must have (or can deduce) fields: chromosome begin end type ref alt"
	}
	set line2 [svgetline $f2]
	#
	# open (output) annotfile	
	set empty [list_fill [llength $outfields] {}]
	set emptyout [join $empty \t]
	set newh [dict get $dbinfo newh]
	set o [open $annotfile.temp w]
	puts $o [join $newh \t]
	if {[gziscompressed $file]} {
		set file "|[gzcat $file] '$file'"
	}
	#
	# go over files
	set ::cchr {}
	set ::cpos 0
	while {[llength $line1] || [llength $line2]} {
		# get overlapping lines from both files in the following format (locfields):
		# {src chr1 begin end type start1 end1 size zyg chr2 start2 end2 (complete line)}
		set list [svmulticompar_getlist $f1 line1 $f2 line2 $maxmargin]
		# putsvars list line1 line2
		# puts [join [list_subindex $list {0 1 2 3 4 5}] \n]
		# puts [join $list \n]\n\n
		set sides [list_remdup [list_subindex $list 0]]
		# next if only annotation (2) in list
		# no output if only annotation
		if {$sides eq "2"} continue
		# output empty for each line if only src (1) in list
		if {$sides eq "1"} {
			# from single source
			foreach line $list {
				puts $o $emptyout
			}
			continue
		}
		# join $list \n
		#
		# multi source, match first
		# split on type
		unset -nocomplain todo
		unset -nocomplain typesa
		foreach line $list {
			set side [lindex $line 0]
			set type [lindex $line 4]
			if {$type eq "trans"} {set type bnd}
			set typesa($type) 1
			lappend todo($type,$side) $line
		}
		set list {}
		set types [lsort [array names typesa]]
		foreach type $types {
			if {![info exists todo($type,1)]} {
				# no output if only annotation
			} elseif {![info exists todo($type,2)]} {
				foreach l1 $todo($type,1) {
					lappend list [list [join [lrange $l1 1 6] \t] $emptyout]
				}
			} else {
				# foreach {list1 list2} [list $todo($type,1) $todo($type,2)] break
				set list1 $todo($type,1)
				set list2 $todo($type,2)
				set plist [annotate_sv_groupdists $list1 $list2 $emptyout $margin $lmargin $tmargin $overlap]
				foreach el $plist {
					lappend list $el
				}
			}
		}
		if {[llength $list] > 1} {
			set list [bsort $list]
		}
		foreach line $list {
			puts $o [lindex $line 1]
		}
	}
	flush $o
	close $o
	svclose $f1
	svclose $f2
	file rename -force -- $annotfile.temp $annotfile
}
