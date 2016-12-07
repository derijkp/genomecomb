#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc tsv_align_compareoverlap {comp1 comp2} {
	foreach {chr1 p11 p12} $comp1 break
	foreach {chr2 p21 p22} $comp2 break
	if {$chr1 ne $chr2} {
		if {$chr2 eq ""} {
			return -1
		} elseif {$chr1 eq ""} {
			return -2
		}
		if {[ssort -natural [list $chr1 $chr2]] eq [list $chr1 $chr2]} {return -1} else {return -2}
	}
	if {![isint $p11]} {
		return -1
	} elseif {![isint $p21]} {
		return -2
	} elseif {$p12 <= $p21} {
		return -1
	} elseif {$p22 <= $p11} {
		return -2
	} else {
		return 0
		return [overlap $p11 $p12 $p21 $p22]
	}
}

proc tsv_align {file1 file2 joinfields1 joinfields2 postfix1 postfix2} {
	catch {close $f1}; catch {close $f2}
	set o stdout
	# open files
	set f1 [open $file1]
	set header1 [split [gets $f1] \t]
	set comparposs1 [list_cor $header1 $joinfields1]
	if {[inlist $comparposs1 -1]} {
		error "$file1 does not contain all joinfields"
	}
	set f2 [open $file2]
	set header2 [split [gets $f2] \t]
	set comparposs2 [list_cor $header2 $joinfields2]
	if {[inlist $comparposs2 -1]} {
		error "$file2 does not contain all joinfields"
	}
	#
	set dummy1 [list_fill [llength $header1] {}]
	set dummy2 [list_fill [llength $header2] {}]
	# start making output files
	# make output header
	set oheader {}
	foreach field $header1 {
		lappend oheader ${field}$postfix1
	}
	# make output header
	foreach field $header2 {
		lappend oheader ${field}$postfix2
	}
	puts $o [join $oheader \t]
	set cur1 [split [gets $f1] \t]
	set comp1 [list_sub $cur1 $comparposs1]
	set cur2 [split [gets $f2] \t]
	set comp2 [list_sub $cur2 $comparposs2]
	set prevcomp1 {}
	set prevcomp2 {}
	# go over the files
	set num 1
	while {![eof $f1] || ![eof $f2]} {
		incr num
		if {![expr {$num % 100000}]} {putslog $num}
		set d [tsv_align_match $comp1 $comp2]
		if {$d == 0} {
			puts $o [join $cur1 \t]\t[join $cur2 \t]
			set cur1 [compare_annot_getline $f1]
			set comp1 [list_sub $cur1 $comparposs1]
			set cur2 [compare_annot_getline $f2]
			set comp2 [list_sub $cur2 $comparposs2]
		} elseif {$d == -1} {
			while {[tsv_align_match $comp1 $comp2] == -1} {
				puts $o [join $cur1 \t]\t[join $dummy2 \t]
				if {[eof $f1]} break
				set cur1 [compare_annot_getline $f1]
				set comp1 [list_sub $cur1 $comparposs1]
				if {![llength $cur1]} break
				incr num
				if {![expr {$num % 100000}]} {putslog $num}
			}
		} elseif {$d == -2} {
			while {[tsv_align_match $comp1 $comp2] == -2} {
				puts $o [join $dummy1 \t]\t[join $cur2 \t]
				if {[eof $f2]} break
				set cur2 [compare_annot_getline $f2]
				set comp2 [list_sub $cur2 $comparposs2]
				if {![llength $cur2]} break
				incr num
				if {![expr {$num % 100000}]} {putslog $num}
			}
		}
	}
	close $f1; close $f2
}

proc cg_tsv_align {args} {
	global scriptname action
	if {[llength $args] < 6} {
		error "format is: $scriptname $action file1 file2 joinfields1 joinfields2 postfix1 postfix2 ?method?"
	}
	foreach {file1 file2 joinfields1 joinfields2 postfix1 postfix2 method} $args break
	if {$method eq ""} {set method match}
	switch $method {
		match {tsv_align $file1 $file2 $joinfields1 $joinfields2 $postfix1 $postfix2}
	}
}

if 0 {

set params {compar-cnvcg.tsv ogtregions.tsv "chr-GS102 begin-GS102 end-GS102" "chromosome begin end" "" "-ogt" overlap}
set params {cnvcg-GS102.tsv cnvcg-GS103.tsv {chr begin end} {chr begin end} -GS102 -GS103 overlap}
set params {allvalsnps3.tsv ../multicompar/compar.tsv {chromosome begin end type} {chromosome begin end type} -val {}}
foreach {file1 file2 joinfields1 joinfields2 postfix1 postfix2 method} $params break

	tsv_align $file1 $file2 $joinfields1 $joinfields2 $postfix1 $postfix2 $method
	
paste acnvcg-GS102.tsv acnvcg-GS103.tsv acnvcg-GS102.tsv > cnvcg-compar.tsv

}

