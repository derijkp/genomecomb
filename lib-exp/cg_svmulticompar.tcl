#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

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

proc svmulticompar_groupdists {dlist type} {
	set sizepos 7
	set curdist [lindex $dlist 0 $sizepos]
	set result {}
	set newlist {}
	set mindiff 35
	if {$type eq "inv"} {
		set mindiff 400
	}
	foreach line $dlist {
		set d [lindex $line $sizepos]
		set diff [expr {abs($d-$curdist)}]
		set sizediff [min [max [expr {round(0.25*$curdist)}] $mindiff] 8000]
		if {$diff > $sizediff} {
			lappend result $newlist
			set curdist [lindex $line $sizepos]
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
	set curstart1 [lindex $plist 0 2]
	set curend1 [lindex $plist 0 3]
	set result {}
	set newlist {}
	foreach line $plist {
		set start1 [lindex $line 2]
		set end1 [lindex $line 3]
		if {[overlap $start1 $end1 $curstart1 $curend1] > 1} {
			lappend newlist $line
		} else {
			lappend result $newlist
			set newlist {}
			set curstart1 $start1
			set curend1 $end1
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
	foreach {src chr1 begin end type start1 end1} $line1 break
	foreach {src chr2 begin end type start2 end2} $line2 break
	if {$chr1 ne $chr2} {
		set nchr1 [chr2num $chr1]
		set nchr2 [chr2num $chr2]
		return [expr {$nchr2 - $nchr1}]
	}
	if {[expr {$end2+20}] < $start1} {
		return -1
	} elseif {[expr {$end1+20}] < $start2} {
		return 1
	} else {
		return 0
	}
}

proc svmulticompar_getline {f poss {type 1}} {
	global cchr cpos
	while 1 {
		set line [split [gets $f] \t]
		if {[llength $line]} break
		if {[eof $f]} {return {}}
	}
	set cur [list_sub $line $poss]
	if {[lindex $cur 3] eq "trans"} {
		set endpos [expr {[lindex $cur 5]+200}]
		lset cur 6 0
		lset cur 9 $endpos
		lset cur 10 [expr {$endpos+400}]
	}
	if {[lindex $cur 0] ne $cchr} {
		putslog "Starting chromosome [lindex $cur 0]"
		set cchr [lindex $cur 0]
		set cpos 1000000
	}
	if {[lindex $cur 2] > $cpos} {
		putslog $cpos
		incr cpos 1000000
	}
	list_concat $type $cur $line
}

proc svmulticompar_write {o id group poss2 dummy1 dummy2 {ddummy1 {}} {ddummy2 {}}} {
	unset -nocomplain todo
	set todo(1) {}; set todo(2) {}
	foreach line $group {
		lappend todo([lindex $line 0]) $line
	}
	if {![llength $todo(1)] || ![llength $todo(2)]} {
		set udummy1 $dummy1; set udummy2 $dummy2
	} else {
		set udummy1 $ddummy1; set udummy2 $ddummy2
	}
	set max [max [llength $todo(1)] [llength $todo(2)]]
	if {$max > 1} {append id -$max}
	# count keeps track of how many matches we have for this sv
	if {[llength $todo(1)]} {
		set count [lindex $todo(1) 0 $::countpos]
	} else {
		set count 0
	}
	if {[llength $todo(2)]} {
		incr count
	}
	foreach l1 $todo(1) l2 $todo(2) {
		if {[llength $l2]} {
			set oline2 [lrange $l2 12 end]
			set cur2 [list_sub $oline2 $poss2]
			set merge [list_sub $oline2 $::mergeget2]
		} else {
			set oline2 $udummy2
		}
		if {[llength $l1]} {
			set oline1 [lrange $l1 12 end]
			set merge [list_sub $oline1 $::mergeget1]
		} else {
			set oline1 $udummy1
			set oline1 [lreplace $oline1 2 12 {*}$cur2]
		}
		lset oline1 0 $id
		lset oline1 1 $count
		set oline1 [list_sub $oline1 -exclude $::mergeposs1]
		set oline2 [list_sub $oline2 -exclude $::mergeposs2]
		puts $o [join [list_concat $oline1 $oline2 $merge] \t]
	}
}

proc svmulticompar_getlist {f1 poss1 len1 line1Var f2 poss2 len2 line2Var} {
	# lines have the following format:
	# {src chr1 begin end type start1 end1 size zyg chr2 start2 end2 (complete line)}
# puts [join [list_subindex [list $line1 $line2] {0 1 2 3 4 5}] \n]
	upvar $line1Var line1
	upvar $line2Var line2
	set startpos 5
	set endpos 6
	set chrpos 1
	set list {}
	set compar [svmulticompar_compar $line1 $line2]
	if {$compar >= 0} {
		set listchr [lindex $line1 1]	
		set liststart [lindex $line1 $startpos]
		set listend [expr {[lindex $line1 $endpos]+20}]
	} elseif {$compar < 0} {
		set listchr [lindex $line2 1]	
		set liststart [lindex $line2 $startpos]
		set listend [expr {[lindex $line2 $endpos]+20}]
	}
	while {![eof $f1] || ![eof $f2]} {
		set match 0
#putsvars listchr liststart listend line1 line2
		if {[llength $line1]} {
			set listchr1 [lindex $line1 $chrpos]
			set liststart1 [lindex $line1 $startpos]
			if {($listchr1 == $listchr) && ($liststart1 < $listend)} {
				lappend list $line1
				set temp [expr {[lindex $line1 $endpos]+20}]
				if {$temp > $listend} {
					set listend $temp
				}
				set line1 [svmulticompar_getline $f1 $poss1 1]
				set match 1
			}
		}
		if {[llength $line2]} {
			set listchr2 [lindex $line2 $chrpos]
			set liststart2 [lindex $line2 $startpos]
			if {($listchr2 == $listchr) && ($liststart2 < $listend)} {
				lappend list $line2
				set temp [expr {[lindex $line2 $endpos]+20}]
				if {$temp > $listend} {
					set listend $temp
				}
				set line2 [svmulticompar_getline $f2 $poss2 2]
				set match 2
			}
		}
		if {!$match} break
	}
	return $list
}

proc svmulticompar {svfile1 svfile2} {
	# set locfields {chr1 start1 end1 type size zyg chr2 start2 end2}
	set mergefields {LeftRepeatClassification RightRepeatClassification LeftGenes RightGenes XRef DeletedTransposableElement KnownUnderrepresentedRepeat FrequencyInBaselineGenomeSet}
	lappend mergefields overlappingGene	knownCNV
	set locfields {chromosome begin end type start1 end1 size zyg chr2 start2 end2}
	set ::countpos [expr {[llength $locfields]+2}]
	if {![file exists $svfile1]} {
		set o [open $svfile1 w]
		puts $o id\tcount\t[join $locfields \t]
		close $o
	}
	catch {close $f1}; catch {close $f2}; catch {close $o}; catch {file delete $tempfile2}
	#file copy -force comparsvcg.tsv.old comparsvcg.tsv
	set o [open $svfile1.temp w]
	#
	# open compar file
	set f1 [gzopen $svfile1]
	set header1 [tsv_open $f1]
	set len1 [llength $header1]
	set poss1 [list_cor $header1 $locfields]
	set mergefields1 [list_common $header1 $mergefields]
	set ::mergeposs1 [list_cor $header1 $mergefields1]
	set dummy1 [list_fill [llength $header1] {}]
	set ddummy1 [list_fill [llength $header1] d]
	#
	# open add file
	set name [file root [file tail $svfile2]]
	set name [lindex [split $name -] end]
	# set tempfile2 [tempfile]
	# cg select -s "chr1 start1" < $svfile2 > $tempfile2
	# set f2 [open $tempfile2]
	set f2 [gzopen $svfile2]
	set header2 [tsv_open $f2]
	set len2 [llength $header2]
	set poss2 [list_cor $header2 $locfields]
	set mergefields2 [list_common $header2 $mergefields]
	set ::mergeposs2 [list_cor $header2 $mergefields2]
	set dummy2 [list_fill [llength $header2] {}]
	set ddummy2 [list_fill [llength $header2] d]
	#
	set finalmerge [list_union $mergefields1 $mergefields2]
	set ::mergeget1 [list_cor $header1 $finalmerge]
	set ::mergeget2 [list_cor $header2 $finalmerge]
	#
	# make new header
	set header [list_sub $header1 -exclude $::mergeposs1]
	foreach field [list_sub $header2 -exclude $::mergeposs2] {append header \t${field}-$name}
	if {[llength $finalmerge]} {append header \t[join $finalmerge \t]}
	puts $o [join $header \t]
	#
	# go over files
	set ::cchr {}
	set ::cpos 0
	set line1 [svmulticompar_getline $f1 $poss1 1]
	set line2 [svmulticompar_getline $f2 $poss2 2]
	set did 1
	set sizepos 7; set typepos 4
	while {![eof $f1] || ![eof $f2]} {
		# get overlapping lines from both files in the following format (locfields):
		# {src chr1 begin end type start1 end1 size zyg chr2 start2 end2 (complete line)}
		set list [svmulticompar_getlist $f1 $poss1 $len1 line1 $f2 $poss2 $len2 line2]
#puts ----
#putsvars list
#puts [join [list_subindex $list {0 1 2 3 4 5}] \n]
#if {[lsearch [list_subindex $list 3] 67759423] != -1} {error STOPPED}
#if {[lindex $list 0 4] eq "trans"} {error STOP}
#puts [join $list \n]\n\n
# join $list \n\n
		if {[llength $list] == 1} {
			svmulticompar_write $o $did $list $poss2 $dummy1 $dummy2 $ddummy1 $ddummy2
			incr did
			continue
		}
		# split on type
		unset -nocomplain todo
		foreach line [lsort -integer -index $sizepos $list] {
			set type [lindex $line $typepos]
			set size [lindex $line $sizepos]
			lappend todo($type) $line
		}
		set list {}
		foreach type [array names todo] {
			if {[llength $todo($type)] < 2} {
				lappend list $todo($type)
				continue
			}
			set plists [svmulticompar_groupdists $todo($type) $type]
			lappend list {*}$plists
#			foreach plist $plists {
#				set plist [lsort -integer -index 3 $plist]
#				set clists [svmulticompar_groupcompatible $plist]
#				foreach clist $clists {
#					lappend list $clist
#				}
#			}
		}
		foreach group $list {
#puts [join [list_subindex $group {0 1 2 3 4 5}] \n]
			svmulticompar_write $o $did $group $poss2 $dummy1 $dummy2 $ddummy1 $ddummy2
			incr did
		}
	}
	flush $o
	close $o
	close $f1
	close $f2
	# file delete $tempfile2
	# cg select -s {chr1 start1} $svfile1.temp $svfile1.temp2
	# file delete $svfile1.temp
	file rename -force $svfile1 $svfile1.old
	file rename $svfile1.temp $svfile1
	putslog "finished adding $name to $svfile1"

}

proc cg_svmulticompar {args} {
	if {[llength $args] < 2} {
		puts "Wrong number of arguments"
		errorformat svmulticompar
		exit 1
	}
	set done {}
	foreach {compar_file} $args break
	if {[file exists $compar_file]} {
		set list [cg select -h $compar_file]
		set poss [list_find -glob $list start1-*]
		set done [list_sub $list $poss]
		set done [list_regsub -all {^start1-} $done {}]
	}
	set files [lrange $args 1 end]
	foreach file $files {
		set name [lindex [split [file tail [file root $file]] -] end]
		if {[inlist $done $name]} {
			putslog "Skipping $file: $name already present"
		} else {
			putslog "Adding $file"
			svmulticompar $compar_file $file
		}
	}
}

if {[info exists argv0] && [file tail [info script]] eq [file tail $argv0]} {
	package require pkgtools
	set appdir [file dir [pkgtools::startdir]]
	lappend auto_path $appdir/lib
	append env(PATH) :[file dir [file dir $appdir]]/bin:$appdir/bin
	package require Extral
	set ::base [file tail [info script]]
	cg_svmulticompar {*}$argv
}


if 0 {
	set svfile1 temp
	set svfile2 cmt71/cmt71_02_a/cgsv-cmt71_02_a.tsv
	set svfile2 cmt71/cmt71_07_b/cgsv-cmt71_07_b.tsv
	rm temp
	cg svmulticompar temp cmt71/cmt71_02_a/cgsv-cmt71_02_a.tsv
	cg svmulticompar temp cmt71/cmt71_07_b/cgsv-cmt71_07_b.tsv

	cd /complgen/projects/ep861
	cg svmulticompar temp ep861.03/cgsv-ep861.03.tsv ep861.04/cgsv-ep861.04.tsv


	# convert old
	set base GS103
	cd /complgen/sv/$base
	set file $base-sv.tsv
	catch {close $o}; catch {close $f}
	set o [open $file w]
	set h {check chr patchstart pos type size zyg problems gapsize/chr2 quality numreads numnontrf weight patchsize slope1 sd1 slope2 sd2 totnum psdiff threads exnum}
	set h {check chr1 start1 end1 type size zyg problems chr2 start2 end2 quality numreads numnontrf weight patchsize slope1 sd1 slope2 sd2 totnum psdiff threads exnum}
	puts $o [join $h \t]
	foreach chr {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X} {
		set f [open $base-$chr-paired-sv.tsv]
		gets $f
		while {![eof $f]} {
			set line [split [gets $f ] \t]
			if {![llength $line]} continue
			set type [lindex $line 4]
			if {$type eq "trans"} {
				set chr2 [lindex $line 8]
				set start2 0
				set end2 0
			} elseif {$type eq "inv"} {
				set chr2 [lindex $line 1]
				set start2 [expr {[lindex $line 3]+abs([lindex $line 5])}]
				set end2 [expr {$start2+abs([lindex $line 13])}]
			} else {
				set chr2 [lindex $line 1]
				set start2 [expr {[lindex $line 3]+abs([lindex $line 8])}]
				set end2 [expr {$start2+abs([lindex $line 13])}]
			}
			lset line 1 chr$chr
			set line [lreplace $line 8 8 chr$chr2 $start2 $end2]
			puts $o [join $line \t]
		}
		close $f
	}
	close $o
	cg select -q {$type != "ins" && !($type == "del" && $size < 250)} $file msv-a$base.tsv
	file rename msv-a$base.tsv /complgen/1.8/svcg/msv-a$base.tsv

}
