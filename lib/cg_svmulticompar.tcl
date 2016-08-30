#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

if 0 {

package require Tclx
signal -restart error SIGINT
lappend auto_path ~/dev/genomecomb/lib ~/dev/genomecomb/lib-exp
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
	set sizepos $::locpos(size)
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

proc svmulticompar_compar {line1 line2 {margin 20}} {
	if {![llength $line1]} {return 1}
	if {![llength $line2]} {return -1}
	# must be adapted on locfields changes
	foreach {src id1 chr1 begin end type start1 end1} $line1 break
	foreach {src id2 chr2 begin end type start2 end2} $line2 break
	if {$id1 ne "" && $id1 eq $id2} {
		return 0
	}
	set chrcomp [chr_compare $chr1 $chr2]
	if {$chrcomp != 0} {
		return $chrcomp
	}
	if {[expr {$end2+$margin}] <= $start1} {
		return 1
	} elseif {[expr {$end1+$margin}] <= $start2} {
		return -1
	} else {
		return 0
	}
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

proc svmulticompar_write {id group poss2 dummy1 dummy2 {ddummy1 {}} {ddummy2 {}}} {
	global locpos
	set resultlist {}
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
			set oline2 [lrange $l2 $locpos(datastart) end]
			set cur2 [list_sub $oline2 $poss2]
			set merge [list_sub $oline2 $::mergeget2]
		} else {
			set oline2 $udummy2
		}
		if {[llength $l1]} {
			set oline1 [lrange $l1 $locpos(datastart) end]
			set merge [list_sub $oline1 $::mergeget1]
		} else {
			set oline1 $udummy1
			set oline1 [lreplace $oline1 1 $locpos(datastart) {*}$cur2]
		}
		lset oline1 0 $count
		lset oline1 1 $id
		if {[llength $::mergeposs1]} {
			set oline1 [list_sub $oline1 -exclude $::mergeposs1]
		}
		if {[llength $::mergeposs2]} {
			set oline2 [list_sub $oline2 -exclude $::mergeposs2]
		}
		lappend resultlist [list_concat $oline1 $oline2 $merge]
	}
	return $resultlist
}

proc svmulticompar_getlist {f1 poss1 len1 line1Var f2 poss2 len2 line2Var} {
	global locpos
	# lines have the following format:
	# {src id chr1 begin end type start1 end1 size zyg chr2 start2 end2 (complete line)}
# puts [join [list_subindex [list $line1 $line2] {0 1 2 3 4 5}] \n]
	upvar $line1Var line1
	upvar $line2Var line2
	set margin 20
	set idpos $locpos(id)
	set chrpos $locpos(chromosome)
	set startpos $locpos(start1)
	set endpos $locpos(end1)
	set list {}
	set compar [svmulticompar_compar $line1 $line2 $margin]
	if {$compar < 0} {
		set listchr [lindex $line1 $chrpos]	
		set liststart [lindex $line1 $startpos]
		set listend [expr {[lindex $line1 $endpos]+$margin}]
		set curid1 [lindex $line1 $idpos]
		set curid2 {}
	} elseif {$compar > 0} {
		set listchr [lindex $line2 $chrpos]	
		set liststart [lindex $line2 $startpos]
		set listend [expr {[lindex $line2 $endpos]+$margin}]
		set curid1 {}
		set curid2 [lindex $line2 $idpos]
	} else {
		set listchr [lindex $line1 $chrpos]	
		set liststart [min [lindex $line1 $startpos] [lindex $line2 $startpos]]
		set listend [expr {[max [lindex $line1 $endpos] [lindex $line2 $endpos]]+$margin}]
		set curid1 [lindex $line1 $idpos]
		set curid2 [lindex $line2 $idpos]
	}
	while {![eof $f1] || ![eof $f2]} {
		set match 0
		if {[llength $line1]} {
			set id1 [lindex $line1 $idpos]
			set listchr1 [lindex $line1 $chrpos]
			set liststart1 [lindex $line1 $startpos]
			if {($id1 ne "" && $id1 eq $curid1) || (($listchr1 == $listchr) && ($liststart1 < $listend))} {
				lappend list $line1
				set curid1 $id1
				set temp [expr {[lindex $line1 $endpos]+$margin}]
				if {$temp > $listend} {
					set listend $temp
				}
				set line1 [svmulticompar_getline $f1 $poss1 1]
				set match 1
			}
		}
		if {[llength $line2]} {
			set id2 [lindex $line2 $idpos]
			set listchr2 [lindex $line2 $chrpos]
			set liststart2 [lindex $line2 $startpos]
			if {($id2 ne "" && $id2 eq $curid2) || (($listchr2 == $listchr) && ($liststart2 < $listend))} {
				lappend list $line2
				set curid2 $id2
				set temp [expr {[lindex $line2 $endpos]+$margin}]
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
	global locpos

	# set locfields {chr1 start1 end1 type size zyg chr2 start2 end2}
	catch {close $f1}; catch {close $f2}; catch {close $o}; catch {file delete $tempfile2}
	set mergefields {LeftRepeatClassification RightRepeatClassification LeftGenes RightGenes XRef DeletedTransposableElement KnownUnderrepresentedRepeat FrequencyInBaselineGenomeSet}
	lappend mergefields overlappingGene	knownCNV
	# set locfields {chromosome begin end type start1 end1 size zyg chr2 start2 end2}
	set locfields {id chromosome begin end type start1 end1 size zyg chr2 start2 end2}
	set locpos(fields) [list_concat src $locfields]
	foreach f {id chromosome begin end type start1 end1 size start2 end2} {
		# add 1, because src will be prepended to lines
		set locpos($f) [expr {[lsearch $locfields $f]+1}]
	}
	set sizepos $locpos(size); set typepos $locpos(type)
	set locpos(datastart) [expr {[llength $locfields]+1}]
	set ::countpos [expr {[llength $locfields]+1}]
	if {![file exists $svfile1]} {
		set o [open $svfile1 w]
		puts $o count\t[join $locfields \t]
		close $o
	}
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
	set templen [expr {[llength $header1]+1}]
	set dummy1 [list_fill $templen {}]
	set ddummy1 [list_fill $templen d]
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
	# make poss2
	set poss2 [list_cor $header2 $locfields]
	set temp [tsv_basicfields $header2 4]
	set poss2 [lreplace $poss2 1 4 {*}$temp]
	set workfields {id chr1 begin end type start1 end1 size zyg chr2 start2 end2 (complete line)}
	# set locfields {id chromosome begin end type start1 end1 size zyg chr2 start2 end2}
	foreach {fld basicpos} {start1 1 end1 1 start2 2 end2 2 chr2 0} {
		set workpos [lsearch $workfields $fld]
		if {[lindex $poss2 $workpos] == -1} {lset poss2 $workpos [lindex $temp $basicpos]}
	}
	set mergefields2 [list_common $header2 $mergefields]
	set ::mergeposs2 [list_cor $header2 $mergefields2]
	set templen [llength $header2]
	set dummy2 [list_fill $templen {}]
	set ddummy2 [list_fill $templen d]
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
			set temp [lindex [svmulticompar_write $did $list $poss2 $dummy1 $dummy2 $ddummy1 $ddummy2] 0]
			puts $o [join $temp \t]
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
		foreach type [lsort [array names todo]] {
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
		set resultlist {}
		foreach group $list {
#puts [join [list_subindex $group {0 1 2 3 4 5}] \n]
			lappend resultlist {*}[svmulticompar_write $did $group $poss2 $dummy1 $dummy2 $ddummy1 $ddummy2]
			incr did
		}
		if {[llength $resultlist] > 1} {
			set resultlist [ssort -index $locpos(id) -natural $resultlist]
			set resultlist [ssort -index $locpos(end) -natural $resultlist]
			set resultlist [lsort -index $locpos(begin) -integer $resultlist]
		}
		foreach temp $resultlist {
			puts $o [join $temp \t]
		}
	}

	flush $o
	close $o
	gzclose $f1
	gzclose $f2
	# file delete $tempfile2
	# cg select -s {chr1 start1} $svfile1.temp $svfile1.temp2
	# file delete $svfile1.temp
	file rename -force $svfile1 $svfile1.old
	file rename -force $svfile1.temp $svfile1
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
	file rename -force msv-a$base.tsv /complgen/1.8/svcg/msv-a$base.tsv

}
