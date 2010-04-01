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

proc svcompare {svfile1 svfile2} {

	catch {close $f1}
	catch {close $f2}
	catch {close $o}
	set tempfile [tempfile]
	set o [open $tempfile w]
	set name1 [lindex [split [file tail $svfile1] -] 0]
	set name2 [lindex [split [file tail $svfile2] -] 0]
	set tempfile1 [tempfile]
	cg select -s patchstart < $svfile1 > $tempfile1
	set f1 [open $tempfile1]
	set tempfile2 [tempfile]
	cg select -s patchstart < $svfile2 > $tempfile2
	set f2 [open $tempfile2]
	set header1 [split [gets $f1] \t]
	set len [llength $header1]
	set patchstartpos [lsearch $header1 patchstart]
	set end1pos [lsearch $header1 pos]
	set header2 [split [gets $f2] \t]
	set header [list_concat match sample $header1]
	foreach f $header2 {lappend header ${f}-2}
	puts $o [join $header \t]
	set poss1 [list_cor $header1 {patchstart pos type size zyg}]
	set poss2 [list_cor $header2 {patchstart pos type size zyg}]
	set list1 {}
	set list2 {}
	set line1 [split [gets $f1] \t]
	foreach {liststart1 listend1} [list_sub $line1 $poss1] break
	set line2 [split [gets $f2] \t]
	while {![eof $f1]} {
		while {![eof $f1]} {
			if {[llength $line1]} {
				foreach {start1 pos1 type1 size1 zyg1} [list_sub $line1 $poss1] break
				if {$start1 > $listend1} break
				lappend list1 $line1
				if {$start1 < $liststart1} {set liststart1 $start1}
				if {$pos1 > $listend1} {set listend1 $pos1}
			}
			set line1 [split [gets $f1] \t]
			set missing [max [expr {$len - [llength $line1]}] 0]
			if {$missing} {
				while {$missing} {lappend line1 {}; incr missing -1}
			}
		}
#if {[inlist [list_subindex $list1 $end1pos] 22603340]} {error STOPPED}
# join $list1 \n
# putsvars liststart1 listend1
		if {![llength $list1]} break
		foreach tline2 $list2 {
			foreach {start2 pos2 type2 size2 zyg2} [list_sub $tline2 $poss2] break
			if {$pos2 < $liststart1} {
				puts $o df\t$name2\t[join $tline2 \t]
				set list2 [list_remove $list2 $tline2]
				flush $o
			}
		}
		while {![eof $f2]} {
			if {[llength $line2]} {
				foreach {start2 pos2} [list_sub $line2 $poss2] break
				if {$start2 >= $listend1} break
				if {$pos2 >= $liststart1} {
					lappend list2 $line2
				} else {
					puts $o df\t$name2\t[join $line2 \t]
				}
			}
			set line2 [split [gets $f2] \t]
		}
# join $list2 \n
#if {[inlist [list_subindex $list2 $end1pos] 37469426]} {error STOPPED}
#if {([llength $list1] > 1) && ([llength $list2] > 1)} {error STOPPED}
		unset -nocomplain a
		if {[llength $list2]} {
			# check lists
			set matrix {}
			set pos1 -1
			# puts [join $list1 \n]\n\n[join $list2 \n]
			foreach l1 $list1 {
				incr pos1
				foreach {start1 end1 type1 size1 zyg1} [list_sub $l1 $poss1] break
				set pos2 -1
				foreach l2 $list2 {
					incr pos2
					foreach {start2 end2 type2 size2 zyg2} [list_sub $l2 $poss2] break
					if {$type2 ne $type1} continue
					set diff [expr {abs($size2-$size1)}]
					set sizediff [min [max [expr {round(0.1*$size1)}] 35] 400]
					set overlap [expr {[overlap $start1 $end1 $start2 $end2]/min(double($end1-$start1+1),double($end2-$start2+1))}]
					set dist [expr {abs($end2 - $end1)}]
					#puts "$pos1\t$pos2\t$type1\t$type2\t$diff\t$overlap\t$dist"
					if {($diff <= $sizediff) && (
						($overlap > 0.5) || (($overlap > 0) && ($dist < 60))
					)} {
						if {($zyg2 eq $zyg1)} {set s sm} else {set s mm}
						lappend matrix [list $pos1 $pos2 $diff $overlap $dist $s]
					}
				}
			}
			set matrix [lsort -integer -index 2 $matrix]
			set remain {}
			# first find best matches
			foreach line $matrix {
				foreach {pos1 pos2 diff overlap dist s} $line break
				if {[info exists a(1,$pos1)] || [info exists a(2,$pos2)]} {
					lappend remain $line
					continue
				}
				set l1 [lindex $list1 $pos1]
				set l2 [lindex $list2 $pos2]
				puts $o $s\t$name1,$name2\t[join $l1 \t]\t[join $l2 \t]
				set a(1,$pos1) 1
				set a(2,$pos2) 1
			}
			# see if there are any doubles left (2 matching 1 other)
			foreach line $remain {
				foreach {pos1 pos2 diff overlap dist s} $line break
				if {[info exists a(1,$pos1)] && [info exists a(2,$pos2)]} continue
				set l1 [lindex $list1 $pos1]
				set l2 [lindex $list2 $pos2]
				puts $o db\t$name1,$name2\t[join $l1 \t]\t[join $l2 \t]
				set a(1,$pos1) 1
				set a(2,$pos2) 1
			}
			set poss {}
			foreach n [array names a 2,*] {
				lappend poss [lindex [split $n ,] end]
			}
			set list2 [list_sub $list2 -exclude $poss]
		}
		# next block
		set pos -1
		foreach l1 $list1 {
			incr pos
			if {[info exists a(1,$pos)]} continue
			puts $o df\t$name1\t[join $l1 \t]
		}
		set list1 {}
		foreach {liststart1 listend1} [list_sub $line1 $poss1] break
	}


	foreach tline2 $list2 {
		puts $o df\t$name2\t[join $tline2 \t]
	}
	while {![eof $f2]} {
		if {![llength $line2]} {
			set line2 [split [gets $f2] \t]
			continue
		}
		puts $o df\t$name2\t[join $line2 \t]
		set line2 [split [gets $f2] \t]
	}
	flush $o
	close $f1
	close $f2
	file delete $tempfile1
	file delete $tempfile2
	close $o
	# cg select -s pos < $tempfile >@ stdout
	set chr [lindex [split [file tail $svfile1] -] 1]
	set outfile svcompar_[lindex [split [file tail $svfile1] -] 0]_[lindex [split [file tail $svfile2] -] 0]-$chr.tsv
	cg select -s pos < $tempfile > $outfile
	putslog "finished $outfile"

}
