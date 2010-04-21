# < 0 if comp1 < comp2
# > 0 if comp1 > comp2
proc comparepos {comp1 comp2} {
	if {![llength $comp1]} {return 1}
	if {![llength $comp2]} {return -1}
	foreach {chr1 pos1} $comp1 break
	foreach {chr2 pos2} $comp2 break
	if {$chr1 eq $chr2} {
		return [expr {$pos1-$pos2}]
	} else {
		set chr1 [chr2num $chr1]
		set chr2 [chr2num $chr2]
		return [expr {$chr1-$chr2}]
	}
}

proc sequenced {r1 comp2} {
	if {$r1 eq 1} {return 1}
	global cache
	if {[llength $comp2] != 2} {return 0}
	foreach {chr pos} $comp2 break
	set chr [chr2num $chr]
	set line $cache($r1)
	while {[llength $line]} {
		foreach {rchr rstart rend} $line break
		set rchr [chr2num $rchr]
		if {$rchr > $chr} break
		if {$rchr == $chr} {
			if {$rstart > $pos} break
			if {($rchr == $chr) && ($pos < $rend) && ($pos >= $rstart)} {
				set cache($r1) $line
				return 1
			}
		}
		set line [gets $r1]
	}
	set cache($r1) $line
	return 0
}

array set compare_annot_join_trans {
	COMPATIBLE 1 MISSENSE 2 DELETE 3 INSERT 4 DELETE+ 5 INSERT+ 6 NONSTOP 7 NONSENSE 8 FRAMESHIFT 9
}

proc compare_annot_join {compar sample cur1 cur2} {
	global comparposs1 mergeposs1 comparposs2 mergeposs2 dummy1 dummy2 restposs1 restposs2
	global compare_annot_join_trans
	if {[inlist {{} -} $cur1]} {
		set region [list_sub $cur2 $comparposs2]
		set merge [list_sub $cur2 $mergeposs2]
		if {$cur1 eq "-"} {
			set cur1 $dummy1
		} else {
			set cur1 [list_change $dummy1 {- {}}]
		}
	} elseif {[inlist {{} -} $cur2]} {
		set region [list_sub $cur1 $comparposs1]
		set merge [list_sub $cur1 $mergeposs1]
		if {$cur2 eq "-"} {
			set cur2 $dummy2
		} else {
			set cur2 [list_change $dummy2 {- {}}]
		}
	} else {
		set region [list_sub $cur1 $comparposs1]
		set merge {}
		foreach el1 [list_sub $cur1 $mergeposs1] el2 [list_sub $cur2 $mergeposs2] {
			lappend merge [list_union $el1 $el2]
		}
	}
	set result [list $compar $sample]
	lappend result {*}$region
	lappend result {*}[list_sub $cur1 $restposs1]
	lappend result {*}[list_sub $cur2 $restposs2]
	lappend result {*}$merge
	return [join $result \t]
}

proc compare_annot_getline {f} {
	set cur1 {}
	while {![eof $f]} {
		set cur1 [split [gets $f] \t]
		if {[llength $cur1]} break
	}
	return $cur1
}

proc compare_annot {id1 file1 regfile1 id2 file2 regfile2 outfile} {
	global cache comparposs1 mergeposs1 comparposs2 mergeposs2 dummy1 dummy2 restposs1 restposs2

	catch {close $f1}; catch {close $f2}; catch {close $r1}; catch {close $r2}; catch {close $o}
	set comparfields {chromosome begin end type}
	set mergefields {xRef trf str segdup selfchain repeat rna checked geneId mrnaAcc proteinAcc orientation exonCategory exon codingRegionKnown aaCategory nucleotidePos proteinPos aaAnnot aaCall aaRef effect neffect}
	set allelefields {alleleSeq1 alleleSeq2}
	set annotvarfields {locus chromosome begin end type alleleSeq1 alleleSeq2 totalScore1 totalScore2 xRef geneId mrnaAcc proteinAcc orientation exonCategory exon codingRegionKnown aaCategory nucleotidePos proteinPos aaAnnot aaCall aaRef}
	set f1 [open $file1]
	set header1 [split [gets $f1] \t]
	set comparposs1 [list_cor $header1 $comparfields]
	set alleleposs1 [list_cor $header1 $allelefields]
	set ralleleposs1 [list_reverse $alleleposs1]
	if {([lsearch $comparposs1 -1] != -1) || ([lsearch $alleleposs1 -1] != -1)} {
		puts stderr "header error in fannot_varfile1 $file1"
		exit 1
	}
	set mergeposs1 [list_remove [list_cor $header1 $mergefields] -1]
	set dummy1 [list_fill [llength $header1] -]
	set f2 [open $file2]
	set header2 [split [gets $f2] \t]
	set comparposs2 [list_cor $header2 $comparfields]
	set alleleposs2 [list_cor $header2 $allelefields]
	set ralleleposs2 [list_reverse $alleleposs2]
	if {([lsearch $comparposs2 -1] != -1) || ([lsearch $alleleposs2 -1] != -1)} {
		puts stderr "header error in fannot_varfile2 $file1"
		exit 1
	}
	set mergeposs2 [list_remove [list_cor $header2 $mergefields] -1]
	set dummy2 [list_fill [llength $header2] -]
	if {$regfile1 ne ""} {
		set r1 [opencgifile $regfile1 header]
		if {[lrange $header 0 2] ne "chromosome begin end"} {
			puts stderr "header error in region_file1 $file1"
			exit 1
		}
		set cache($r1) [split [gets $r1] \t]
	} else {
		puts "no region file for file1"
		set r1 1
	}
	if {$regfile2 ne ""} {
		set r2 [opencgifile $regfile2 header]
		if {[lrange $header 0 2] ne "chromosome begin end"} {
			puts stderr "header error in region_file2 $file2"
			exit 1
		}
		set cache($r2) [split [gets $r2] \t]
	} else {
		puts "no region file for file2"
		set r2 1
	}
	# start
	set o [open $outfile w]
	# make output header
	set oheader [list_concat {compar sample} $comparfields]
	set restfields1 [list_lremove [list_lremove $header1 $mergefields] $comparfields]
	set restposs1 [list_cor $header1 $restfields1]
	foreach field $restfields1 {
		lappend oheader ${field}-1
	}
	set restfields2 [list_lremove [list_lremove $header2 $mergefields] $comparfields]
	set restposs2 [list_cor $header2 $restfields2]
	foreach field $restfields2 {
		lappend oheader ${field}-2
	}
	set oheader [list_concat $oheader $mergefields]
	puts $o [join $oheader \t]
	set cur1 [split [gets $f1] \t]
	set comp1 [list_sub $cur1 $comparposs1]
	set cur2 [split [gets $f2] \t]
	set comp2 [list_sub $cur2 $comparposs2]
	set num 1
	while {![eof $f1] || ![eof $f2]} {
		incr num
		if {![expr {$num % 100000}]} {putslog $num}
		set d [comparepos $comp1 $comp2]
		if {$d == 0} {
			set s1 [sequenced $r1 $comp1]
			set s2 [sequenced $r2 $comp2]
			if {!$s1 || !$s2} {
				puts $o [compare_annot_join fl $id1,$id2 $cur1 $cur2]
			} else {
				#set type [list_remove [split [lindex $cur1 4] _] ref-consistent ref-inconsistent =]
				#lset cur1 4 $type
				#set type [list_remove [split [lindex $cur2 4] _] ref-consistent ref-inconsistent =]
				#lset cur2 4 $type
				if {($comp1 eq $comp2) && (
					([list_sub $cur1 $alleleposs1] eq [list_sub $cur2 $alleleposs2])
					|| ([list_sub $cur1 $ralleleposs1] eq [list_sub $cur2 $alleleposs2])
					)} {
					puts $o [compare_annot_join sm $id1,$id2 $cur1 $cur2]
				} else {
					puts $o [compare_annot_join mm $id1,$id2 $cur1 $cur2]
				}
			}
			set cur1 [compare_annot_getline $f1]
			set comp1 [list_sub $cur1 $comparposs1]
			set cur2 [compare_annot_getline $f2]
			set comp2 [list_sub $cur2 $comparposs2]
		} elseif {$d < 0} {
			while {[comparepos $comp1 $comp2] < 0} {
				set s [sequenced $r1 $comp1]
				if {!$s} {
					puts $o [compare_annot_join fl $id1 $cur1 -]
				} else {
					set s [sequenced $r2 $comp1]
					if {$s} {
						puts $o [compare_annot_join df $id1 $cur1 {}]
					} else {
						puts $o [compare_annot_join un $id1 $cur1 -]
					}
				}
				if {[eof $f1]} break
				set cur1 [compare_annot_getline $f1]
				set comp1 [list_sub $cur1 $comparposs1]
				if {[llength $cur1]} break
			}
		} else {
			while {[comparepos $comp1 $comp2] > 0} {
				set s [sequenced $r2 $comp2]
				if {!$s} {
					puts $o [compare_annot_join fl $id2 - $cur2]
				} else {
					set s [sequenced $r1 $comp2]
					if {$s} {
						puts $o [compare_annot_join df $id2 {} $cur2]
					} else {
						puts $o [compare_annot_join un $id2 - $cur2]
					}
				}
				if {[eof $f2]} break
				set cur2 [compare_annot_getline $f2]
				set comp2 [list_sub $cur2 $comparposs2]
				if {[llength $cur2]} break
			}
		}
	}

	close $f1; close $f2; close $r1; close $r2; close $o
}

proc lset_always {varName pos value} {
	upvar $varName var
	if {[catch {lset var $pos $value}]} {
		set var [list_concat $var [list_fill [expr {$pos - [llength $var]}] {}]]
		lappend var $value
	}
}

array set chrtrans {M 95 X 96 Y 97}
proc chr2num {chr} {
	if {[isint $chr]} {return $chr}
	regsub ^chr $chr {} chr
	set nchr [get ::chrtrans($chr) $chr]
	return $nchr
}

proc annot_compare_region {compar_file reg_file field tvalue fvalue} {
	set f1 [open $compar_file]
	set header1 [gets $f1]
	set poss1 [list_cor $header1 {chromosome begin end}]
	if {[lsearch $poss1 -1] != -1} {
		set poss1 [list_cor $header1 {chr patchstart pos}]
	}
	if {[lsearch $poss1 -1] != -1} {
		error "Wrong header"
	}
	set f2 [open $reg_file]
	set poss2 [open_region $f2]
	set num 0
	set line1 [split [gets $f1] \t]
	foreach {chr1 start1 end1} [list_sub $line1 $poss1] break
	set nchr1 [chr2num $chr1]
	set line2 [split [gets $f2] \t]
	foreach {chr2 start2 end2} [list_sub $line2 $poss2] break
	set nchr2 [chr2num $chr2]
	set ipos [lsearch $header1 $field]
	if {$ipos == -1} {
		lappend header1 $field
	}
	set o stdout
	puts $o [join $header1 \t]
	while 1 {
		incr num
		if {![expr $num%100000]} {putslog $num}
		set cur $fvalue
		while 1 {
			# putsvars chr1 chr2 start1 end1 start2 end2
			if {$nchr2 > $nchr1} break
			if {$nchr2 == $nchr1} {
				if {$start2 > $end1} break
				if {$start2 == $end1} {
					if {($start2 == $end2) || ($start1 == $end1)} {set cur $tvalue}
					break
				}
				if {$start2 > $start1} {
					set cur $tvalue
					break
				} elseif {$end2 >= $start1} {
					set cur $tvalue
					break
				}
			}
			set line2 [getnotempty $f2]
			if {![llength $line2]} break
			foreach {chr2 start2 end2} [list_sub $line2 $poss2] break
			set nchr2 [chr2num $chr2]
		}
		if {$ipos == -1} {
			lappend line1 $cur
		} elseif {$cur eq $tvalue} {
			lset_always line1 $ipos $cur
		}
		# puts $o #x\tx\tx\t[join $line2 \t]
		puts $o [join $line1 \t]
		if {[eof $f1]} break
		set line1 [getnotempty $f1]
		if {![llength $line1]} break
		foreach {chr1 start1 end1} [list_sub $line1 $poss1] break
		set nchr1 [chr2num $chr1]
	}
	close $f1
	close $f2
}

proc reannot_compare {compar_file dir1 dir2 outfile} {
	puts "reannotating $compar_file -> $outfile"
	set name1 [file tail $dir1]
	set name2 [file tail $dir2]
	set f [open $compar_file]
	set header [split [gets $f] \t]
	set todo {}
	lappend todo [list [lsearch $header refcons-1] rc $dir1/reg_refcons-$name1.tsv]
	lappend todo [list [lsearch $header refcons-2] rc $dir2/reg_refcons-$name2.tsv]
	lappend todo [list [lsearch $header cluster-1] cl $dir1/reg_cluster-$name1.tsv]
	lappend todo [list [lsearch $header cluster-2] cl $dir2/reg_cluster-$name2.tsv]
	list_foreach {field value regfile} $todo {
		annot_region_init $regfile
	}
	annot_coverage_init $dir1
	annot_coverage_init $dir2
	set o [open $outfile w]
	set poss [list_cor $header {chromosome begin end refscore-1 coverage-1 refscore-2 coverage-2}]
	foreach {var p}  {r1pos 3 c1pos 4 r2pos 5 c2pos 6} {
		set $var [lindex $poss $p]
	}
	puts $o [join $header \t]
	set num 0
	while {![eof $f]} {
		incr num
		if {![expr $num%10000]} {putslog $num}
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		foreach {chr begin end r1 c1 r2 c2} [list_sub $line $poss] break
		if {$r1 eq ""} {
			foreach {r c} [annot_coverage_get $dir1 $chr $begin] break
			lset line $r1pos $r
			lset line $c1pos $c
		}
		if {$r2 eq ""} {
			foreach {r c} [annot_coverage_get $dir2 $chr $begin] break
			lset line $r2pos $r
			lset line $c2pos $c
		}
		list_foreach {field value regfile} $todo {
			if {[lindex $line $field] == "-"} continue
			# if {[lindex $line $field] != ""} continue
			set r [annot_region_get $regfile $chr $begin $end]
			if {$r} {lset line $field $value} else {lset line $field {}}
		}
		puts $o [join $line \t]
	}
	close $o
	close $f
	list_foreach {field value regfile} $todo {
		annot_region_close $regfile
	}
	annot_coverage_close $dir1
	annot_coverage_close $dir2
}

if 0 {

lappend auto_path ~/dev/completegenomics/lib
package require Extral
package require Tclx
signal -restart error SIGINT

set compar_file /complgen/testcompar_GS102_GS103/compar_GS102_GS103.tsv
cd [file dir $compar_file]
set dir1 /complgen/GS102
set dir2 /complgen/GS103
reannot_compare $compar_file $dir1 $dir2 atemp

set compar_file compar/78vs79_compar.tsv
set reg_file GS00102/reg-refcons-GS000000078-ASM.tsv
set field refcons
set tvalue rc
set fvalue ""

set compar_file compar/78vs79_compar.tsv
set reg_file GS00103/reg-refcons-GS000000079-ASM.tsv
set field refcons
set tvalue rc
set fvalue ""
set poss1 {3 4 5}
set poss2 {0 1 2}

set compar_file compar/78vs79_compar-filter-rc2.tsv
set reg_file /data/db/_data_db_ucsc-simple_repeats.tsv
set field trf
set tvalue trf
set fvalue ""


set o [open compar/test w]

cd /complgen/
cd /media/passport/complgen

set file1 GS00102/annotvar-GS000000078-ASM.tsv
set file2 GS00103/annotvar-GS000000079-ASM.tsv
set id1 78
set id2 79
set regfile1 GS00102/reg-GS000000078-ASM.tsv
set regfile2 GS00103/reg-GS000000079-ASM.tsv
set outfile compar/78vs79_compar.tsv

catch {close $f1}
catch {close $f2}
catch {close $r1}
catch {close $r2}
catch {close $o}
unset -nocomplain cache

set f [open compar/78vs79_compar.tsv]
set line [split [gets $f] \t]
llength $line

}
