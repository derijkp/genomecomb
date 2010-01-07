proc comparepos {cur1 cur2} {
	array set trans {X 100 Y 101 M 102}
	set chr1 [lindex $cur1 1]
	set chr2 [lindex $cur2 1]
	if {$chr1 eq $chr2} {
		set pos1 [lindex $cur1 2]
		set pos2 [lindex $cur2 2]
		return [expr {$pos1-$pos2}]
	} else {
		set chr1 [get trans($chr1) $chr1]
		set chr2 [get trans($chr2) $chr2]
		return [expr {$chr1-$chr2}]
	}
}

proc sequenced {r1 cur2} {
	if {$r1 eq 1} {return 1}
	array set trans {X 100 Y 101 M 102}
	global cache
	set chr [lindex $cur2 1]
	set chr [get trans($chr) $chr]
	set pos [lindex $cur2 2]
	set line $cache($r1)
	while {[llength $line]} {
		foreach {rchr rstart rend} $line break
		set rchr [get trans($rchr) $rchr]
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

proc readonevar {f} {
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {[llength $line]} break
	}
	return $line
}

array set compare_annot_join_trans {
	COMPATIBLE 1 MISSENSE 2 DELETE 3 INSERT 4 DELETE+ 5 INSERT+ 6 NONSTOP 7 NONSENSE 8 FRAMESHIFT 9
}

proc compare_annot_join {compar sample cur1 cur2} {
	global compare_annot_join_trans
	set id [lindex $cur1 0]
	set type [lindex $cur1 4]
	set a1 [lindex $cur1 5]
	set a2 [lindex $cur1 6]
	set oloc [lindex $cur1 14]
	set oeffect [lindex $cur1 17]
	set xref [lindex $cur1 9]
	set score1 [lindex $cur1 7]
	set score2 [lindex $cur1 8]
	set hz {}
	if {[llength $cur1] < 23} {
		lset_always cur1 22 {}
	}
	if {![llength $cur2]} {
		if {$a1 eq $a2} {set hz hz}
		set cur2 {{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}}
		set type2 $type
		set a12 $a1
		set a22 $a2
		set oloc2 $oloc
		set oeffect2 $oeffect
		set xref2 $xref
		set score12 $score1
		set score22 $score2
	} else {
		if {[llength $cur2] < 23} {
			lset_always cur2 22 {}
		}
		set type2 [lindex $cur2 4]
		set a12 [lindex $cur2 5]
		set a22 [lindex $cur2 6]
		set oloc2 [lindex $cur2 14]
		set oeffect2 [lindex $cur2 17]
		set xref2 [lindex $cur2 9]
		set score12 [lindex $cur2 7]
		set score22 [lindex $cur2 8]
	}
	if {[regexp dbsnp $xref]||[regexp dbsnp $xref2]} {set dbsnp dbsnp} else {set dbsnp ""}
	if {[regexp {\?} $a1$a2$a12$a22] || [regexp {N} $a1$a2$a12$a22]} {set ns n?} else {set ns ""}
	set minscore [lmath_min [list_remove [list $score1 $score2 $score12 $score22] {}]]
	if {$minscore < 60} {set lowscore ls} else {set lowscore ""}
	if {$type eq "snp"} {
		set snp 1
	} else {
		set snp 0
	}
	set loc ""
	foreach temp {EXON UTR BEGIN END INTRON} {
		if {[regexp $temp $oloc]} {
			set loc $temp
			break
		}
	}
	set effect OTHER
	foreach temp {FRAMESHIFT NONSENSE NONSTOP INSERT+ DELETE+ INSERT DELETE MISSENSE COMPATIBLE UNDEFINED} {
		if {[regexp $temp $oeffect]} {
			set effect $temp
			break
		}
		if {[regexp $temp $oeffect2]} {
			set effect $temp
			break
		}
	}
	set neffect [get compare_annot_join_trans($effect) 0]
	# set cur1 {locus chromosome begin end type alleleSeq1 alleleSeq2 totalScore1 totalScore2 xRef geneId mrnaAcc proteinAcc orientation exonCategory exon codingRegionKnown aaCategory nucleotidePos proteinPos aaAnnot aaCall aaRef}
	# set cur2 {locus2 chromosome2 begin2 end2 type2 alleleSeq12 alleleSeq22 totalScore12 totalScore22 xRef2 geneId2 mrnaAcc2 proteinAcc2 orientation2 exonCategory2 exon2 codingRegionKnown2 aaCategory2 nucleotidePos2 proteinPos2 aaAnnot2 aaCall2 aaRef2}
	set result [list $compar $sample]
	list_append result [lrange $cur1 0 8]
	list_append result [lrange $cur2 5 8]
	lappend result $dbsnp $ns $lowscore $minscore $loc $effect $neffect $hz
	list_append result [lrange $cur1 9 22]
	list_append result [lrange $cur2 9 22]
	return [join $result \t]
}

proc compare_annot {id1 file1 regfile1 id2 file2 regfile2 outfile} {
	global cache
	set f1 [open $file1]
	set header [split [gets $f1] \t]
	set annotvarfields {locus chromosome begin end type alleleSeq1 alleleSeq2 totalScore1 totalScore2 xRef geneId mrnaAcc proteinAcc orientation exonCategory exon codingRegionKnown aaCategory nucleotidePos proteinPos aaAnnot aaCall aaRef}
	if {$header ne $annotvarfields} {
		puts stderr "header error in annot_varfile1 $file1"
		exit 1
	}
	set f2 [open $file2]
	set header [split [gets $f2] \t]
	if {$header ne $annotvarfields} {
		puts stderr "header error in annot_varfile2 $file1"
		exit 1
	}
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
	puts $o [join {
		compar sample
		locus chromosome begin end type
		alleleSeq1 alleleSeq2 totalScore1 totalScore2 
		alleleSeq1-2 alleleSeq2-2 totalScore1-2 totalScore2-2 
		dbsnp ns lowscore minscore loc effect neffect hz
		xRef geneId mrnaAcc proteinAcc orientation exonCategory exon codingRegionKnown aaCategory nucleotidePos proteinPos aaAnnot aaCall aaRef
		xRef2 geneId2 mrnaAcc2 proteinAcc2 orientation2 exonCategory2 exon2 codingRegionKnown2 aaCategory2 nucleotidePos2 proteinPos2 aaAnnot2 aaCall2 aaRef2
	} \t]
	set cur1 [readonevar $f1]
	set cur2 [readonevar $f2]
	set num 1
	while {![eof $f1] || ![eof $f2]} {
		incr num
		if {![expr {$num % 100000}]} {puts stderr $num}
		set d [comparepos $cur1 $cur2]
		if {$d == 0} {
			set s1 [sequenced $r1 $cur1]
			set s2 [sequenced $r2 $cur2]
			if {!$s1 || !$s2} {
				puts $o [compare_annot_join fl $id1,$id2 $cur1 $cur2]
			} else {
				set type [list_remove [split [lindex $cur1 4] _] ref-consistent ref-inconsistent =]
				lset cur1 4 $type
				set type [list_remove [split [lindex $cur2 4] _] ref-consistent ref-inconsistent =]
				lset cur2 4 $type
				if {([lrange $cur1 1 4] eq [lrange $cur2 1 4]) && (
					([list_sub $cur1 {5 6}] eq [list_sub $cur2 {5 6}])
					|| ([list_sub $cur1 {6 5}] eq [list_sub $cur2 {5 6}])
					)} {
					puts $o [compare_annot_join sm $id1,$id2 $cur1 $cur2]
				} else {
					puts $o [compare_annot_join mm $id1,$id2 $cur1 $cur2]
				}
			}
			set cur1 [readonevar $f1]
			set cur2 [readonevar $f2]
			continue
		} elseif {$d < 0} {
			while {[comparepos $cur1 $cur2] < 0} {
				set s [sequenced $r1 $cur1]
				if {!$s} {
					puts $o [compare_annot_join fl $id1 $cur1 {}]
				} else {
					set s [sequenced $r2 $cur1]
					if {$s} {
						puts $o [compare_annot_join df $id1 $cur1 {}]
					} else {
						puts $o [compare_annot_join un $id1 $cur1 {}]
					}
				}
				set cur1 [readonevar $f1]
			}
			continue
		} else {
			while {[comparepos $cur1 $cur2] > 0} {
				set s [sequenced $r2 $cur2]
				if {!$s} {
					puts $o [compare_annot_join fl $id2 $cur2 {}]
				} else {
					set s [sequenced $r1 $cur2]
					if {$s} {
						puts $o [compare_annot_join df $id2 $cur2 {}]
					} else {
						puts $o [compare_annot_join un $id2 $cur2 {}]
					}
				}
				set cur2 [readonevar $f2]
			}
		}
	}
	close $o
}

proc lset_always {varName pos value} {
	upvar $varName var
	if {[catch {lset var $pos $value}]} {
		set var [list_concat $var [list_fill [expr {$pos - [llength $var]}] {}]]
		lappend var $value
	}
}

array set chrtrans {X 90 Y 91 M 92}
proc chr2num {chr} {
	regsub ^chr $chr {} chr
	set nchr [get ::chrtrans($chr) $chr]
	return $nchr
}

proc annot_compare_region {compar_file reg_file field tvalue fvalue} {
	set f1 [open $compar_file]
	set header1 [gets $f1]
	set poss1 [list_cor $header1 {chromosome begin end}]
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
		if {![expr $num%100000]} {puts stderr $num}
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
			set line2 [split [gets $f2] \t]
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
		set line1 [split [gets $f1] \t]
		if {![llength $line1]} break
		foreach {chr1 start1 end1} [list_sub $line1 $poss1] break
		set nchr1 [chr2num $chr1]
	}
	close $f1
	close $f2
}

if 0 {

lappend auto_path ~/dev/completegenomics/lib
package require Extral
package require Tclx
signal -restart error SIGINT

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
