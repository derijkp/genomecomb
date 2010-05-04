# < 0 if comp1 < comp2
# > 0 if comp1 > comp2
proc comparepos {comp1 comp2} {
	if {![llength $comp1]} {return 1}
	if {![llength $comp2]} {return -1}
	foreach {chr1 pos1 end1 type1} $comp1 break
	foreach {chr2 pos2 end2 type2} $comp2 break
	if {$chr1 ne $chr2} {
		set chr1 [chr2num $chr1]
		set chr2 [chr2num $chr2]
		return [expr {$chr1-$chr2}]
	} elseif {$pos1 != $pos2} {
		return [expr {$pos1-$pos2}]
	} elseif {$end1 != $end2} {
		return [expr {$end1-$end2}]
	} elseif {$type1 ne $type2} {
		if {$type1 < $type2} {return -1} else {return 1}
	} else {
		return 0
	}
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

	catch {close $f1}; catch {close $f2}; catch {annot_region_close $regfile1}; catch {annot_region_close $regfile2}; catch {close $o}
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
		annot_region_init $regfile1
	} else {
		puts "no region file for file1"
	}
	if {$regfile2 ne ""} {
		annot_region_init $regfile2
	} else {
		puts "no region file for file2"
	}
	# start
	set o [open $outfile.temp w]
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
			set s1 [annot_region_in $regfile1 {*}[lrange $comp1 0 2]]
			set s2 [annot_region_in $regfile2 {*}[lrange $comp2 0 2]]
			if {!$s1 || !$s2} {
				error "variant \"$cur1\" and \"$cur2\" not sequenced"
				# puts $o [compare_annot_join fl $id1,$id2 $cur1 $cur2]
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
				set s [annot_region_in $regfile1 {*}[lrange $comp1 0 2]]
				if {!$s} {
					error "variant1 \"$cur1\" not sequenced"
					puts $o [compare_annot_join fl $id1 $cur1 -]
				} else {
					set s [annot_region_in $regfile2 {*}[lrange $comp1 0 2]]
					if {$s} {
						puts $o [compare_annot_join df $id1 $cur1 {}]
					} else {
						puts $o [compare_annot_join un $id1 $cur1 -]
					}
				}
				if {[eof $f1]} break
				set cur1 [compare_annot_getline $f1]
				set comp1 [list_sub $cur1 $comparposs1]
				if {![llength $cur1]} break
			}
		} else {
			while {[comparepos $comp1 $comp2] > 0} {
				set s [annot_region_in $regfile2 {*}[lrange $comp2 0 2]]
				if {!$s} {
					error "variant2 \"$cur2\" not sequenced"
					puts $o [compare_annot_join fl $id2 - $cur2]
				} else {
					set s [annot_region_in $regfile1 {*}[lrange $comp2 0 2]]
					if {$s} {
						puts $o [compare_annot_join df $id2 {} $cur2]
					} else {
						puts $o [compare_annot_join un $id2 - $cur2]
					}
				}
				if {[eof $f2]} break
				set cur2 [compare_annot_getline $f2]
				set comp2 [list_sub $cur2 $comparposs2]
				if {![llength $cur2]} break
			}
		}
	}

	annot_region_close $regfile1
	annot_region_close $regfile2
	close $f1; close $f2; close $o
	file rename $outfile.temp $outfile
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
	set o [open $outfile.temp w]
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
	file rename $outfile.temp $outfile
}

if 0 {

	lappend auto_path ~/dev/completegenomics/lib
	lappend auto_path ~/bin/complgen/apps/cg/lib
	package require Extral
	package require Tclx
	signal -restart error SIGINT
	
	set basedir /media/passport/complgen
	set basedir /complgen
	set dbdir /complgen/refseq
	set dir1 $basedir/GS101A01
	set dir2 $basedir/GS101A02
	set resultsdir $basedir/testcompar_GS101A01_GS101A02
	set force 0
	set name1 [file tail $dir1]
	set id1 $name1
	set id2 $name2
	set name2 [file tail $dir2]
	set file1 $dir1/fannotvar-$name1.tsv
	set file2 $dir2/fannotvar-$name2.tsv
	set regfile1 $dir1/sreg-$name1.tsv
	set regfile2 $dir2/sreg-$name2.tsv
	set outfile $resultsdir/compar_${name1}_${name2}.tsv
	file mkdir $resultsdir
	cd $resultsdir
	
	
	compare_annot $id1 $file1 $regfile1 $id2 $file2 $regfile2 $outfile
	
	reannot_compare $compar_file $dir1 $dir2 atemp
	
}
