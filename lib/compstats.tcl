proc count {var filters snps dbsnps locs effects ids} {
	upvar $var a
	set result 0
	foreach filter $filters {
		foreach snp $snps {
			foreach dbsnp $dbsnps {
				foreach loc $locs {
					foreach effect $effects {
						foreach id $ids {
							set code $filter,$snp,$dbsnp,$loc,$effect,$id
							incr result [get a($code) 0]
						}
					}
				}
			}
		}
	}
	return $result
}

proc comstats {ofile file} {
	array set trans {
		COMPATIBLE 1 MISSENSE 2 DELETE 3 INSERT 4 DELETE+ 5 INSERT+ 6 NONSTOP 7 NONSENSE 8 FRAMESHIFT 9
	}
	unset -nocomplain a
	set filtered {}
	set f [open $file]
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set id [lindex $line 0]
		set ida($id) 1
		set type [lindex $line 5]
		set a1 [lindex $line 6]
		set a2 [lindex $line 7]
		set oloc [lindex $line 15]
		set oeffect [lindex $line 18]
		set xref [lindex $line 10]
		set score1 [lindex $line 8]
		set score2 [lindex $line 9]
		if {[llength $line] > 25} {
			set type2 [lindex $line 29]
			set a12 [lindex $line 30]
			set a22 [lindex $line 31]
			set oloc2 [lindex $line 39]
			set oeffect2 [lindex $line 42]
			set xref2 [lindex $line 34]
			set score12 [lindex $line 32]
			set score22 [lindex $line 33]
		} else {
			set type2 $type
			set a12 $a1
			set a22 $a2
			set oloc2 $oloc
			set oeffect2 $oeffect
			set xref2 $xref
			set score12 $score1
			set score22 $score2
		}
		if {[regexp dbsnp $xref]} {set dbsnp 1} else {set dbsnp 0}
		if {[regexp dbsnp $xref2]} {set dbsnp2 1} else {set dbsnp2 0}
		set filter 0
		if {![regexp {ref-inconsistent} $type] && ![regexp {ref-inconsistent} $type2]} {
			set filter 1
			if {![regexp {ref-consistent} $type] && ![regexp {ref-consistent} $type2]} {
				set filter 2
				if {![regexp {\?} $a1$a2] && ![regexp {N} $a1$a2] && ![regexp {\?} $a12$a22] && ![regexp {N} $a12$a22]} {
					set filter 3
					if {([isint $score1] && ($score1 >= 60)) && ([isint $score2] && ($score2 >= 60))
						&&([isint $score12] && ($score12 >= 60)) && ([isint $score22] && ($score22 >= 60))} {
						set filter 4
					}
				}
			}
		}
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
		set neffect [get trans($effect) 0]
		set code $filter,$snp,$dbsnp,$loc,$effect,$id
		if {![info exists a($code)]} {
			set a($code) 1
		} else {
			incr a($code)
		}
#		if {($filter == 4) && ($loc eq "EXON") && !$dbsnp} {
#			lappend filtered [list_concat $neffect $line]
#		}
		if {($filter == 4) && ([inlist {FRAMESHIFT NONSENSE NONSTOP INSERT+ DELETE+ INSERT DELETE MISSENSE} $effect])
			&& (!$dbsnp || !$dbsnp2)
			&& [regexp 79 $id] 
			&& ($type ne "delins") && ($type2 ne "delins")} {
			lappend filtered [list_concat $neffect $line]
		}
	}
	close $f
	set idlist [array names ida]
	set filterlist {0 1 2 3}
	set loclist {EXON UTR BEGIN END INTRON OTHER ""}
	set effectlist {FRAMESHIFT NONSENSE NONSTOP INSERT+ DELETE+ INSERT DELETE MISSENSE COMPATIBLE UNDEFINED OTHER}
	set o [open $ofile a]
	puts $o "\n------------------------- $file -------------------------"
	set checkall 0
	foreach name [array names a] {incr checkall $a($name)}
	puts $checkall
	foreach {title filters} {
		nofilter {0 1 2 3 4}
		"filter 1 (removed ref-inconsistent)" {1 2 3 4}
		"filter 2 (removed ref-inconsistent and ref-consistent)" {2 3 4}
		"filter 3 (removed ref-inconsistent and ref-consistent and ?/N in alleles)" {3 4}
		"filter 4 (removed ref-inconsistent and ref-consistent and ?/N in alleles and scores < 60)" {4}
	} {
		puts $o "\n ----- $file: $title -----"
		set all [count a $filters {0 1} {0 1} $loclist $effectlist $idlist]
		puts $o "total:\t$all"
		set snps [count a $filters 1 {0 1} $loclist $effectlist $idlist]
		puts $o "snps:\t$snps ([expr {$all-$snps}] other)"
		puts $o "\nloc\ttotal\tnew\tdbsnp\t\t[join $idlist \t]"
		foreach l $loclist {
			set temp1 [count a $filters {0 1} 0 $l $effectlist $idlist]
			set temp2 [count a $filters {0 1} 1 $l $effectlist $idlist]
			set oline "$l:\t[expr {$temp1+$temp2}]\t${temp1}\t${temp2}\t"
			foreach id $idlist {
				append oline \t[count a $filters {0 1} {0 1} $l $effectlist $id]
			}
			puts $o $oline
		}
		set temp "\neffect\ttotal\tnew\tdbsnp\t"
		foreach id $idlist {
			append temp \t$id\t$id\(dbsnp\)
		}
		puts $o $temp
		foreach e $effectlist {
			set temp1 [count a $filters {0 1} 0 EXON $e $idlist]
			set temp2 [count a $filters {0 1} 1 EXON $e $idlist]
			set oline "$e\t[expr {$temp1+$temp2}]\t$temp1\t$temp2\t"
			foreach id $idlist {
				append oline \t[count a $filters {0 1} 0 EXON $e $id]
				append oline \t[count a $filters {0 1} 1 EXON $e $id]
			}
			puts $o $oline
		}
		set temp1 [count a $filters {0 1} 0 EXON {FRAMESHIFT NONSENSE NONSTOP INSERT+ DELETE+ INSERT DELETE} $idlist]
		set temp2 [count a $filters {0 1} 1 EXON {FRAMESHIFT NONSENSE NONSTOP INSERT+ DELETE+ INSERT DELETE} $idlist]
		set oline "*BAD*\t[expr {$temp1+$temp2}]\t$temp1\t$temp2\t"
		foreach id $idlist {
			append oline \t[count a $filters {0 1} 0 EXON {FRAMESHIFT NONSENSE NONSTOP INSERT+ DELETE+ INSERT DELETE} $id]
			append oline \t[count a $filters {0 1} 1 EXON {FRAMESHIFT NONSENSE NONSTOP INSERT+ DELETE+ INSERT DELETE} $id]
		}
		puts $o $oline
	}
	close $o
	set o [open [file root $file]-filtered.tsv w]
	set len [llength [lindex $filtered 0]]
	set header {locus chromosome begin end type alleleSeq1 alleleSeq2 totalScore1 totalScore2 xRef geneId mrnaAcc proteinAcc orientation exonCategory exon codingRegionKnown aaCategory nucleotidePos proteinPos aaAnnot aaCall aaRef}
	if {$len > 26} {
		set header [list_concat $header id2 $header]
	}	
	puts $o effect\tid\t[join $header \t]
	foreach line [lsort -integer -index 0 -decreasing $filtered] {
		puts $o [join $line \t]
	}
	close $o
}

if 0 {




lappend auto_path ~/dev/completegenomics/lib
package require Extral

cd /complgen/compar
set file 78vs79_mismatch.tsv
set ofile summary.tsv
catch {file delete summary.tsv}
comstats summary.tsv diff.tsv
comstats summary.tsv same.tsv
comstats summary.tsv uniqueseq.tsv

}
