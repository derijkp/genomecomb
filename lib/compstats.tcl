#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc compare_pvt {} {
	set f stdin
	# set f [open $compar_file]
	unset -nocomplain a
	set header [split [gets $f] \t]
	set fields {compar sample chromosome type trf str repeat segdup selfchain}
	set poss [list_cor $header $fields]
	set ufields {
		type xRef
		alleleSeq1-1 alleleSeq2-1 alleleSeq1-2 alleleSeq2-2
		totalScore1-1 totalScore2-1 totalScore1-2 totalScore2-2
		refscore-1 coverage-1 refcons-1 cluster-1
		refscore-2 coverage-2 refcons-2 cluster-2
		exonCategory aaCategory
	}
	set uposs [list_cor $header $ufields]

	# set fields {compar sample chromosome type dbsnp ns lowscore loc effect refcons cluster trf str repeat segdup selfchain a100 a70 b15 b20 b30}
	set nf [list_find $poss -1]
	if {[llength $nf]} {
		set fields [list_sub $fields -exclude $nf]
		set poss [list_sub $poss -exclude $nf]
	}
	set num 0
	while {![eof $f]} {
		incr num
		if {![expr {$num%100000}]} {putsprogress $num}
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		foreach $ufields [list_sub $line $uposs] break
		set line [list_sub $line $poss]
		# dbsnp
		if {[regexp dbsnp $xRef]} {set dbsnp dbsnp} else {set dbsnp {}}
		# ns
		set alleles [list ${alleleSeq1-1} ${alleleSeq2-1} ${alleleSeq1-2} ${alleleSeq2-2}]
		if {[regexp {[N?]} $alleles]} {set ns ns} else {set ns ""}
		# lowscore
		set scores [list_remove [list ${totalScore1-1} ${totalScore2-1} ${totalScore1-2} ${totalScore2-2}] {} - ?]
		set minscore [lmath_min $scores]
		if {$minscore < 60} {set lowscore ls} else {set lowscore ""}
		# loc
		set loc ""
		foreach temp {EXON UTR BEGIN END INTRON} {
			if {[regexp $temp $exonCategory]} {
				set loc $temp
				break
			}
		}
		# effect
		set effect OTHER
		foreach temp {FRAMESHIFT NONSENSE NONSTOP INSERT+ DELETE+ INSERT DELETE MISSENSE COMPATIBLE UNDEFINED} {
			if {[regexp $temp $aaCategory]} {
				set effect $temp
				break
			}
		}
		set refcons [lindex [list_remove [list ${refcons-1} ${refcons-2}] {} - ?] 0]
		set cluster [lindex [list_remove [list ${cluster-1} ${cluster-2}] {} - ?] 0]
		set coverages [lsort -integer [list_remove [list ${coverage-1} ${coverage-2}] {} - ?]]
		set mincov [lindex $coverages 0]
		set maxcov [lindex $coverages end]
		if {$mincov < 15} {set b15 b15} else {set b15 {}}
		if {$mincov < 20} {set b20 b20} else {set b20 {}}
		if {$mincov < 30} {set b30 b30} else {set b30 {}}
		if {$maxcov > 70} {set a70 a70} else {set a70 {}}
		if {$maxcov > 100} {set a100 a100} else {set a100 {}}
		lappend line $dbsnp $ns $lowscore $loc $effect $refcons $cluster $a100 $a70 $b15 $b20 $b30
		if {![info exists a($line)]} {
			set a($line) 1
		} else {
			incr a($line)
		}
	}
	set o stdout
	# set o [open pvtcompar_${name1}_${name2}.tsv w]
	set oheader $fields
	lappend oheader dbsnp ns lowscore loc effect refcons cluster a100 a70 b15 b20 b30 tnum
	puts $o [join $oheader \t]
	set list [lsort -dict [array names a]]
	foreach line $list {
		puts $o [join $line \t]\t$a($line)
	}
}

proc comstats_filter {slist cfilter} {
	global header
	if {![llength $cfilter]} {
		return $slist
	}
	set fields [list_unmerge $cfilter 1 values]
	set poss [list_cor $header $fields]
	set result {}
	foreach line $slist {
		set add 1
		foreach lval [list_sub $line $poss] val $values {
			if {![inlist $val $lval] && ($val ne $lval)} {
				set add 0
				break
			}
		}
		if {$add} {
			lappend result $line
		}
	}
	return $result
}

proc count {slist {cfilter {}}} {
	if {[llength $cfilter]} {
		set slist [comstats_filter $slist $cfilter]
	}
	if {![llength $slist]} {return 0}
	set end [expr {[llength [lindex $slist 0]]-1}]
	set result [format %0.0f [lmath_sum [list_subindex $slist $end]]]
	return $result
}

proc compare_loadtable {file headerVar} {
	upvar $headerVar header
	if {$file ne "stdin"} {
		set f [gzopen $file]
	} else {
		set f stdin
	}
	set header [split [gets $f] \t]
	set num 0
	set table {}
	while {![eof $f]} {
		incr num
		if {![expr {$num%100000}]} {putsprogress $num}
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		lappend table $line
	}
	if {$file ne "stdin"} {
		close $f
	}
	return $table
}

proc compare_pvtsummary_out1 {o title header flist} {
	set loclist {EXON UTR BEGIN END INTRON OTHER ""}
	set effectlist {FRAMESHIFT NONSENSE NONSTOP INSERT+ DELETE+ INSERT DELETE MISSENSE COMPATIBLE UNDEFINED OTHER}
	set pos [lsearch $header sample]
	set idlist [list_remdup [split [list_remdup [list_subindex $flist $pos]] " ,"]]
	set len [llength $header]
	set end [expr {$len-1}]
	puts $o "\n ==================== $title ===================="
	foreach compar {df mm sm} {
		putslog $compar
		puts $o "\n---------- $compar: $title ----------"
		set list [comstats_filter $flist [list compar $compar]]
		set checkall [format %0.0f [lmath_sum [list_subindex $list $end]]]
		puts $checkall
		set all [count $list]
		set allnew [count $list {dbsnp ""}]
		puts $o "total:\t$all"
		set snps [count $list {type snp}]
		set snpsnew [count $list {type snp dbsnp ""}]
		puts $o "snps:\t$snps ([expr {$all-$snps}] other)\tnewsnps\t$snpsnew ([expr {$allnew-$snpsnew}] other)"
		if {$compar eq "df"} {
			set idlistheader \t\t[join $idlist \t]
		} else {
			set idlistheader ""
		}
		puts $o "\nloc\ttotal\tnew\tdbsnp$idlistheader"
		foreach l $loclist {
			set temp1 [count $list [list dbsnp "" loc $l]]
			set temp2 [count $list [list dbsnp dbsnp loc $l]]
			set oline "$l:\t[expr {$temp1+$temp2}]\t${temp1}\t${temp2}\t"
			if {$compar eq "df"} {
				foreach id $idlist {
					append oline \t[count $list [list loc $l sample $id]]
				}
			}
			puts $o $oline
		}
		set temp "\neffect\ttotal\tnew\tdbsnp\t"
		if {$compar eq "df"} {
			foreach id $idlist {
				append temp \t$id\t$id\(dbsnp\)
			}
		}
		puts $o $temp
		foreach e $effectlist {
			set temp1 [count $list [list dbsnp "" loc EXON effect $e]]
			set temp2 [count $list [list dbsnp dbsnp loc EXON effect $e]]
			set oline "$e\t[expr {$temp1+$temp2}]\t$temp1\t$temp2\t"
			if {$compar eq "df"} {
				foreach id $idlist {
					append oline \t[count $list [list dbsnp "" loc EXON effect $e sample $id]]
					append oline \t[count $list [list dbsnp dbsnp loc EXON effect $e sample $id]]
				}
			}
			puts $o $oline
		}
		set temp1 [count $list [list dbsnp "" loc EXON effect {FRAMESHIFT NONSENSE NONSTOP INSERT+ DELETE+ INSERT DELETE}]]
		set temp2 [count $list [list dbsnp dbsnp loc EXON effect {FRAMESHIFT NONSENSE NONSTOP INSERT+ DELETE+ INSERT DELETE}]]
		set oline "*BAD*\t[expr {$temp1+$temp2}]\t$temp1\t$temp2\t"
		if {$compar eq "df"} {
			foreach id $idlist {
				append oline \t[count $list [list dbsnp "" loc EXON effect {FRAMESHIFT NONSENSE NONSTOP INSERT+ DELETE+ INSERT DELETE} sample $id]]
				append oline \t[count $list [list dbsnp dbsnp loc EXON effect {FRAMESHIFT NONSENSE NONSTOP INSERT+ DELETE+ INSERT DELETE} sample $id]]
			}
		}
		puts $o $oline
	}
}

proc compare_pvtsummary_out2 {o title header flist} {
	set len [llength $header]
	set end [expr {$len-1}]
	foreach compar {df mm sm un} {
		set list [comstats_filter $flist [list compar $compar]]
		set checkall [format %0.0f [lmath_sum [list_subindex $list $end]]]
		set all [count $list]
		set allnew [count $list {dbsnp ""}]
		set snps [count $list {type snp}]
		set indels [expr {$all-$snps}]
		set snpsnew [count $list {type snp dbsnp ""}]
		lappend oline $all $allnew $snps $indels $snpsnew
	}
	puts $o $title\t[join $oline \t]
}

proc compare_pvtsummary {{type 1}} {
	global header
	set file stdin
	set o stdout
	# set o [open /complgen/multicompar/summary_GS102_GS103.tsv w]
	set table [compare_loadtable $file header]
	# selectivity
	if {$type == 1} {
		compare_selectivity $o $table
	}
	# filters
	array set trans {
		COMPATIBLE 1 MISSENSE 2 DELETE 3 INSERT 4 DELETE+ 5 INSERT+ 6 NONSTOP 7 NONSENSE 8 FRAMESHIFT 9
	}
	# set o [open $ofile a]
	set o stdout
	set filterlist {
		nofilter {}
		"filter 1 (removed ref-(in)consistent)" {refcons ""}
	}
	set ctitle {removed ref-(in)consistent}
	set cfilter {refcons ""}
	foreach {name q} {
		{?/N} {ns ""}
		{low scores (< 60 for cg, posterior < 4 for rtg)} {lowscore ""}
		cluster {cluster ""}
		{tandem repeats} {trf ""}
		microsatelite {str ""}
		{coverage between 20 and 100} {a100 "" b20 ""}
		segdup {segdup ""}
		selfchain {selfchain ""}
		repeat {repeat ""}
	} {
		append ctitle ", $name"
		eval lappend cfilter $q
		lappend filterlist $ctitle $cfilter
	}

	if {$type == 2} {
		set wheader filter
		foreach compar {df mm sm un} {
			lappend wheader ${compar}-all ${compar}-allnew ${compar}-snps ${compar}-indels ${compar}-snpsnew
		}
		puts $o [join $wheader \t]
	}
	foreach {title filters} $filterlist {
		putslog "===== $title ====="
		set flist [comstats_filter $table $filters]
		compare_pvtsummary_out$type $o $title $header $flist
		flush $o
	}
}

proc compare_selectivity {o table} {
	global header
	set len [llength [lindex $table 0]]
	set end [expr {$len-1}]
	set dflist [comstats_filter $table [list compar df]]
	set mmlist [comstats_filter $table [list compar mm]]
	set smlist [comstats_filter $table [list compar sm]]
	set dftotal [lmath_sum [list_subindex $dflist $end]]
	set mmtotal [lmath_sum [list_subindex $mmlist $end]]
	set smtotal [lmath_sum [list_subindex $smlist $end]]
	set newdflist [comstats_filter $dflist [list dbsnp ""]]
	set newmmlist [comstats_filter $mmlist [list dbsnp ""]]
	set newsmlist [comstats_filter $smlist [list dbsnp ""]]
	set newdftotal [lmath_sum [list_subindex $newdflist $end]]
	set newmmtotal [lmath_sum [list_subindex $newmmlist $end]]
	set newsmtotal [lmath_sum [list_subindex $newsmlist $end]]
	puts $o "---------- total ----------"
	puts $o "df\t[format %0.0f $dftotal]\tnewdf\t[format %0.0f $newdftotal]"
	puts $o "mm\t[format %0.0f $mmtotal]\tnewmm\t[format %0.0f $newmmtotal]"
	puts $o "sm\t[format %0.0f $smtotal]\tnewsm\t[format %0.0f $newsmtotal]"
	foreach {sfield svalue} {
		type snp dbsnp dbsnp ns ns lowscore ls 
		refcons rc trf trf str str repeat rp segdup sd selfchain sc cluster cl
		a100 a100 b20 b20
	} {
		set list [comstats_filter $dflist [list $sfield $svalue]]
		set dfselnum [lmath_sum [list_subindex $list $end]]
		set list [comstats_filter $mmlist [list $sfield $svalue]]
		set mmselnum [lmath_sum [list_subindex $list $end]]
		set list [comstats_filter $smlist [list $sfield $svalue]]
		set smselnum [lmath_sum [list_subindex $list $end]]
		#new
		set list [comstats_filter $newdflist [list $sfield $svalue]]
		set newdfselnum [lmath_sum [list_subindex $list $end]]
		set list [comstats_filter $newmmlist [list $sfield $svalue]]
		set newmmselnum [lmath_sum [list_subindex $list $end]]
		set list [comstats_filter $newsmlist [list $sfield $svalue]]
		set newsmselnum [lmath_sum [list_subindex $list $end]]
		puts $o "---------- $sfield == $svalue ----------"
		puts $o "df\t[format %0.0f $dfselnum]\t[format %.2f [expr {100.0*$dfselnum/$dftotal}]]\tnewdf\t[format %0.0f $newdfselnum]\t[format %.2f [expr {100.0*$newdfselnum/$newdftotal}]]"
		puts $o "mm\t[format %0.0f $mmselnum]\t[format %.2f [expr {100.0*$mmselnum/$mmtotal}]]\tnewmm\t[format %0.0f $newmmselnum]\t[format %.2f [expr {100.0*$newmmselnum/$newmmtotal}]]"
		puts $o "sm\t[format %0.0f $smselnum]\t[format %.2f [expr {100.0*$smselnum/$smtotal}]]\tnewsm\t[format %0.0f $newsmselnum]\t[format %.2f [expr {100.0*$newsmselnum/$newsmtotal}]]"
	}
}

if 0 {
}

if 0 {


package require Tclx
signal -restart error SIGINT
lappend auto_path ~/dev/completegenomics/lib
package require Extral

set compar_file /complgen/multicompar/compar.tsv.rz
set file /complgen/multicompar/pvtcompar_GS102_GS103.tsv
set sample1 GS102
set sample2 GS103

cd /complgen/compar_GS102_GS103
set file pvtcompar_GS102_GS103.tsv
set ofile summarycompar2_GS102_GS103.tsv
}
