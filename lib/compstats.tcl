proc compare_pvt {} {
	set f stdin
	# set f [open $compar_file]
	unset -nocomplain a
	set header [split [gets $f] \t]
	set fields {compar sample chromosome type dbsnp ns lowscore loc effect refcons cluster trf str repeat segdup selfchain a100 a70 b20 b30}
	set poss [list_cor $header $fields]
	set nf [list_find $poss -1]
	if {[llength $nf]} {
		set fields [list_sub $fields -exclude $nf]
		set poss [list_sub $poss -exclude $nf]
	}
	set num 0
	while {![eof $f]} {
		incr num
		if {![expr {$num%100000}]} {puts stderr $num}
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set line [list_sub $line $poss]
		if {![info exists a($line)]} {
			set a($line) 1
		} else {
			incr a($line)
		}
	}
	set o stdout
	puts $o [join $fields \t]\tnum
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

proc compare_pvtsummary {} {
	global header
	set f stdin
	# set f [open $file]
	set header [split [gets $f] \t]
	set num 0
	set table {}
	while {![eof $f]} {
		incr num
		if {![expr {$num%100000}]} {puts stderr $num}
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		lappend table $line
	}
	set len [llength $header]
	set end [expr {$len-1}]
	set o stdout
	# selectivity
	compare_selectivity $o $table
	# filters
	array set trans {
		COMPATIBLE 1 MISSENSE 2 DELETE 3 INSERT 4 DELETE+ 5 INSERT+ 6 NONSTOP 7 NONSENSE 8 FRAMESHIFT 9
	}
	set pos [lsearch $header sample]
	set idlist [list_remdup [split [list_remdup [list_subindex $table $pos]] " ,"]]
	set filterlist {0 1 2 3}
	set loclist {EXON UTR BEGIN END INTRON OTHER ""}
	set effectlist {FRAMESHIFT NONSENSE NONSTOP INSERT+ DELETE+ INSERT DELETE MISSENSE COMPATIBLE UNDEFINED OTHER}
	# set o [open $ofile a]
	set o stdout
	foreach {title filters} {
		nofilter {}
		"filter 1 (removed ref-(in)consistent)" {refcons ""}
		"filter 3 (removed ref-(in)consistent, ?/N)" {refcons "" ns ""}
		"filter 4 (removed ref-(in)consistent, ?/N,  scores < 60)" {refcons "" ns "" lowscore ""}
		"filter 5 (removed ref-(in)consistent, ?/N, scores < 60, tandem repeats)" {refcons "" ns "" lowscore "" trf ""}
		"filter 5b (removed ref-(in)consistent, ?/N, coverage between 20 and 100, tandem repeats)" {refcons "" ns "" a100 "" b20 "" trf ""}
		"filter 6 (removed ref-(in)consistent, ?/N, scores < 60, tandem repeats, str)" {refcons "" ns "" lowscore "" trf "" str ""}
		"filter 7 (removed ref-(in)consistent, ?/N, scores < 60, tandem repeats, str, repeats)" {refcons "" ns "" lowscore "" trf "" str "" repeat ""}
		"filter 8 (removed ref-(in)consistent, ?/N, scores < 60, tandem repeats, str, repeats, segdup)" {refcons "" ns "" lowscore "" trf "" str "" repeat "" segdup ""}
		"filter 9 (removed ref-(in)consistent, ?/N, scores < 60, tandem repeats, str, repeats, segdup, self chained)" {refcons "" ns "" lowscore "" trf "" str "" repeat "" segdup "" selfchain ""}
		"filter 9b (removed ref-(in)consistent, ?/N, coverage between 20 and 100, tandem repeats, str, repeats, segdup, self chained)" {refcons "" ns "" a100 "" b20 "" trf "" str "" repeat "" segdup "" selfchain ""}
		"filter 10 (removed ref-(in)consistent, ?/N, scores < 60, coverage between 20 and 100, tandem repeats, str, repeats, segdup, self chained)" {refcons "" ns "" lowscore "" a100 "" b20 "" trf "" str "" repeat "" segdup "" selfchain ""}
		"filter 11 (removed ref-(in)consistent, ?/N, scores < 60, coverage between 20 and 100, tandem repeats, str, repeats, segdup, self chained, cluster)" {refcons "" ns "" lowscore "" a100 "" b20 "" trf "" str "" repeat "" segdup "" selfchain "" cluster ""}
	} {
		puts stderr "===== $title ====="
		puts $o "\n ==================== $title ===================="
		set flist [comstats_filter $table $filters]
		foreach compar {df mm sm} {
			puts stderr $compar
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
						append oline \t[count $list [list dbsnp "" loc EXON effect $e sample $id]]
					}
				}
				puts $o $oline
			}
			set temp1 [count $list [list dbsnp "" loc EXON effect {FRAMESHIFT NONSENSE NONSTOP INSERT+ DELETE+ INSERT DELETE}]]
			set temp2 [count $list [list dbsnp dbsnp loc EXON effect {FRAMESHIFT NONSENSE NONSTOP INSERT+ DELETE+ INSERT DELETE}]]
			set oline "*BAD*\t[expr {$temp1+$temp2}]\t$temp1\t$temp2\t"
			if {$compar eq "df"} {
				foreach id $idlist {
					append oline \t[count $filters [list dbsnp "" loc EXON effect {FRAMESHIFT NONSENSE NONSTOP INSERT+ DELETE+ INSERT DELETE} sample $id]]
					append oline \t[count $filters [list dbsnp dbsnp loc EXON effect {FRAMESHIFT NONSENSE NONSTOP INSERT+ DELETE+ INSERT DELETE} sample $id]]
				}
			}
			puts $o $oline
		}
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
		type snp dbsnp dbsnp ns n? lowscore ls 
		refcons rc trf trf str str repeat rp segdup sd selfchain sc cluster cl
		a100 a100 a70 a70 b20 b20 b30 b30
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


set compar_file compar/78vs79_compar-filter-rna.tsv
set pvtfile 

lappend auto_path ~/dev/completegenomics/lib
package require Extral

cd /complgen/compar
set file 78vs79_mismatch.tsv
set ofile summary.tsv
catch {file delete summary.tsv}

set file compar/78vs79_compar_pvt.tsv
set file compar_GS102_GS103/pvtcompar_GS102_GS103.tsv

}
