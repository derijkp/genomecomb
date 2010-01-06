proc pvtcomstats {} {
	set f stdin
	# set f [open $compar_file]
	unset -nocomplain a
	set header [split [gets $f] \t]
	set fields {compar sample chromosome type dbsnp ns lowscore loc effect refcons trf str rp sd sc}
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

proc count {var slist {cfilter {}}} {
	upvar $var a
	if {[llength $cfilter]} {
		set slist [comstats_filter $slist $cfilter]
	}
	set result 0
	foreach line $slist {
		incr result [get a($line) 0]
	}
	return $result
}

proc comstats {} {
	global header
	set f stdin
	# set f [open $file]
	unset -nocomplain a
	set header [split [gets $f] \t]
	set num 0
	while {![eof $f]} {
		incr num
		if {![expr {$num%100000}]} {puts stderr $num}
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set num [list_pop line]
		set a($line) $num
	}
	set o stdout
	array set trans {
		COMPATIBLE 1 MISSENSE 2 DELETE 3 INSERT 4 DELETE+ 5 INSERT+ 6 NONSTOP 7 NONSENSE 8 FRAMESHIFT 9
	}
	set pos [lsearch $header sample]
	# set idlist [list_remdup [list_subindex [array names a] $pos]]
	set idlist [list_remdup [split [list_remdup [list_subindex [array names a] $pos]] " ,"]]
	set filterlist {0 1 2 3}
	set loclist {EXON UTR BEGIN END INTRON OTHER ""}
	set effectlist {FRAMESHIFT NONSENSE NONSTOP INSERT+ DELETE+ INSERT DELETE MISSENSE COMPATIBLE UNDEFINED OTHER}
	# set o [open $ofile a]
	set o stdout
	set alllist [array names a]
	foreach {title filters} {
		nofilter {}
		"filter 1 (removed ref-(in)consistent)" {refcons ""}
		"filter 3 (removed ref-(in)consistent, ?/N)" {refcons "" ns ""}
		"filter 4 (removed ref-(in)consistent, ?/N,  scores < 60)" {refcons "" ns "" lowscore ""}
		"filter 5 (removed ref-(in)consistent, ?/N, scores < 60, tandem repeats)" {refcons "" ns "" lowscore "" trf ""}
		"filter 6 (removed ref-(in)consistent, ?/N, scores < 60, tandem repeats, str)" {refcons "" ns "" lowscore "" trf "" str ""}
		"filter 7 (removed ref-(in)consistent, ?/N, scores < 60, tandem repeats, str, repeats)" {refcons "" ns "" lowscore "" trf "" str "" rp ""}
		"filter 8 (removed ref-(in)consistent, ?/N, scores < 60, tandem repeats, str, repeats, segdup)" {refcons "" ns "" lowscore "" trf "" str "" rp "" sd ""}
		"filter 9 (removed ref-(in)consistent, ?/N, scores < 60, tandem repeats, str, repeats, segdup, self chained)" {refcons "" ns "" lowscore "" trf "" str "" rp "" sd "" sc ""}
	} {
		puts stderr "===== $title ====="
		puts $o "\n ==================== $title ===================="
		set flist [comstats_filter $alllist $filters]
		foreach compar {df mm sm} {
			puts stderr $compar
			puts $o "\n---------- $compar: $title ----------"
			set list [comstats_filter $flist [list compar $compar]]
			set checkall 0
			foreach name $list {incr checkall $a($name)}
			puts $checkall
			set all [count a $list]
			puts $o "total:\t$all"
			set snps [count a $list {type snp}]
			puts $o "snps:\t$snps ([expr {$all-$snps}] other)"
			if {$compar eq "df"} {
				set idlistheader \t\t[join $idlist \t]
			} else {
				set idlistheader ""
			}
			puts $o "\nloc\ttotal\tnew\tdbsnp$idlistheader"
			foreach l $loclist {
				set temp1 [count a $list [list dbsnp "" loc $l]]
				set temp2 [count a $list [list dbsnp dbsnp loc $l]]
				set oline "$l:\t[expr {$temp1+$temp2}]\t${temp1}\t${temp2}\t"
				if {$compar eq "df"} {
					foreach id $idlist {
						append oline \t[count a $list [list loc $l sample $id]]
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
				set temp1 [count a $list [list dbsnp "" loc EXON effect $e]]
				set temp2 [count a $list [list dbsnp dbsnp loc EXON effect $e]]
				set oline "$e\t[expr {$temp1+$temp2}]\t$temp1\t$temp2\t"
				if {$compar eq "df"} {
					foreach id $idlist {
						append oline \t[count a $list [list dbsnp "" loc EXON effect $e sample $id]]
						append oline \t[count a $list [list dbsnp "" loc EXON effect $e sample $id]]
					}
				}
				puts $o $oline
			}
			set temp1 [count a $list [list dbsnp "" loc EXON effect {FRAMESHIFT NONSENSE NONSTOP INSERT+ DELETE+ INSERT DELETE}]]
			set temp2 [count a $list [list dbsnp dbsnp loc EXON effect {FRAMESHIFT NONSENSE NONSTOP INSERT+ DELETE+ INSERT DELETE}]]
			set oline "*BAD*\t[expr {$temp1+$temp2}]\t$temp1\t$temp2\t"
			if {$compar eq "df"} {
				foreach id $idlist {
					append oline \t[count a $filters [list dbsnp "" loc EXON effect {FRAMESHIFT NONSENSE NONSTOP INSERT+ DELETE+ INSERT DELETE} sample $id]]
					append oline \t[count a $filters [list dbsnp dbsnp loc EXON effect {FRAMESHIFT NONSENSE NONSTOP INSERT+ DELETE+ INSERT DELETE} sample $id]]
				}
			}
			puts $o $oline
		}
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
comstats summary.tsv diff.tsv
comstats summary.tsv same.tsv
comstats summary.tsv uniqueseq.tsv

set file compar/78vs79_compar_pvt.tsv

}
