if 0 {
set multicomparfile /complgen/multicompar/compar.tsv
set mpvtfile /complgen/multicompar/mpvt_GS102_GS103_rtg102_rtg103.tsv
}

proc compare_mpvt {multicomparfile mpvtfile} {
	catch {close $f}
	set f [rzopen $multicomparfile]
	set samples [lrange [split [file root [file tail $mpvtfile]] _] 1 end]
	unset -nocomplain a
	set header [split [gets $f] \t]
	set fields {type trf str repeat segdup selfchain checked}
	set poss [list_cor $header $fields]
	set ufields {type xRef exonCategory aaCategory}
	set uposs [list_cor $header $ufields]
	unset -nocomplain possa
	foreach field {
		sequenced alleleSeq1 alleleSeq2 totalScore1 totalScore2
		refscore coverage refcons cluster posterior
	} {
		foreach sample $samples {
			set possa($field,$sample) [lsearch $header ${field}-$sample]
		}
	}
	foreach sample $samples {
		set possa(alleles,$sample) $possa(alleleSeq1,$sample)
		lappend possa(alleles,$sample) $possa(alleleSeq2,$sample)
		set possa(scores,$sample) $possa(totalScore1,$sample)
		lappend possa(scores,$sample) $possa(totalScore2,$sample)
		lappend possa(sequenced) $possa(sequenced,$sample)
	}
	# make header
	set oheader {}
	set pos 1
	foreach s1 [lrange $samples 0 end-1] {
		foreach s2 [lrange $samples $pos end] {
			lappend oheader compar_${s1}_${s2}
		}
		incr pos
	}
	lappend oheader {*}$fields
	lappend oheader dbsnp loc effect
	foreach sample $samples {
		foreach field {ns lowscore refcons cluster coverage sequenced} {
			lappend oheader ${field}-$sample
		}
	}
	lappend oheader tnum
#set o [open /complgen/multicompar/test.tsv w]
#puts $o "[join $header \t]"
#puts $o "chromosome\tbegin\tend\t[join $oheader \t]"
	set num 0
	while {![eof $f]} {
		incr num
		if {![expr {$num%100000}]} {putslog $num}
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		foreach $ufields [list_sub $line $uposs] break
		set result {}
		# do comparisons
		set sequenced [list_sub $line $possa(sequenced)]
		# if {![inlist $sequenced v]} continue
		set pos 1
#set keepcompar 0
		foreach s1 [lrange $samples 0 end-1] seq1 [lrange $sequenced 0 end-1] {
			foreach s2 [lrange $samples $pos end] seq2 [lrange $sequenced $pos end] {
				set alleles1 [list_sub $line $possa(alleles,$s1)]
				set alleles2 [list_sub $line $possa(alleles,$s2)]
				set compar [mcompar $seq1 $seq2 {*}$alleles1 {*}$alleles2]
				lappend result $compar
#if {$s1 eq "GS102" && $s2 eq "rtg102"} {set keepcompar $compar}
			}
			incr pos
		}
		# add data to result
		lappend result {*}[list_sub $line $poss]
		# dbsnp
		if {[regexp dbsnp $xRef]} {set dbsnp dbsnp} else {set dbsnp {}}
		lappend result $dbsnp
		# loc
		set loc ""
		foreach temp {EXON UTR BEGIN END INTRON} {
			if {[regexp $temp $exonCategory]} {
				set loc $temp
				break
			}
		}
		lappend result $loc
		# effect
		set effect OTHER
		foreach temp {FRAMESHIFT NONSENSE NONSTOP INSERT+ DELETE+ INSERT DELETE MISSENSE COMPATIBLE UNDEFINED} {
			if {[regexp $temp $aaCategory]} {
				set effect $temp
				break
			}
		}
		lappend result $effect
		# per sample
		foreach sample $samples {
			# ns
			set alleles [list_sub $line $possa(alleles,$sample)]
			if {[regexp {[N?]} $alleles]} {set ns ns} else {set ns ""}
			# lowscore
			if {[llength $possa(scores,$sample)] == 2} {
				set scores [list_sub $line $possa(scores,$sample)]
				set minscore [lindex [lsort -dict $scores] 0]
				if {$minscore < 60} {set lowscore ls} else {set lowscore ""}
			} else {
				# lowpost
				set minpost [lindex $line $possa(posterior,$sample)]
				if {[isdouble $minpost]} {
					if {$minpost < 4} {set lowscore ls} else {set lowscore ""}
				} else {
					set lowscore ""
				}
			}
			set refcons [lindex $line $possa(refcons,$sample)]
			set cluster [lindex $line $possa(cluster,$sample)]
			set coverage [lindex $line $possa(coverage,$sample)]
			if {![isint $coverage]} {
				set coverage vl
			} elseif {$coverage < 10} {
				set coverage vl
			} elseif {$coverage < 20} {
				set coverage l
			} elseif {$coverage > 100} {
				set coverage h
			} else {
				set coverage n
			}
			set sequenced [lindex $line $possa(sequenced,$sample)]
			lappend result $ns $lowscore $refcons $cluster $coverage $sequenced
		}
#if {$keepcompar eq "sm" || $keepcompar eq "i"} {
#	puts $o [join $line \t]
#}
		if {![info exists a($result)]} {
			set a($result) 1
		} else {
			incr a($result)
		}
	}
#close $o
	set o [open $mpvtfile w]
	puts $o [join $oheader \t]
	set list [lsort -dict [array names a]]
	foreach line $list {
		puts $o [join $line \t]\t$a($line)
	}
}

proc mpvt_select {query} {
	global mpvt
	set mpvtfile $mpvt(file)
	set num [cg select -f tnum -q $query $mpvtfile | nsum]
}

proc mpvt_totals {name query} {
	global mpvt
	set mpvtfile $mpvt(file)
	set changelist $mpvt(changelist)
	set name [string_change $name $changelist]
	set query [string_change $query $changelist]
	set totals {}
	foreach compar {sm df mm un} totlist $totals {
		set uquery [string trim [string_change $query [list COMPAR "\"$compar\"" \n {}]]]
		set temp {}
		lappend temp ${compar}_$name [cg select -f tnum -q $uquery $mpvtfile | nsum]
		lappend totals $temp
	}
	return $totals
}

proc mpvt_comparinfo {{label {}} query totals} {
	global mpvt
	set mpvtfile $mpvt(file)
	set o $mpvt(o)
	set changelist $mpvt(changelist)
	set label [string_change $label $changelist]
	set query [string_change $query $changelist]
	puts $o "---------- $label ----------"
	puts $o "           [string trim $query]"
	foreach compar {sm df mm un} totlist $totals {
		set uquery [string trim [string_change $query [list COMPAR "\"$compar\"" \n {}]]]
		set num [cg select -f tnum -q $uquery $mpvtfile | nsum]
		set result [list $compar $num]
		foreach {vs tot} $totlist {
			if {![isint $tot] || ($tot == 0)} {
				lappend result ""
			} else {
				lappend result "[format %.2f [expr {100.0*$num/$tot}]]% ($vs)"
			}
		}
		puts $o [join $result ]
	}
}

proc mpvt_comparinfoperc {num tot} {
	if {![isdouble $num]} {return {}}
	if {$tot == 0} {return ""}
	format %.2f [expr {100.0*$num/$tot}]
}

proc mpvt_comparinfo {{label {}} query totals} {
	global mpvt
	set mpvtfile $mpvt(file)
	set o $mpvt(o)
	set changelist $mpvt(changelist)
	set label [string_change $label $changelist]
	set query [string_change $query $changelist]
	putslog $label
	unset -nocomplain a
	set result {}
	lappend result [string trim $label]
	set totals [lrange $totals 0 2]
	foreach compar {sm df mm} totlist $totals {
		set uquery [string trim [string_change $query [list COMPAR "\"$compar\"" \n {}]]]
		set num [cg select -f tnum -q $uquery $mpvtfile | nsum]
		set a($compar) $num
		set tot {}
		foreach {vs tot} $totlist {
			set a($compar,tot) $tot
		}
		lappend result $num [mpvt_comparinfoperc $num $tot]%
	}
	set diffs [expr {$a(df)+$a(mm)}]
	set pdiffs [mpvt_comparinfoperc $diffs [expr {$a(df,tot)+$a(mm,tot)}]]
	lappend result $diffs [format %.2f $pdiffs]
	if {$a(sm,tot) != 0} {
		set psm [expr {100.0*$a(sm)/$a(sm,tot)}]
		if {$psm == 0} {
			lappend result #
		} else {
			lappend result [format %.2f [expr {$pdiffs/$psm}]]
		}
	} else {
		lappend result #
	}
	regsub -all \n [string trim $query] { } query
	regsub -all \t [string trim $query] { } query
	lappend result $query
	puts $o [join $result \t]
}

proc mpvt_summary {mpvtfile mpvtsummaryfile} {
	global filters mpvt

	set samples [lrange [split [file root [file tail $mpvtsummaryfile]] _] 1 end]
	if {[llength $samples] != 4} {
		error "file $mpvtsummaryfile does not contain four samples (in the form prefix_sample1_sample2_rtgsample1_rtgsample2)"
	}
	foreach {sample1 sample2 rtgsample1 rtgsample2} $samples break
	set mpvt(file) $mpvtfile
	catch {close $o}
	set o [open $mpvtsummaryfile w]
	set mpvt(o) $o
	set changelist [list @SAMPLE1@ $sample1 @SAMPLE2@ $sample2 @RTGSAMPLE1@ $rtgsample1 @RTGSAMPLE2@ $rtgsample2]
	set mpvt(changelist) $changelist
	# set o stdout
	set filterstemplate {
		refcons {$refcons-@SAMPLE1@ == "" && $refcons-@SAMPLE2@ == ""}
		lowscore {$lowscore-@SAMPLE1@ == "" && $lowscore-@SAMPLE2@ == ""}
		cluster {$cluster-@SAMPLE1@ == "" && $cluster-@SAMPLE2@ == ""}
		trf {$trf == ""}
		str {$str == ""}
		coverage {$coverage-@SAMPLE1@ == "n" && $coverage-@SAMPLE2@ == "n"}
		coverage10 {($coverage-@SAMPLE1@ == "l" || $coverage-@SAMPLE1@ == "n") && ($coverage-@SAMPLE2@ == "l" || $coverage-@SAMPLE2@ == "n")}
		segdup {$segdup == ""}
		selfchain {$selfchain == ""}
		repeat {$repeat == ""}
		rtg {$compar_@SAMPLE1@_@RTGSAMPLE1@ ~ /sm|i/ && $compar_@SAMPLE2@_@RTGSAMPLE2@ ~ /sm|i/}
	}
	array set filters [string_change $filterstemplate $changelist]
	set oheader [join {
		label
		"sm" "% sm" "df" "% df" "mm" "% mm" "df+mm" "% df+mm"
		"Ratio % df+mm / % sm" query
	} \t]
	
	set totalsfull [mpvt_totals @SAMPLE1@_@SAMPLE2@ {$compar_@SAMPLE1@_@SAMPLE2@ == COMPAR}]
	set totalssnp [mpvt_totals @SAMPLE1@_@SAMPLE2@ {$compar_@SAMPLE1@_@SAMPLE2@ == COMPAR && $type == "snp"}]
	set totalsbasic [mpvt_totals basic_@SAMPLE1@_@SAMPLE2@ {$compar_@SAMPLE1@_@SAMPLE2@ == COMPAR && $refcons-@SAMPLE1@ == "" && $refcons-@SAMPLE2@ == "" && $type == "snp"}]
	set rtgtotals [mpvt_totals @RTGSAMPLE1@_@RTGSAMPLE2@ {$compar_@RTGSAMPLE1@_@RTGSAMPLE2@ == COMPAR}]
	
	
	puts $o "\n-------------------- Totals --------------------\n"
	puts $o $oheader
	mpvt_comparinfo "total" {
		$compar_@SAMPLE1@_@SAMPLE2@ == COMPAR
	} $totalsfull
	
	mpvt_comparinfo "only snps" {
		$compar_@SAMPLE1@_@SAMPLE2@ == COMPAR && $type == "snp"
	} $totalsfull
	
	flush $o
	
	puts $o "\n-------------------- indiv filter quality: numbers filtered, snps only, % versus totalsnps (per type) --------------------\n"
	puts $o $oheader
	
	mpvt_comparinfo "indels" {
		$compar_@SAMPLE1@_@SAMPLE2@ == COMPAR
		&& $type != "snp"
	} $totalsfull
	
	foreach {fname} {
		refcons lowscore cluster trf str
		coverage coverage10 segdup selfchain repeat rtg
	} {
		set query [subst {
			\$compar_@SAMPLE1@_@SAMPLE2@ == COMPAR && \$type == "snp"
			&& !( $filters($fname) )
		}]
		mpvt_comparinfo $fname $query $totalssnp
	}
	
	flush $o
	
	puts $o "\n-------------------- ${sample1}_${sample2} filters (% vs basic filter) --------------------\n"
	puts $o $oheader
	
	set query {$compar_@SAMPLE1@_@SAMPLE2@ == COMPAR && $refcons-@SAMPLE1@ == "" && $refcons-@SAMPLE2@ == "" && $type == "snp"}
	mpvt_comparinfo "total snps - refcons" $query $totalsbasic
	foreach {fname} {
		lowscore cluster trf str
		coverage segdup selfchain repeat
	} {
		set name " - $fname"
		append query " && $filters($fname)"
		mpvt_comparinfo $name $query $totalsbasic
	}
	
	flush $o
	
	puts $o "\n-------------------- ${sample1}_${sample2} filters - not same geno in rtg--------------------\n"
	puts $o $oheader
	
	set query {$compar_@SAMPLE1@_@SAMPLE2@ == COMPAR && $refcons-@SAMPLE1@ == "" && $refcons-@SAMPLE2@ == "" && $type == "snp"}
	mpvt_comparinfo "total - refcons" $query $totalsbasic
	foreach {fname} {
		rtg lowscore cluster trf str
		coverage segdup selfchain repeat
	} {
		set name " - $fname"
		append query " && $filters($fname)"
		mpvt_comparinfo $name $query $totalsbasic
	}
	
	puts $o "\n-------------------- ${sample1}_${sample2} filters - other order --------------------\n"
	puts $o $oheader
	
	set query {$compar_@SAMPLE1@_@SAMPLE2@ == COMPAR && $refcons-@SAMPLE1@ == "" && $refcons-@SAMPLE2@ == "" && $type == "snp"}
	mpvt_comparinfo "total - refcons" $query $totalsbasic
	foreach {fname} {
		str trf cluster rtg coverage lowscore 
		segdup selfchain repeat
	} {
		set name " - $fname"
		append query " && $filters($fname)"
		mpvt_comparinfo $name $query $totalsbasic
	}
	
	puts $o "\n-------------------- ${sample1}_${sample2} filters - paper order coverage10 --------------------\n"
	puts $o $oheader
	
	set query {$compar_@SAMPLE1@_@SAMPLE2@ == COMPAR && $refcons-@SAMPLE1@ == "" && $refcons-@SAMPLE2@ == "" && $type == "snp"}
	mpvt_comparinfo "total - refcons" $query $totalsbasic
	foreach {fname} {
		coverage10 lowscore cluster str trf segdup rtg 
	} {
		set name " - $fname"
		append query " && $filters($fname)"
		mpvt_comparinfo $name $query $totalsbasic
	}
	
	puts $o "\n-------------------- ${sample1}_${sample2} filters - paper order coverage20 --------------------\n"
	puts $o $oheader
	
	set query {$compar_@SAMPLE1@_@SAMPLE2@ == COMPAR && $refcons-@SAMPLE1@ == "" && $refcons-@SAMPLE2@ == "" && $type == "snp"}
	mpvt_comparinfo "total - refcons" $query $totalsbasic
	foreach {fname} {
		coverage lowscore cluster str trf segdup rtg 
	} {
		set name " - $fname"
		append query " && $filters($fname)"
		mpvt_comparinfo $name $query $totalsbasic
	}
	
	puts $o "\n-------------------- ${sample1}_${sample2} extra filters --------------------\n"
	puts $o $oheader
	
	mpvt_comparinfo "total - refcons - not same geno in rtg - trf" {
		$compar_@SAMPLE1@_@SAMPLE2@ == COMPAR && $refcons-@SAMPLE1@ == "" && $refcons-@SAMPLE2@ == "" && $type == "snp"
		&& $compar_@SAMPLE1@_@RTGSAMPLE1@ ~ /sm|i/ && $compar_@SAMPLE2@_@RTGSAMPLE2@ ~ /sm|i/
		&& $trf == ""
	} $totalsbasic
	
	mpvt_comparinfo "total - refcons - not same geno in rtg - trf - str" {
		$compar_@SAMPLE1@_@SAMPLE2@ == COMPAR && $refcons-@SAMPLE1@ == "" && $refcons-@SAMPLE2@ == "" && $type == "snp"
		&& $compar_@SAMPLE1@_@RTGSAMPLE1@ ~ /sm|i/ && $compar_@SAMPLE2@_@RTGSAMPLE2@ ~ /sm|i/
		&& $trf == "" && $str == ""
	} $totalsbasic
	
	mpvt_comparinfo "total - refcons - not same geno in rtg - lowscore - cluster - trf - str - coverage - segdup - notchecked" {
		$compar_@SAMPLE1@_@SAMPLE2@ == COMPAR && $refcons-@SAMPLE1@ == "" && $refcons-@SAMPLE2@ == "" && $type == "snp"
		&& $compar_@SAMPLE1@_@RTGSAMPLE1@ ~ /sm|i/ && $compar_@SAMPLE2@_@RTGSAMPLE2@ ~ /sm|i/
		&& $lowscore-@SAMPLE1@ == "" && $lowscore-@SAMPLE2@ == "" && $cluster-@SAMPLE1@ == "" && $cluster-@SAMPLE2@ == ""
		&& $trf == "" && $str == "" && $coverage-@SAMPLE1@ == "n" && $coverage-@SAMPLE2@ == "n"
		&& $segdup == "" && $checked == ""
	} $totalsbasic
	
	flush $o
	
	puts $o "\n-------------------- ${rtgsample1}_$rtgsample2 --------------------\n"
	puts $o $oheader
	
	mpvt_comparinfo "total (@RTGSAMPLE1@_@RTGSAMPLE2@)" {$compar_@RTGSAMPLE1@_@RTGSAMPLE2@ == COMPAR} $rtgtotals
	
	mpvt_comparinfo "total (@RTGSAMPLE1@_@RTGSAMPLE2@) - !lowcoverage" {
		$compar_@RTGSAMPLE1@_@RTGSAMPLE2@ == COMPAR && $coverage-@RTGSAMPLE1@ != "l" && $coverage-@RTGSAMPLE2@ != "l"
	} $rtgtotals
	
	mpvt_comparinfo "total (@RTGSAMPLE1@_@RTGSAMPLE2@) - lowscore" {
		$compar_@RTGSAMPLE1@_@RTGSAMPLE2@ == COMPAR && $lowscore-@RTGSAMPLE1@ == "" && $lowscore-@RTGSAMPLE2@ == ""
	} $rtgtotals


	
	set changelist [list @SAMPLE1@ $rtgsample1 @SAMPLE2@ $rtgsample2 @RTGSAMPLE1@ $rtgsample1 @RTGSAMPLE2@ $rtgsample2]
	set mpvt(changelist) $changelist
	array set filters [string_change $filterstemplate $changelist]
	set query {$compar_@SAMPLE1@_@SAMPLE2@ == COMPAR && $type == "snp"}
	set totalsrtg [mpvt_totals @SAMPLE1@_@SAMPLE2@ {$compar_@SAMPLE1@_@SAMPLE2@ == COMPAR && $type == "snp"}]
	mpvt_comparinfo "total" $query $totalsbasic
	foreach {fname} {
		str trf cluster coverage lowscore 
		segdup selfchain repeat
	} {
		set name " - $fname"
		append query " && $filters($fname)"
		mpvt_comparinfo $name $query $totalsrtg
	}


	flush $o
	
	close $o
}
