proc scoredist_incrset {var compar} {
	global compara
	upvar $var line
	set pos $compara($compar)
	if {![info exists line]} {
		set line {0 0 0 0 0 0 0 0}
	}
	set v [lindex $line $pos]
	incr v
	lset line $pos $v
}

array set compara {sm 0 df 1 mm 2 un 3 sm-new 4 df-new 5 mm-new 6 un-new 7}

proc mcompar {sequenced1 sequenced2 allele11 allele21 allele12 allele22} {
	if {$sequenced1 eq "u" || $sequenced2 eq "u"} {return un}
	if {$sequenced1 ne "v" && $sequenced2 ne "v"} {return i}
	foreach {s1 s2} [lsort [list $sequenced1 $sequenced2]] break
	if {$s1 eq "v"} {
		set alleles1 [lsort [list $allele11 $allele21]]
		set alleles2 [lsort [list $allele12 $allele22]]
		if {$alleles1 eq $alleles2} {
			return sm
		} else {
			return mm
		}
	} else {
		return df
	}
}

proc scoredist {file sample1 sample2} {
	global compara scorea refa cova mcova

	catch {close $f}
	set f [rzopen $file]
	set header [split [gets $f] \t]
#	set comparpos [lsearch $header compar]
	set possa(score) [list_cor $header [list totalScore1-$sample1 totalScore2-$sample1 totalScore1-$sample2 totalScore2-$sample2]]
	foreach field {refscore coverage posterior alleleSeq1 alleleSeq2 sequenced refcons} {
		set possa($field) [list_cor $header [list $field-$sample1 $field-$sample2]]
	}
	set alleleposs [list_merge $possa(alleleSeq1) $possa(alleleSeq2)]
	set alleleposs1 [lrange $alleleposs 0 1]
	set alleleposs2 [lrange $alleleposs 2 3]
	set alleleposs2 [lrange $alleleposs 2 3]
	set dbsnppos [lsearch $header xRef]
	unset -nocomplain sa
	set num 0
	while {![eof $f]} {
		incr num
		if {![expr {$num%10000}]} {puts $num}
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set sequenced [list_sub $line $possa(sequenced)]
		if {[lsearch $sequenced v] == -1} continue
		set compar [mcompar {*}$sequenced {*}[list_sub $line $alleleposs1] {*}[list_sub $line $alleleposs2]]
		if {![regexp dbsnp [lindex $line $dbsnppos]]} {set new 1} else {set new 0}
		foreach {field take} {score min refscore max coverage min posterior min} {
			set list [list_sub $line $possa($field)]
			if {![llength [list_remove $list {}]]} continue
			if {$take eq "max"} {
				set score [max {*}$list]
			} else {
				set score [min {*}$list]
			}
			scoredist_incrset sa($field,$score) $compar
			if {$new} {
				scoredist_incrset sa($field,$score) ${compar}-new
			}
			set pos [lsearch $list $score]
			if {$field eq "score"} {
				if {$pos < 2} {set sample $sample1} else {set sample $sample2}
			} else {
				if {$pos < 1} {set sample $sample1} else {set sample $sample2}
			}
			scoredist_incrset sa($field-$sample,$score) $compar
			if {$new} {
				scoredist_incrset sa($field-$sample,$score) ${compar}-new
			}
			set rc [list_sub $line $possa(refcons)]
			if {![inlist $rc rc]} {
				scoredist_incrset sa(${field}-rc,$score) $compar
				if {$new} {
					scoredist_incrset sa(${field}-rc,$score) ${compar}-new
				}
			}
		}
	}
	close $f

#	# totals
#	set totals {0 0 0 0 0 0 0 0}
#	foreach score [array names scorea] {
#		set totals [lmath_calc $totals + $scorea($score)]
#	}
#	set temp {}
#	foreach el $totals {
#		lappend temp [expr {double($el)}]
#	}
#	set totals $temp
#	set totals [linsert $totals 3 [expr {[lindex $totals 1]+[lindex $totals 2]}]]
#	set totals [linsert $totals 5 {}]
#	set totals [linsert $totals 9 [expr {[lindex $totals 7]+[lindex $totals 8]}]]
	foreach {field name} [list \
		score scoredist \
		refscore refscoredist \
		coverage coveragedist \
		posterior posteriordist \
		score-rc scoredist-rc \
		refscore-rc refscoredist-rc \
		coverage-rc coveragedist-rc \
		posterior-rc posteriordist-rc \
		score-$sample1 scoredist-$sample1 \
		refscore-$sample1 refscoredist-$sample1 \
		coverage-$sample1 coveragedist-$sample1 \
		posterior-$sample1 posteriordist-$sample1 \
		score-$sample2 scoredist-$sample2 \
		refscore-$sample2 refscoredist-$sample2 \
		coverage-$sample2 coveragedist-$sample2 \
		posterior-$sample2 posteriordist-$sample2 \
	] {
		set base [file tail [file root [rzroot $file]]]
		set ofile [file join [file dir $file] $name-$base-$sample1-$sample2.tsv]
		set names [array names sa $field,*]
		regsub -all $field, $names {} scores
		set scores [lsort -real [list_remove $scores {} - ?]]
		if {![llength $scores]} continue
		set o [open $ofile w]
		puts $o [join [list $field sm df mm df+mm un {} sm-new df-new mm-new df+mm-new un-new] \t]
		foreach score $scores {
			set line $sa($field,$score)
			if {![llength $line]} continue
			set line [linsert $line 3 [expr {[lindex $line 1]+[lindex $line 2]}]]
			set line [linsert $line 5 {}]
			set line [linsert $line 9 [expr {[lindex $line 7]+[lindex $line 8]}]]
#			lappend line {}
#			foreach v [lmath_calc [lrange $line 0 4] / [lrange $totals 0 4]] {
#				lappend line [format %.3f $v]
#			}
#			lappend line {}
#			foreach v [lmath_calc [lrange $line 6 10] / [lrange $totals 6 10]] {
#				lappend line [format %.3f $v]
#			}
			puts $o $score\t[join $line \t]
		}
		close $o
		puts "finished $ofile"
	}

}

if 0 {

	lappend auto_path ~/dev/completegenomics/lib
	package require Extral
	package require Tclx
	signal -restart error SIGINT
	set file /complgen/multicompar/compar.tsv
	set sample1 GS102
	set sample2 GS103

	set sample1 rtg102
	set sample2 rtg103
	
	set files [glob /complgen/compar_*/fcompar_*.tsv.gz]
	foreach file $files {
		puts $file
		scoredist $file
	}
	

	set file /complgen/compar_GS100_GS102/fcompar_GS100_GS102.tsv.gz
	set file /complgen/compar_GS100_GS101A01/fcompar_GS100_GS101A01.tsv.gz
	set file /complgen/compar_GS102_GS103/fcompar_GS102_GS103.tsv.gz

	cg select -q '$compar == "sm" && min($totalScore1-1,$totalScore1-1,$totalScore1-2,$totalScore2-2) == 20' < compar_GS102_GS103/fcompar_GS102_GS103.tsv | wc


}