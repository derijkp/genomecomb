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

proc scoredist file {
	global compara scorea refa cova mcova

	catch {close $f}
	set f [rzopen $file]
	set header [split [gets $f] \t]
	set scoreposs [list_cor $header {totalScore1-1 totalScore2-1 totalScore1-2 totalScore2-2}]
	set refscoreposs [list_cor $header {refscore-1 refscore-2}]
	set coverageposs [list_cor $header {coverage-1 coverage-2}]
	set comparpos [lsearch $header compar]
	set dbsnppos [lsearch $header xRef]
	unset -nocomplain scorea
	unset -nocomplain refa
	unset -nocomplain cova
	unset -nocomplain mcova
	set num 0
	while {![eof $f]} {
		incr num
		if {![expr {$num%10000}]} {puts $num}
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set score [lmath_min [list_remove [list_sub $line $scoreposs] {} -]]
		set compar [lindex $line $comparpos]
		set refscore [lmath_max [list_remove [list_sub $line $refscoreposs] {} -]]
		set coverage [lmath_min [list_remove [list_sub $line $coverageposs] {} -]]
		set mcoverage [lmath_max [list_remove [list_sub $line $coverageposs] {} -]]
		scoredist_incrset scorea($score) $compar
		scoredist_incrset refa($refscore) $compar
		scoredist_incrset cova($coverage) $compar
		scoredist_incrset mcova($mcoverage) $compar
		if {![regexp dbsnp [lindex $line $dbsnppos]]} {
			scoredist_incrset scorea($score) ${compar}-new
			scoredist_incrset refa($refscore) ${compar}-new
			scoredist_incrset cova($coverage) ${compar}-new
			scoredist_incrset mcova($mcoverage) ${compar}-new
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

	foreach {head name var} {
		score scoredist scorea
		refscore refscoredist refa
		coverage coveragedist cova
		maxcoverage maxcoveragedist mcova
	} {
		upvar #0 $var a
		set ofile [file join [file dir $file] $name-[file tail [file root [rzroot $file]]].tsv]
		set o [open $ofile w]
		puts $o [join [list $head sm df mm df+mm un {} sm-new df-new mm-new df+mm-new un-new] \t]
		set scores [lsort -integer [list_remove [array names a] {}]]
		foreach score $scores {
			set line $a($score)
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