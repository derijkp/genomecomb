#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

if 0 {
	set path /complgen/refseq/hg18
	set build hg18
	set dbname snp130
}

package require BioTcl

proc download_dbsnp_convline {line} {
	foreach {chr begin end class ref observed name alleleFreqs} $line break
	set rest [lrange $line 8 end]
	set observed [split $observed /]
	set alleleFreqs [split [string trimright $alleleFreqs ,] ,]
	set rp [lsearch $observed $ref]
	if {$rp == -1} {
		set complement [seq_complement $observed]
		set rp [lsearch $complement $ref]
		if {$rp != -1} {
			set observed $complement
		}
	}
	set alts {}
	set freqs {}
	set rlen [string length $ref]
	set alen {}
	if {[llength $alleleFreqs] > [llength $observed]} {
		set alleleFreqs [lrange $alleleFreqs 0 [expr {[llength $observed]-1}]]
	} elseif {[llength $alleleFreqs] < [llength $observed]} {
		lappend alleleFreqs {*}[list_fill [expr {[llength $observed]-[llength $alleleFreqs]}] 0]
	}
	foreach el $observed freq $alleleFreqs {
		if {$el eq $ref} continue
		if {$el eq "-"} {set el {}}
		lappend alts $el
		lappend freqs [format "%.3f" [expr {100.0*$freq}]]
		lappend alen [string length $el]
	}
	if {$rlen == 0} {
		set type ins
	} elseif {[inlist $alen 0]} {
		set type del
	} elseif {$rlen == 1 && [inlist $alen 1]} {
		set type snp
	} else {
		set type sub
	}
	set len [llength $alts]
	set names [list_fill $len $name]
	list $chr $begin $end $type $ref $alts $names $freqs {*}$rest
}

proc cg_download_dbsnp {resultfile build dbname} {
	set ufilename $resultfile.ucsc
	puts "Making $resultfile"
	if {[file exists $resultfile]} {
		puts "file $resultfile exists: skipping"
		return
	}
	catch {file mkdir [file dir $resultfile]}
	if {![file exists $ufilename]} {
		puts "Downloading $ufilename"
		cg_download_ucsc $ufilename $build $dbname
	}
	puts "Converting $ufilename"
	catch {close $f} ; catch {close $o}
	set f [open $ufilename]
	set o [open $resultfile.temp w]
	set header [split [gets $f] \t]
	set poss [list_cor $header {chrom start end class refNCBI observed name alleleFreqs avHet avHetSE strand molType valid func weight exceptions submitterCount submitters alleleFreqCount alleles alleleNs bitfields}]
	puts $o [join {chrom start end type ref alt name freqp avHet avHetSE strand molType valid func weight exceptions submitterCount submitters alleleFreqCount alleles alleleNs bitfields} \t]
	set pline [list_sub [split [gets $f] \t] $poss]
	set pline [download_dbsnp_convline $pline]
	set num 0 ; set next 100000
	while 1 {
		incr num; if {$num >= $next} {puts $num; incr next 100000}
		set line [split [gets $f] \t]
		set line [list_sub $line $poss]
		set line [download_dbsnp_convline $line]
		if {[lrange $pline 0 3] eq [lrange $line 0 3]} {
			foreach {alt name freq} [list_sub $pline {5 6 7}] break
			lappend alt {*}[lindex $line 5]
			lappend name {*}[lindex $line 6]
			lappend freq {*}[lindex $line 7]
			lset pline 5 $alt
			lset pline 6 $name
			lset pline 7 $freq
		} else {
			unset -nocomplain a
			foreach {alts names freqs} [list_sub $pline {5 6 7}] break
			foreach alt $alts name $names freq $freqs {
				lappend a($alt) [list $name $freq]
			}
			set alts {}
			set names {}
			set freqs {}
			foreach alt [lsort [array names a]] {
				lappend alts $alt
				lappend names [join [list_subindex $a($alt) 0] \;]
				set fr [list_remove [list_subindex $a($alt) 1] {}]
				if {[llength $fr] > 1} {set fr [lmath_max $fr]}
				lappend freqs $fr
			}
			lset pline 5 [join $alts ,]
			lset pline 6 [join $names ,]
			lset pline 7 [join $freqs ,]
			puts $o [join $pline \t]
			set pline $line
		}
		if {[eof $f]} break
	}
	close $f ; close $o
	file rename -force $resultfile.ucsc.info $resultfile.info
	puts "Sorting $resultfile"
	cg select -s - $resultfile.temp $resultfile.temp2
	file rename -force $resultfile.temp2 $resultfile
	file delete $resultfile.temp $resultfile.ucsc
}
