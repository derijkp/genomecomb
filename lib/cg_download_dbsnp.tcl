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
	set ufilename $resultfile.ucsc.zst
	puts "Making $resultfile"
	if {[file exists $resultfile]} {
		puts "file $resultfile exists: skipping"
		return
	}
	catch {file mkdir [file dir $resultfile]}
	if {![file exists $ufilename]} {
		puts "Downloading $ufilename"
		cg_download_ucsc -cl 1 $ufilename $build $dbname
	}
	puts "Converting $ufilename"
	catch {close $f} ; catch {close $o}
	set f [gzopen $ufilename]
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
	snp/snp/file rename -force -- $resultfile.ucsc.info [gzroot $resultfile].info
	puts "Sorting $resultfile"
	cg select -s - $resultfile.temp $resultfile.temp2[file extension $resultfile]
	file rename -force -- $resultfile.temp2[file extension $resultfile] $resultfile
	file delete $resultfile.temp $resultfile.ucsc
}

proc cg_download_dbsnp_new {resultfile build dbname} {
	regsub ^snp $dbname {} dbsnpversion
	file_write [gzroot $resultfile].info [subst [deindent {
		= dbsnp =
		
		== Download info ==
		dbname	dbsnp
		version	$dbsnpversion
		source	https://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp$dbsnpversion.bb
		time	[timestamp]
		
		== Description ==
		
		This track shows short genetic variants (up to approximately 50 base pairs) 
		from dbSNP build $dbsnpversion: single-nucleotide variants (SNVs), small insertions, deletions, 
		and complex deletion/insertions (indels), relative to the reference genome assembly. 
		Most variants in dbSNP are rare, not true polymorphisms, and some variants are known 
		to be pathogenic. 
		
		== Category ==
		Annotation
	}]]
	wgetfile https://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp$dbsnpversion.bb
	wgetfile http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bigBedToBed
	chmod u+x bigBedToBed
	exec ./bigBedToBed dbSnp$dbsnpversion.bb dbSnp$dbsnpversion.bb.bed
	file delete bigBedToBed dbSnp$dbsnpversion.bb
	#
	catch {gzclose $f} ; catch {gzclose $o}
	set f [open dbSnp$dbsnpversion.bb.bed]
	set o [wgzopen $resultfile.temp.zst]
	puts $o [join {chromosome begin end type ref alt name} \t]
	array set typea {snv snp mnv sub}
	set num 0 ; set fnum 0
	while {[gets $f line] != -1} {
		if {[incr num] >= 1000000} {
			incr fnum 1
			set num 0
			puts $fnum
		}
		foreach {chromosome begin end name ref numalt alt shiftBases freqSourceCount minorAlleleFreq majorAllele minorAllele maxFuncImpact type ucscNotes} [split $line \t] break
		set alt [string trimright $alt ,]
		if {$type eq "identity"} continue
		if {$type eq "snv"} {
			set type snp
		} elseif {$type eq "mnv"} {
			set type sub
		} elseif {$type eq "del"} {
			set len [string length $ref]
			if {$len > 1} {set ref $len}
		}
#		if {$numalt > 1} {
#			set name [join [list_fill $numalt $name] ,]
#		}
		puts $o [join [list $chromosome $begin $end $type $ref $alt $name] \t]
	}
	gzclose $f ; gzclose $o
	file rename $resultfile.temp.zst $resultfile
	# file rename -force -- $resultfile.ucsc.info [gzroot $resultfile].info
	puts "Sorting $resultfile"
	cg select -s - $resultfile.temp.zst $resultfile.temp2[file extension $resultfile]
	file rename -force -- $resultfile.temp2[file extension $resultfile] $resultfile
	file delete $resultfile.temp.zst dbSnp$dbsnpversion.bb.bed
}
