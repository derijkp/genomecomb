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

proc downloaddb_dbsnp_convline {line} {
	foreach {chr begin end class ref observed name freq} $line break
	set rest [lrange $line 8 end]
	set observed [split $observed /]
	set rp [lsearch $observed $ref]
	if {$rp == -1} {
		set complement [seq_complement $observed]
		set rp [lsearch $complement $ref]
		if {$rp != -1} {
			set observed $complement
		}
	}
	set observed [list_remove $observed $ref]
	if {$ref eq "-"} {set ref ""}
	set observed [list_change $observed {- {}}]
	set rlen [string length $ref]
	set alen {}
	foreach el $observed {
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
	set len [llength $observed]
	set name [list_fill $len $name]
	set freq [list_fill $len $freq]
	list $chr $begin $end $type $ref $observed $name $freq {*}$rest
}

proc downloaddb_dbsnp {path build dbname} {
	set ufilename $path/tmp/$build/ucsc_${build}_$dbname.tsv
	set filename $path/$build/var_${build}_$dbname.tsv
	puts "Making $filename"
	if {[file exists $filename]} {
		puts "file $filename exists: skipping"
		return
	}
	catch {file mkdir [file dir $filename]}
	if {![file exists $ufilename]} {
		downloaddb $path/tmp $build $dbname
	}
	puts "Converting $ufilename"
	catch {close $f} ; catch {close $o}
	set f [open $ufilename]
	set o [open $filename.temp w]
	set header [split [gets $f] \t]
	set poss [list_cor $header {chrom start end class refNCBI observed name avHet avHetSE strand molType valid func weight exceptions submitterCount submitters alleleFreqCount alleles alleleNs alleleFreqs bitfields}]
	puts $o [join {chrom start end type ref alt name freq avHetSE strand molType valid func weight exceptions submitterCount submitters alleleFreqCount alleles alleleNs alleleFreqs bitfields} \t]
	set pline [list_sub [split [gets $f] \t] $poss]
	set pline [downloaddb_dbsnp_convline $pline]
	set num 0 ; set next 100000
	while 1 {
		incr num; if {$num >= $next} {puts $num; incr next 100000}
		set line [split [gets $f] \t]
		set line [list_sub $line $poss]
		set line [downloaddb_dbsnp_convline $line]
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
	puts "Sorting $filename"
	cg select -s - $filename.temp $filename.temp2
	file rename -force $filename.temp2 $filename
	file delete $filename.temp
}
