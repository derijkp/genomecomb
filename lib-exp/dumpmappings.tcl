#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

if 0 {

cd /data/peter/complgen
set libfile CGI-DATA-sample/GS00028-DNA-C01/MAP/GS000005323-MAP-sample/lib_DNB_GS00433-CLS.tsv
set file CGI-DATA-sample/GS00028-DNA-C01/MAP/GS000005323-MAP-sample/mapping-100000.tsv
set rfile CGI-DATA-sample/GS00028-DNA-C01/MAP/GS000005323-MAP-sample/reads-100000.tsv
set resultsdir mappings 

lappend auto_path ~/dev/genomecomb/lib
package require Extral

proc getmategap {libfile} {
	set f [opencgifile $libfile header]
	while {![eof $f]} {
		set line [gets $f]
		foreach {id type armId indArm objArm min max} $line break
		if {$type eq "mategap"} {return $min}
	}
	return {}
}

set mategap [getmategap $libfile]
# majority = 419, lmath_majoritydev = 8
# expected (cgi): 450
set gapmin 300
set gapmax 600

set f [opencgifile $file header]
set os [open $resultsdir/single.tsv w]
set oc [open $resultsdir/chromosome.tsv w]
set od [open $resultsdir/dir.tsv w]
set og [open $resultsdir/gap.tsv w]
array set numa {1 0 2 0 m 0 c 0 d 0 g 0}
set diffs {}
set list {}
while {![eof $f]} {
	set line [gets $f]
	foreach {flags chromosome offsetInChr gap1 gap2 gap3 weight mateRec} [split $line \t] break
	set end [expr {$offsetInChr + 35 + $gap1 + $gap2 + $gap3}]
	binary scan $weight c weight
	incr weight -33
	set lastdnbrecord [expr {$flags & 0x01}]
	if {[expr {$flags & 0x02}]} {set side r} else {set side l}
	if {[expr {$flags & 0x04}]} {set strand -} else {set strand +}
	lappend list [list $flags $chromosome $offsetInChr $end $strand $side $gap1 $gap2 $gap3 $weight $mateRec]
	if {$lastdnbrecord} {
		set len [llength $list]
		if {$len > 2} {incr numa(m)} else {incr numa($len)}
		if {$len == 1} {
			puts $os [join [list_sub [lindex $list 0] {1 2 3 4}] \t]
		}
		if {($len != 2) || ([lindex $list 0 5] eq [lindex $list 1 5])} {
			set list {}
			continue
		}
		set list [lsort -integer -index 2 $list]
		set result [list_concat [list_sub [lindex $list 0] {1 2 3 4}] [list_sub [lindex $list 1] {1 2 3 4}]]
		foreach {chr1 s1 e1 str1 chr2 s2 e2 str2} $result break
		set diff [expr {$s2 - $e1}]
		lappend diffs $diff
		if {$chr1 ne $chr2} {
			incr numa(c)
			puts $oc [join $result \t]
		} elseif {$str1 ne $str2} {
			incr numa(d)
			puts $od [join $result \t]
		} elseif {($diff < $gapmin) || ($diff > $gapmax)} {
			incr numa(g)
			puts $og [join $result \t]
		}
		set list {}
	}
}
close $f
close $os
close $oc
close $od
close $og

foreach temp {chromosome dir gap single} {
	exec sort mappings/$temp.tsv > mappings/$temp-s.tsv
}

# GASV output
#Left Chromosome:	The chromosome where the left read mapped. (Acceptable formats: CHR1, chr1, 1, X, chrX, chrY)
#Left Mapped Start:	The coordinate of the start of the left read.
#Left Mapped End:	The coordinate of the end of the left read.
#Left End Orientation:	The orientation of the left read. (Acceptable formats: include +/-, Plus/Minus)
#Right Chromosome:	Chromosome where the right read mapped.
#Right Mapped Start: 	Coordinate of the start of the right read.
#Right Mapped End:	Coordinate of the end of the right read.
#Right End Orientation:	The orientation of the right read.
}
