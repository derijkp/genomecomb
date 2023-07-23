#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

package require Extral

proc cg_bed2genetsv {args} {
	if {([llength $args] < 0) || ([llength $args] > 2)} {
		errorformat bed2genetsv
	}
	if {[llength $args] > 0} {
		set filename [lindex $args 0]
		set f [gzopen $filename]
	} else {
		set f stdin
	}
	if {[llength $args] > 1} {
		set outfile [lindex $args 1]
		set o [wgzopen $outfile]
	} else {
		set o stdout
	}
	while {![eof $f]} {
		set line [gets $f]
		if {[string index $line 0] eq "\#"} {
			puts $o $line
		} elseif {[inlist {browser track} [lindex $line 0]]} {
			puts $o #$line
		} else {
			break
		}
	}
	set extrafields {
		name2 cdsStartStat cdsEndStat exonFrames txId transcripttype geneName genetype geneid transcript
	}
	set bedfields {
		chromosome begin end name score strand thickStart thickEnd itemRgb
		blockCount blockSizes blockStarts
	}
	lappend bedfields {*}$extrafields
	if {[eof $f] && ![llength $line]} {
		puts $o [join {chromosome begin end} \t]
	} else {
		puts $o [join {chromosome begin end strand exonStarts exonEnds cdsStart cdsEnd transcript geneid gene transcripttype genetype} \t]
		while 1 {
			set line [split $line \t]
			set size [expr {[llength $line]-1}]
			foreach $bedfields $line break
			set exonStarts {}
			set exonEnds {}
			foreach start [split $blockStarts ,] size [split $blockSizes ,] {
				if {$start eq ""} continue
				lappend exonStarts [expr {$begin+$start}]
				lappend exonEnds [expr {$begin+$start+$size}]
			}
			set cdsStart $thickStart
			set cdsEnd $thickEnd
			if {$transcript eq ""} {set transcript $name}
			if {$name2 ne ""} {set gene $name2}
			if {$gene eq ""} {aset gene $geneid}
			puts $o [join [list $chromosome $begin $end $strand [join $exonStarts ,] [join $exonEnds ,] $cdsStart $cdsEnd $transcript $geneid $gene $transcripttype $genetype] \t]
			if {[gets $f line] == -1} break
		}
	}
	if {$o ne "stdout"} {catch {close $o}}
	if {$f ne "stdin"} {catch {gzclose $f}}
}

proc cg_bed2sft {args} {
	cg_bed2genetsv {*}$args
}

