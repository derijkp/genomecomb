#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

package require Extral

proc cg_bed2tsv {args} {
	if {([llength $args] < 0) || ([llength $args] > 2)} {
		errorformat bed2sft
	}
	if {[llength $args] > 0} {
		set filename [lindex $args 0]
		set f [gzopen $filename]
	} else {
		set f stdin
	}
	if {[llength $args] > 1} {
		set outfile [lindex $args 1]
		set o [open $outfile w]
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
	set bedfields {chromosome begin end name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts}
	if {[eof $f] && ![llength $line]} {
		puts $o [join {chromosome begin end} \t]
	} else {
		set size [expr {[llength $line]-1}]
		puts $o [join [lrange $bedfields 0 $size] \t]
		puts $o $line
		fconfigure $f -translation binary
		fconfigure $o -translation binary
		fcopy $f $o
	}
	if {$o ne "stdout"} {catch {close $o}}
	if {$f ne "stdin"} {catch {gzclose $f}}
}

proc cg_bed2sft {args} {
	cg_bed2tsv {*}$args
}
