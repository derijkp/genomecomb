#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

package require Extral

proc cg_csv2tsv {args} {
	set csvfile {}
	set tsvfile {}
	cg_options csv2tsv args {
	} {csvfile tsvfile} 0 2 {
		Converts comma-separated value (csv) data to tab-sparated (tsv)
	}
	if {$csvfile ne ""} {
		set f [gzopen $csvfile]
	} else {
		set f stdin
	}
	if {$tsvfile ne ""} {
		set o [open $tsvfile w]
	} else {
		set o stdout
	}
	while {![eof $f]} {
		set line [gets $f]
		if {[string index $line 0] eq "\#"} {
			puts $o $line
		} else {
			puts $o [join [csv_split $line] \t]
		}
	}
	if {$o ne "stdout"} {catch {close $o}}
	if {$f ne "stdin"} {catch {gzclose $f}}
}
