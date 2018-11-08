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
	set cmd [list | csv2tsv]
	if {$tsvfile ne ""} {
		lappend cmd > $tsvfile
	} else {
		lappend cmd >@ stdout
	}
	set o [open $cmd w]
	fcopy $f $o
	if {$o ne "stdout"} {catch {close $o}}
	if {$f ne "stdin"} {catch {gzclose $f}}
}
