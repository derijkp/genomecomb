#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

package require Extral

proc cg_tsv2csv {args} {
	set tsvfile {}
	set csvfile {}
	cg_options tsv2csv args {
	} {tsvfile csvfile} 0 2 {
		Converts comma-separated value (csv) data to tab-sparated (tsv)
	}
	if {$tsvfile ne ""} {
		set f [gzopen $tsvfile]
	} else {
		set f stdin
	}
	set cmd [list | tsv2csv]
	if {$csvfile ne ""} {
		lappend cmd > $csvfile
	} else {
		lappend cmd >@ stdout
	}
	set o [open $cmd w]
	fcopy $f $o
	if {$o ne "stdout"} {catch {close $o}}
	if {$f ne "stdin"} {catch {gzclose $f}}
}
