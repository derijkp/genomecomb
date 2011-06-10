#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc ucsc2region {ucsc_file} {
	catch {close $f}
	set f [open $ucsc_file]
	set temp [string range [gets $f] 1 end]
	set header [split $temp \t]
	set poss [list_cor $header {chrom chromStart chromEnd name}]
	puts [join [list chromosome start end name] \t]
	while {![eof $f]} {
		set line [getnotempty $f]
		foreach {chrom chromStart chromEnd name} [list_sub $line $poss] break
		puts [join [list $chrom $chromStart $chromEnd $name] \t]
	}
	close $f
}

proc cg_ucsc2region {args} {
	global scriptname action
	if {[llength $args] != 1} {
		puts stderr "format is: $scriptname $action ucsc_file"
		puts stderr "convert ucsc tab file to regionfile for further use in filtering"
		exit 1
	}
	foreach {ucsc_file} $args break
	ucsc2region $ucsc_file
}

if 0 {
# ------------------------------------------------------------------------------
lappend auto_path ~/dev/completegenomics/lib
package require Extral
package require Tclx
signal -restart error SIGINT

cd /data/db
set ucsc_file _data_db_ucsc-simple_repeats.tsv

}
