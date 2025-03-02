#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg_sortvcf {args} {
	set split 0
	set refseq {}
	set dbdir {}
	set threads 1
	cg_options vcf2tsv args {
		-threads {set threads $value}
	} {infile outfile} 0 2 {
		sort a vcf file, chromosome will be in natural sort order
	}
	if {[info exists infile]} {
		set f [gzopen $infile]
	} else {
		set f stdin
	}
	if {[info exists outfile]} {
		set o [wgzopen $outfile]
	} else {
		set o stdout
	}
	set header [tsv_open $f comment]
	set c [split [string range $comment 0 end-2] \n]
	set poss [list_find -regexp $c {^##contig=}]
	set list {}
	foreach line [list_sub $c $poss] {
		if {![regexp {contig=<ID=([^\t,]+)} $line temp chr]} {
			error "format error in contig line: $line"
		}
		lappend list [list $chr $line]
	}
	set list [bsort -sortchromosome $list]
	set neworder [list_subindex $list 1]
	set c [list_sub $c -exclude $poss]
	set insert [lindex $poss 0]
	if {$insert == [llength $c]} {
		lappend c {*}$neworder
	} else {
		set c [lreplace $c $insert -1 {*}$neworder]
	}
	puts $o [join $c \n]
	puts $o \#[join $header \t]
	chanexec $f $o "gnusort8 --parallel $threads -T \"[scratchdir]\" --buffer-size=500M --compress-program=zstd-mt-1 -t \\t -s --chromosome-sort"
	if {$o ne "stdout"} {catch {close $o}}
	if {$f ne "stdin"} {catch {gzclose $f}}
}
