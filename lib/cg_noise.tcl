#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg_noise {args} {
	set resultfile -
	set mindepth 20
	set refseq {}
	set filtered 0
	set region {}
	set aa 0
	set debug -1
	cg_options noise args {
		-refseq {
			set refseq $value
		}
		-mindepth {set mindepth $value}
		-q {set q $value}
		-Q {set Q $value}
		-all {set aa $value}
		-f - -filtered {
			set filtered 1
			if {![info exists q]} {set q 20}
			if {![info exists Q]} {set Q 20}
		}
		-debug {
			set debug $value
		}
	} {bamfile resultfile} 1 ... {
		measure noise level (percentage alt allele) in bam file
	}
	set informat [ext2format $bamfile bam {bam cram sam}]
	set refseq [refseq $refseq]
	set samtoolsargs {}
	if {[info exists q]} {lappend samtoolsargs -q $q}
	if {[info exists Q]} {lappend samtoolsargs -Q $Q}
	if {$aa} {lappend samtoolsargs -aa}
	set o [wgzopen $resultfile]
	catch_exec samtools mpileup -f $refseq -A --ignore-overlaps {*}$samtoolsargs $bamfile | noise $mindepth $debug >@ $o
	gzclose $o
}
