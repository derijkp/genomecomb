#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg_estimate_error_rate {args} {
	set samtoolsargs {}
	set resultfile -
	set mindepth 20
	set perpos 0
	set refseq {}
	set filtered 0
	set region {}
	set q 0
	set Q 0
	set aa 0
	set debug 0
	set ignoreoverlaps 1
	set countorphans 1
	set progress 1000000
	set typediffs all
	cg_options estimate_error_rate args {
		-refseq {
			set refseq $value
		}
		-mindepth {set mindepth $value}
		-maxdepth {lappend samtoolsargs -d $value}
		-perpos {set perpos [true $value]}
		-q - -min-MQ {set q $value}
		-Q - -min-BQ {set Q $value}
		-f - -filtered {
			set filtered 1
			if {![info exists q]} {set q 20}
			if {![info exists Q]} {set Q 20}
		}
		-exclude {
			set exclude $value
		}
		-include {
			set include $value
		}
		-ignore-overlaps {
			set ignoreoverlaps $value
		}
		-count-orphans {
			set countorphans $value
		}
		-debug {
			set debug $value
		}
		-progress {
			set progress $value
		}
	} {bamfile regionfile varfile resultfile} 3 ... {
		estimate error rate (percentage) in bam file
	}
	set informat [ext2format $bamfile bam {bam cram sam}]
	set refseq [refseq $refseq]
	if {[info exists exclude] || [info exists include]} {
		if {[info exists exclude]} {
			set exclude [list_remove [split $exclude ",; "] {}]
		} else {
			set exclude {UNMAP SECONDARY QCFAIL DUP}
		}
		set exclude [list_lremove $exclude [split $include ",; "]]
		set excludeflag [sam_filter $exclude]
		lappend samtoolsargs --excl-flags $excludeflag
	}
	if {$ignoreoverlaps} {
		lappend samtoolsargs --ignore-overlaps
	}
	if {$countorphans} {
		lappend samtoolsargs --count-orphans
	}
	# lappend samtoolsargs --output-extra FLAG
	set regfile [tempfile].tsv.zst
	set tempvar [tempfile].tsv.zst
	cg select -s -  $varfile $tempvar
	cg regsubtract $regionfile $tempvar | cg zst -c 1 > $regfile
	set o [wgzopen $resultfile]
	catchchildkilled_exec samtools mpileup {*}$samtoolsargs --no-BAQ -f $refseq  \
		-q $q -Q $Q $bamfile | estimate_error_rate $regfile $debug $progress >@ $o 2>@ stderr
	gzclose $o
}
