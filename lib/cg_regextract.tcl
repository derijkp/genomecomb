#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg_regextract {args} {
	if {[llength $args] < 4} {
		errorformat regextract
		exit 1
	}
	set qfields {coverage uniqueSequenceCoverage}
	set posfields {offset pos position begin start}
	set above 0; set shift 0
	set pos 0
	foreach {key value} $args {
		switch -- $key {
			-above {set above $value}
			-shift {set shift $value}
			-qfields {set qfields $value}
			-posfields {set posfields $value}
			-- break
			default {break}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	foreach {cutoff} $args break
	set files [lrange $args 1 end]
	set files [lsort -dict $files]
	set o stdout
	puts $o "chromosome\tbegin\tend"
	foreach file $files {
		putslog "Processing $file"
		set chr [lindex [file root [split [gzroot $file] -]] 1]
		set f [gzopen $file]
		set header [tsv_open $f]
		catch {close $f}
		foreach field {chromosome chrom} {
			set chrcol [lsearch $header $field]
			if {$chrcol != -1} break
		}
		if {$chrcol != -1} {
			set line [split [get $f] \t]
			set chr [lindex $line $chrcol]
			catch {close $f}
			set f [gzopen $file]
			set header [tsv_open $f]
		}
		foreach field $posfields {
			set poscol [lsearch $header $field]
			if {$poscol != -1} break
		}
		if {$poscol == -1} {
			error "no position column (one of [join $posfields ,]) found in $file"
		}
		foreach field $qfields {
			set qcol [lsearch $header $field]
			if {$qcol != -1} break
		}
		if {$qcol == -1} {
			error "no query column (one of [join $qfields ,]) found in $file"
		}
		if {[inlist {.rz .gz .bgz} [file extension $file]]} {set cat zcat} else {set cat cat}
		set error [catch {
			exec $cat $file | getregions $chr $poscol $qcol $cutoff $above $shift >@ $o
		} errmessage]
		if {$error && ![regexp {decompression OK, trailing garbage ignored} $errmessage]} {
			error $errmessage
		}
	}
}
