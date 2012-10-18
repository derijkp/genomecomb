#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg_regextract {args} {
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
	if {[llength $args] < 2} {
		errorformat regextract
		exit 1
	}
	foreach {cutoff} $args break
	set files [lrange $args 1 end]
	set files [lsort -dict $files]
	set o stdout
	puts $o "chromosome\tbegin\tend"
	foreach file $files {
		putslog "Processing $file"
		set chr [lindex [file root [split [gzroot $file] -]] 1]
		set ext [file extension $file]
		if {$ext eq ".bcol"} {
			set bcol [bcol_open $file]
			set start [lindex [dict get $bcol table] 0 0]
			set max [dict get $bcol max]
			set type [dict get $bcol type]
			set file [gzfile $file.bin]
			if {[inlist {.rz .gz .bgz} [file extension $file]]} {set cat zcat} else {set cat cat}
			set error [catch {
				# puts "$cat $file | getregionsbcol $chr $type $start $cutoff $above $shift"
				exec $cat $file | getregionsbcol $chr $type $start $cutoff $above $shift >@ $o
			} errmessage]
			if {$error} {
				set errmessage [split [string trim $errmessage] \n]
				set errmessage [list_sub $errmessage -exclude [list_find -glob $errmessage {*decompression OK, trailing garbage ignored*}]]
				if {[llength $errmessage]} {
					error [join $errmessage \n]
				}
			}
		} elseif {$ext eq ".bam"} {
			set chrcol 0
			set poscol 1
			set valuecol 2
			exec samtools depth $file | getregions "unkown" $chrcol $poscol $valuecol $cutoff $above $shift 0 >@ $o
		} else {
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
				exec $cat $file | getregions $chr $chrcol $poscol $qcol $cutoff $above $shift 1 >@ $o
			} errmessage]
			if {$error && ![regexp {decompression OK, trailing garbage ignored} $errmessage]} {
				error $errmessage
			}
		}
	}
}
