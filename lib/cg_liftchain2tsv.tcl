#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg_liftchain2tsv {args} {
	global scriptname action
	set pos 0
	foreach {key value} $args {
		switch -- $key {
			-ref {
				set ref $value
			}
			-destref {
				set destref $value
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	if {[llength $args] > 2} {
		errorformat liftchain2tsv
	}
	foreach {srcfile destfile} $args break
	if {![info exists ref] || ![info exists destref]} {
		set liftoverfilebase [lindex [split [file tail $srcfile] .] 0]
		if {![regexp {^(.*)(To|2)(.*)} $liftoverfilebase temp oldrefname temp newrefname]} {
			error "Cannot deduce source and destination genome from filename (expected to be <srcgenome>To<destgenome>.*).\nUse -ref and -destref options"
		}
		if {![info exists ref]} {
			set ref [string tolower $oldrefname]
		}
		if {![info exists destref]} {
			set destref [string tolower $newrefname]
		}
	}
	if {$srcfile eq ""} {
		set f stdin
	} else {
		set f [gzopen $srcfile]
	}
	if {$destfile eq ""} {
		set o [open "| cg select -s - >@ stdout" w]
	} else {
		set o [open $destfile.temp w]
	}
	puts $o "#filetype\ttsv/liftover"
	puts $o "#fileversion\t[fileversion]"
	puts $o "#ref\t$ref"
	puts $o "#destref\t$destref"
	puts $o [join {chromosome begin end strand destchromosome destbegin destend deststrand} \t]
	while 1 {
		while 1 {
			if {[gets $f line] == -1} break
			if {[llength $line]} break
		}
		if {[eof $f]} break
		foreach {temp score tname tsize tstrand tstart tend qname qsize qstrand qstart qend id} $line break
		set tname [chr_clip $tname]
		set qname [chr_clip $tname]
		if {$temp ne "chain"} {error "error in chain format: expected chain at start of: $line"}
		if {$tstrand eq "-"} {
			error "error: did not expect negative strand in header for reference sequence"
		}
		if {$qstrand eq "+"} {
			while 1 {
				if {[gets $f line] == -1} break
				foreach {size dt dq} $line break
				puts $o [join [list $tname $tstart [expr {$tstart+$size}] $tstrand $qname $qstart [expr {$qstart+$size}] $qstrand] \t]
				if {[llength $line] == 1} {
					break
				}
				set tstart [expr {$tstart+$size+$dt}]
				set qstart [expr {$qstart+$size+$dq}]
			}
		} else {
			set qrstart [expr {$qsize-$qend}]
			set qrend [expr {$qsize-$qstart}]
			while 1 {
				if {[gets $f line] == -1} break
				foreach {size dt dq} $line break
				puts $o [join [list $tname $tstart [expr {$tstart+$size}] $tstrand $qname [expr {$qrend-$size}] $qrend $qstrand] \t]
				if {[llength $line] == 1} {
					break
				}
				set tstart [expr {$tstart+$size+$dt}]
				set qrend [expr {$qrend-$size-$dq}]
			}
		}
	}
	if {$srcfile ne ""} {
		gzclose $f
	}
	if {$destfile ne ""} {
		close $o
		cg select -s - $destfile.temp $destfile.temp2
		file rename -force $destfile.temp2 $destfile
		file delete $destfile.temp
	}
}

