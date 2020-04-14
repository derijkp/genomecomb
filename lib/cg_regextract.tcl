#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg_regextract {args} {
	set qfields {coverage uniqueSequenceCoverage}
	set posfields {offset pos position begin start}
	set above 0; set shift {}
	set filtered 0
	set region {}
	set aa 0
	set refseq {}
	set min 0
	set max {}
	cg_options regextract args {
		-min {set min $value}
		-max {set max $value}
		-shift {set shift $value}
		-qfields {set qfields $value}
		-posfields {set posfields $value}
		-q {set q $value}
		-Q {set Q $value}
		-all {set aa 1}
		-f - -filtered {
			set filtered 1
			if {![info exists q]} {set q 20}
			if {![info exists Q]} {set Q 20}
		}
		-region {
			set region $value
		}
		-refseq {
			set refseq $value
		}
	} {} 1
	set samtoolsargs {}
	if {[info exists q]} {lappend samtoolsargs -q $q}
	if {[info exists Q]} {lappend samtoolsargs -Q $Q}
	if {$aa} {lappend samtoolsargs -aa}
	set files $args
	set files [bsort $files]
	set o stdout
	puts $o "chromosome\tbegin\tend"
	foreach file $files {
		putslog "Processing $file"
		set ext [file extension $file]
		if {$shift ne ""} {set useshift $shift} else {
			if {$ext in ".bam .cram .sam"} {set useshift -1} else {set useshift 0}
		}
		if {$ext eq ".bcol"} {
			if {$region ne ""} {error "option -region only supported for bam, cram and sam files"}
			set bcol [bcol_open $file]
			if {[dict get $bcol version] == 0} {
				set chr [lindex [file root [split [gzroot $file] -]] 1]
				set start [lindex [dict get $bcol table] 0 0]
				if {$max eq ""} {
					set max [dict get $bcol max]
				}
				set type [dict get $bcol type]
				set file [gzfile $file.bin]
				set cat [gzcat $file]
				set error [catch {
					# puts "$cat $file | getregionsbcol $chr $type $start $min $max $useshift"
					exec {*}$cat $file | getregionsbcol $chr $type $start $min $max $useshift >@ $o
				} errmessage]
				if {$error} {
					set errmessage [split [string trim $errmessage] \n]
					set errmessage [list_sub $errmessage -exclude [list_find -glob $errmessage {*decompression OK, trailing garbage ignored*}]]
					if {[llength $errmessage]} {
						error [join $errmessage \n]
					}
				}
			} else {
				# puts "getregionsbcol2 $file $min $max $useshift"
				exec getregionsbcol2 $file $min $max $useshift >@ $o
			}
		} elseif {$ext in ".bam .cram .sam"} {
			set chrcol 0
			set poscol 1
			if {$region ne ""} {
				if {$ext ni ".bam .cram"} {error "option -region only supported for bam or cram files"}
			}
			set index $file.[indexext $file]
			if {![file exists $index]} {
				exec samtools index $file
			}
			set regions [samregions $region $refseq]
			if {![llength $regions]} {set regions {{}}}
			foreach region $regions {
				set opts $samtoolsargs
				if {$region ne ""} {
					lappend opts -r $region
				}
				if {[catch {
					if {!$filtered} {
						set valuecol 2
						catch_exec samtools depth {*}$opts $file | getregions "unkown" $chrcol $poscol $valuecol $min $max $useshift 0 >@ $o
					} else {
						set valuecol 3
						catch_exec samtools mpileup --ignore-overlaps {*}$opts $file | getregions "unkown" $chrcol $poscol $valuecol $min $max $useshift 0 >@ $o
					}
				} msg]} {
					if {$region ne "" && [regexp "samtools depth: can't parse region \"$region\"" $msg]} {
						putslog "warning: $msg"
					} else {
						error $msg
					}
				}
			}
		} else {
			if {$region ne ""} {error "option -region only supported for bam, cram and sam files"}
			set f [gzopen $file]
			set header [tsv_open $f]
			catch {gzclose $f}
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
			} else {
				set chr [lindex [file root [split [gzroot $file] -]] 1]
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
			set cat [gzcat $file]
			set error [catch {
				exec {*}$cat $file | getregions $chr $chrcol $poscol $qcol $min $max $useshift 1 >@ $o
			} errmessage]
			if {$error && ![regexp {decompression OK, trailing garbage ignored} $errmessage] && ![regexp {Successfully decoded} $errmessage]} {
				error $errmessage
			}
		}
	}
}
