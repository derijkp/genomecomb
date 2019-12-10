#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

package require Extral

proc collapseoverlap_join {cur scorepos {numpos -1}} {
	if {[llength $cur] == 1} {return [lindex $cur 0]}
	if {$scorepos != -1} {
		set cur [ssort -natural -decreasing -index $scorepos $cur]
	}
	set result {}
	set len [llength [lindex $cur 0]]
	for {set i 0} {$i < $len} {incr i} {
		if {$i == $scorepos} {
			lappend result [lindex $cur 0 $i]
		} elseif {$i == $numpos} {
			lappend result [format %0.f [lmath_sum [list_subindex $cur $i]]]
		} else {
			set temp [list_subindex $cur $i]
			set nodup [lsort -dict [list_remdup $temp]]
			if {[llength $nodup] < 2} {
				lappend result $nodup
			} else {
				lappend result [join $nodup ,]
			}
		}
	}
	return $result
}

proc collapseoverlap {{infile stdin} {resultfile stdout} {scorefield score} {numfield num}} {
	catch {close $f} ; catch {close $o}
	if {$infile ne "stdin"} {
		putslog "making $resultfile"
		if {[catch {open $infile} f]} {
			error "Could not open $infile"
		}
	} else {
		set f stdin
	}
	set cor [open_region $f header]
	foreach {chrpos startpos endpos} $cor break
	set scorepos [lsearch $header $scorefield]
	set numpos [lsearch $header $numfield]
	if {$resultfile ne "stdout"} {
		putslog "making $resultfile"
		if {[catch {open $resultfile.temp w} o]} {
			error "Could not write outputfile $resultfile"
		}
	} else {
		set o stdout
	}
	puts $o [join $header \t]
	set line [split [gets $f] "\t"]
	foreach {chr start end} [list_sub $line $cor] break
	set cur [list $line]
	set num 0; set next 100000
	set prev {}
	while {1} {
		incr num; if {$num >= $next} {putslog $num; incr next 100000}
		set read [gets $f line]
		set line [split $line "\t"]
		set cchr {}
		foreach {cchr cstart cend} [list_sub $line $cor] break
		set newchr [expr {$cchr ne $chr}]
		if {$newchr || ($cstart > $start)} {
			# write preceeding
			set stop 0
			set keepend -1
			foreach l $cur {
				set end [lindex $l $endpos]
				if {!$newchr && ($end > $cstart)} {
					set end $cstart
					set stop 1
				}
				if {$end != $keepend} {
					set joined [collapseoverlap_join $cur $scorepos $numpos]
					lset joined $startpos $start
					lset joined $endpos $end
					if {[llength $prev]} {
						if {($start != [lindex $prev $endpos]) || ([list_sub $joined -exclude $cor] ne $prevsub)} {
							puts $o [join $prev \t]
							set prev $joined
							set prevsub [list_sub $prev -exclude $cor]
						} else {
							lset prev $endpos $end
						}
					} else {
						set prev $joined
						set prevsub [list_sub $prev -exclude $cor]
					}
					# puts $o [join $joined \t]
					set start $end
				}
				if {$stop} break
				set keepend $end
				list_shift cur
			}
		}
		if {$read == -1} break
		# check sort
		if {!$newchr && $cstart < $start} {error "file incorrectly sorted: $chr: $cstart < $start"}
		# add new line
		lappend cur $line
		if {[llength $cur] > 1} {
			set cur [lsort -real -index $endpos $cur]
		}
		set chr $cchr
		set start $cstart
	}
	puts $o [join $prev \t]
	close $f
	if {$resultfile ne "stdout"} {
		close $o
		file rename -force -- $resultfile.temp $resultfile
		putslog "Finished $resultfile"
	}
}

proc cg_regcollapse {args} {
	set pos 0
	set resultfile stdout
	set scorefield score
	set numfields num
	cg_options regcollapse args {
		-s - -scorefield {
			set scorefield $value
		}
		-n - -numfields {
			set numfields $value
		}
		-o - -output {
			set resultfile $value
		}
	}
	if {$resultfile ne "stdout"} {
		if {[file exists $resultfile]} {
			puts "$resultfile already exists"
			exit 1
		}
	}
	putslog "Collapsing file(s) to $resultfile "
	set tempfile [tempfile]
	set tempfile2 [tempfile]
	if {[llength $args] > 1} {
		putslog "Concatenating files"
		cg cat -m {*}$args > $tempfile
		putslog "Sorting"
		set f [gzopen $tempfile]
		set header [tsv_open $f]
		set bposs [tsv_basicfields $header 3]
		gzclose $f
		set sfields [list_sub $header $bposs]
		lappend sfields {*}[list_sub $header -exclude $bposs]
		cg select -s $sfields $tempfile $tempfile2
	} elseif {[llength $args] == 0} {
		set tempfile2 stdin
	} else {
		set tempfile2 [lindex $args 0]
	}
	putslog "Collapsing"
	collapseoverlap $tempfile2 $resultfile $scorefield $numfields
	file delete $tempfile
	if {$tempfile2 ne [lindex $args 0]} {
		file delete $tempfile2
	}
}
