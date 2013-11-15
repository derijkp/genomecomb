#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

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
			set nodup [list_remdup $temp]
			if {[llength $nodup] < 2} {
				lappend result $nodup
			} else {
				lappend result [join $nodup ,]
			}
		}
	}
	return $result
}

proc collapseoverlap {file {resultfile stdout} {scorefield score} {numfield num}} {
	catch {close $f} ; catch {close $o}
	if {[catch {open $file} f]} {
		error "Could not open file $file"
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
	while {![eof $f]} {
		incr num; if {$num >= $next} {putslog $num; incr next 100000}
		set line [split [gets $f] "\t"]
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
		file rename -force $resultfile.temp $resultfile
		putslog "Finished $resultfile"
	}
}

proc cg_collapseoverlap {args} {
	set pos 0
	set resultfile {}
	foreach {key value} $args {
		switch -- $key {
			-o {
				set resultfile $value
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	if {([llength $args] < 1)} {
		puts stderr "format is: $::base file ..."
		puts stderr " - Collapses overlapping regions in a region file."
		puts stderr " - makes a new file with reg_ prepended to the original filename"
		puts stderr " - Removal of overlap can be done by taking only the highest"
		puts stderr " - scoring region (this is always done when score is available)"
		puts stderr " - or taking all regions in 1 line (if score is not available)"
		puts stderr " - for a field with then name num all values will be added"
		puts stderr " - use the -o option to collapse multiple files into one new file (the filename given with the option -o)"
		exit 1
	}
	if {$resultfile ne ""} {
		if {[file exists $resultfile]} {
			puts "Skipping $resultfile: already exists"
			exit 0
		}
		puts "Collapsing multiple files into $resultfile"
		set ffile [lindex $args 0]
		set args [lrange $args 1 end]
		set f [gzopen $ffile]
		set header [tsv_open $f]
		set bposs [tsv_basicfields $header 3]
		close $f
		foreach file $args {
			set f [gzopen $file]
			set tempheader [tsv_open $f]
			close $f
			if {$tempheader ne $header} {
				error "for collapsing multiple files into one (-o option), all files must have the same columns"
			}
		}
		file copy -force $ffile $resultfile.temp1
		foreach file $args {
			exec tail -n +2 $file >> $resultfile.temp1
		}
		puts "Sorting"
		set sfields [list_sub $header $bposs]
		cg select -s $sfields $resultfile.temp1 $resultfile.temp2
		puts "Collapsing"
		collapseoverlap $resultfile.temp2 $resultfile
		file delete $resultfile.temp1 $resultfile.temp2
	} else {
		foreach {path} $args break
		puts "----------------------------------------------------"
		foreach file $args {
			set path [file dir $file]
			set tail [file tail $file]
			if {[string range $tail 0 4] eq "ucsc_"} {set tail [string range $tail 5 end]}
			set resultfile ${path}/reg_$tail
			if {[file exists $resultfile]} {
				puts "Skipping $resultfile: already exists"
				continue
			}
			collapseoverlap $file $resultfile
		}
	}
}

proc cg_regcollapse {args} {
	set pos 0
	set resultfile stdout
	set scorefield score
	set numfields num
	foreach {key value} $args {
		switch -- $key {
			-s {
				set scorefield $value
			}
			-n {
				set numfields $value
			}
			-o {
				set resultfile $value
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	if {([llength $args] < 1)} {
		errorformat reg_collapse
		exit 1
	}
	if {$resultfile ne "stdout"} {
		if {[file exists $resultfile]} {
			puts "$resultfile already exists"
			exit 0
		}
	}
	putslog "Collapsing file(s) to $resultfile"
	set tempfile [tempfile]
	set tempfile2 [tempfile]
	if {[llength $args]} {
		putslog "Concatenating files"
		cg cat -m {*}$args > $tempfile
		putslog "Sorting"
		set f [gzopen $tempfile]
		set header [tsv_open $f]
		set bposs [tsv_basicfields $header 3]
		close $f
		set sfields [list_sub $header $bposs]
		cg select -s $sfields $tempfile $tempfile2
	} else {
		set tempfile2 [lindex $args 0]
	}
	putslog "Collapsing"
	collapseoverlap $tempfile2 $resultfile $scorefield $numfields
	file delete $tempfile $tempfile2
}

if {[info exists argv0] && [file tail [info script]] eq [file tail $argv0]} {
	package require pkgtools
	set appdir [file dir [pkgtools::startdir]]
	lappend auto_path $appdir/lib
	append env(PATH) :[file dir [file dir $appdir]]/bin:$appdir/bin
	package require Extral
	set ::base $scriptname
	cg_collapseoverlap {*}$argv
}


