#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc multitranscript_open {isoformfiles aVar} {
	upvar $aVar a
	# prepare: open files, parse headers, load first line
	set fnum -1
	foreach file $isoformfiles {
		incr fnum
		catch {close $a(f,$fnum)}
	}
	unset -nocomplain a
	set header1 {}
	set common {}
	set basefields {chromosome begin end strand exonStarts exonEnds cdsStart cdsEnd transcript gene geneid}
	set fnum -1
	foreach file $isoformfiles {
		incr fnum
		set a(f,$fnum) [gzopen $file]
		set header [tsv_open $a(f,$fnum) a(comment,$fnum)]
		set poss [list_sub [tsv_basicfields $header 14 0] {0 1 2 6 7 8 9 10 11 12 13}]
		foreach pos [lrange $poss 0 2] field {chromosome begin end} {
			if {$pos == -1} continue
			lset header $pos $field
		}
		if {[inlist [lrange $poss 0 5] -1]} {
			set missing [list_sub $basefields [lrange [list_find $poss -1] 0 11]]
			error "file $file is missing essential fields: $missing"
		}
		if {[lindex $header [lindex $poss 9]] eq "gene_name" && [lindex $header [lindex $poss 10]] eq "geneid"} {
			# for isoquant where gene_name is sometimes empty, but gene isn't
			set pos [lsearch $header gene]
			if {$pos != -1} {
				lset poss 9 $pos
			}
		}
		set a(h,$fnum) $header
		set a(id,$fnum) $poss
		set a(data,$fnum) [list_find -glob $header *-*]
		set a(empty,$fnum) [list_fill [llength $a(data,$fnum)] 0.0]
		set a(status,$fnum) [gets $a(f,$fnum) a(curline,$fnum)]
		set a(curline,$fnum) [split $a(curline,$fnum) \t]
		set a(curid,$fnum) [list_sub $a(curline,$fnum) $poss]
		set a(prev,$fnum) [list_sub $a(curline,$fnum) $poss]
		set other [list_sub $header -exclude [list {*}$poss {*}$a(data,$fnum)]]
		set other [list_lremove $other $basefields]
		if {$header1 eq ""} {
			set header1 [list_sub $header -exclude $a(data,$fnum)]
			set common $other
		} else {
			set common [list_common $common $other]
		}
	}
	# find what fields need to go in the header
	set fnum -1
	foreach file $isoformfiles {
		incr fnum
		set a(common,$fnum) [list_cor $a(h,$fnum) $common]
	}
	set header $basefields
	lappend header {*}$common
	set fnum -1
	foreach file $isoformfiles {
		incr fnum
		lappend header {*}[list_sub $a(h,$fnum) $a(data,$fnum)]
	}
	return $header
}

proc multitranscript_directmatch {aVar o isoformfiles} {
	upvar $aVar a
	while 1 {
		set curids {}
		set fnum -1
		foreach file $isoformfiles {
			incr fnum
			lappend curids [lrange $a(curid,$fnum) 0 end-5]
		}
		set curid [lindex [bsort $curids] 0]
		if {$curid eq {{} {} {} {} {} {}}} break
		set pos [lsearch $curids $curid]
		set mfile [lindex $isoformfiles $pos]
		set line $a(curid,$pos)
		lappend line {*}[list_sub $a(curline,$pos) $a(common,$pos)]
		set fnum -1
		foreach file $isoformfiles {
			incr fnum
			if {[lrange $a(curid,$fnum) 0 end-5] eq $curid} {
				lappend line {*}[list_sub $a(curline,$fnum) $a(data,$fnum)]
				set a(status,$fnum) [gets $a(f,$fnum) a(curline,$fnum)]
				set a(curline,$fnum) [split $a(curline,$fnum) \t]
				set temp [list_sub $a(curline,$fnum) $a(id,$fnum)]
				set ctemp [lrange $temp 0 end-5]
				if {[lindex [bsort [list $curid $ctemp]] 1] ne $ctemp} {
					error "file $file not sorted correctly; should be sorted on: chromosome begin end strand exonStarts exonEnds"
				}
				set a(curid,$fnum) $temp
			} else {
				lappend line {*}$a(empty,$fnum)
			}
		}
		puts $o [join $line \t]
	}
}

proc multitranscript_approxmatch {aVar o isoformfiles match} {
	upvar $aVar a
	while 1 {
		# find earliest transcript to start region
		set curids {}
		set fnum -1
		foreach file $isoformfiles {
			incr fnum
			if {$a(status,$fnum) == -1} continue
			lappend curids [lrange $a(curid,$fnum) 0 end-5]
		}
		if {![llength $curids]} break
		set curid [lindex [bsort $curids] 0]
		foreach {curchr curbegin curend} $curid break
		# load all transcripts (from all files) that overlap
		# loop until no files have overlapping transcripts
		unset -nocomplain seta
		unset -nocomplain setma
		while 1 {
			set num 0
			set fnum -1
			foreach file $isoformfiles {
				incr fnum
				if {$a(status,$fnum) == -1} continue
				while 1 {
					set line $a(curid,$fnum)
					foreach {chromosome begin end strand exonStarts exonEnds cdsStart cdsEnd transcript gene_name geneid} $line break
					if {$chromosome ne $curchr} break
					if {$begin >= $curend} break
					incr num
					if {$end > $curend} {set curend $end}
					set common [list_sub $a(curline,$fnum) $a(common,$fnum)]
					set data [list_sub $a(curline,$fnum) $a(data,$fnum)]
					lappend line $common $data $fnum
#putsvars line transcript
##if {[regexp transcript $transcript]} error
#if {$transcript eq "transcript_chr9_27560423-e1226i731e96i3054e60i1086e488i6266e56"} {error tr}
#if {$exonStarts eq "27560423,27562380,27565530,27566676,27573430,"} {error starts}
#
					if {[regexp $match $transcript]} {
						set id [list $strand]
						set starts [split [string trim $exonStarts ,] ,]
						if {[llength $starts] == 1} {
							set id single
						} else {
							foreach s [lrange $starts 1 end] \
								e [lrange [split [string trim $exonEnds ,] ,] 0 end-1] {
								lappend id $e $s
							}
						}
						lappend setma($id) $line
					} else {
						set id [list $begin $end $strand $exonStarts $exonEnds]
						lappend seta($id) $line
					}
					set a(status,$fnum) [gets $a(f,$fnum) a(curline,$fnum)]
					set a(curline,$fnum) [split $a(curline,$fnum) \t]
					set temp [list_sub $a(curline,$fnum) $a(id,$fnum)]
					set ctemp [lrange $temp 0 end-5]
					if {[lindex [bsort [list $curid $ctemp]] 1] ne $ctemp} {
						error "file $file not sorted correctly; should be sorted on: chromosome begin end strand exonStarts exonEnds"
					}
					set a(curid,$fnum) $temp
					if {$a(status,$fnum) == -1} break
				}
			}
			if {!$num} break
		}
#set sid "- 27561649 27562380 27562476 27565530 27565590 27566676 27567164 27573430"
#if {[info exists setma($sid)]} {puts "in setma"}
		if {[info exists setma(single)]} {
			set list [bsort $setma(single)]
			if {[llength $list] == 1} {
				set line [lindex $list 0]
				foreach {chromosome begin end strand exonStarts exonEnds} $line break
				set id [list $begin $end $strand $exonStarts $exonEnds]
				lappend seta($id) $line
			} else {
				set line [lindex $list 0]
				set result [list $line]
				foreach {curchr curbegin curend} $line break
				set todo [lrange $list 1 end]
				lappend todo {}
				foreach line $todo {
					foreach {chr begin end} $line break
					if {$begin >= $curend || ![llength $line]} {
						set id [list $curbegin $curend $strand $curbegin, $curend,]
						foreach l $result {
							lset l 1 $curbegin
							lset l 2 $curend
							lset l 4 $curbegin,
							lset l 5 $curend,
							lappend seta($id) $l
						}
						set result [list $line]
						foreach {curchr curbegin curend} $line break
					} else {
						if {$end > $curend} {set curend $end}
						lappend result $line 
					}
				}
			}
			unset setma(single)
		}
		foreach id [array names setma] {
			set list $setma($id)
			if {[llength $list] > 1} {
				set begin [lmath_min [list_subindex $list 1]]
				set end [lmath_max [list_subindex $list 2]]
				foreach {strand starts ends} [list_sub [lindex $list 0] {3 4 5}] break
				set starts [join [lreplace [split $starts ,] 0 0 $begin] ,]
				set ends [join [lreplace [split $ends ,] end-1 end-1 $end] ,]
				lset list 0 4 $starts
				lset list 0 5 $ends
				set id [list $begin $end $strand $starts $ends]
				set seta($id) $list
			} else {
				set line [lindex $list 0]
				set id [list_sub $line {1 2 3 4 5}]
				lappend seta($id) $line
			}
		}
		set ids [array names seta]
		if {![llength $ids]} break
		foreach id [bsort $ids] {
			foreach {begin end strand starts ends} $id break
#set sid {27560423 27573766 - 27560423,27562380,27565530,27566676,27573430, 27561649,27562476,27565590,27567164,27573766,}
#if {$id eq $sid} {error seta}
			# if {$starts eq "1253911,1256044,1256991,1257207,1263345,1267861,1273665,"} {error error2}
			set ts $seta($id)
			unset -nocomplain va
			foreach line $ts {
				foreach {common data fnum} [lrange $line end-2 end] break
				set va($fnum) $data
			}
			set line [lindex $ts 0]
			set common [lindex $line end-2]
			set line [lrange $line 0 end-3]
			lset line 1 $begin
			lset line 2 $end
			lset line 4 $starts
			lset line 5 $ends
			lappend line {*}$common
			set fnum -1
			foreach file $isoformfiles {
				incr fnum
				if {[info exists va($fnum)]} {
					lappend line {*}$va($fnum)
				} else {
					lappend line {*}$a(empty,$fnum)
				}
			}
			puts $o [join $line \t]
		}
	}
}

proc cg_multitranscript {args} {
	set match {}
	cg_options multitranscript args {
		-match {
			set match $value
		}
	} compar_file 2
	set isoformfiles $args

	set header [multitranscript_open $isoformfiles a]
	# open result file
	catch {close $o}
	set o [open $compar_file.temp w]
	if {$a(comment,0) ne ""} {puts -nonewline $o $a(comment,0)}
	puts $o [join $header \t]
	if {$match eq ""} {
		# plain, direct matching (transcripts have to be identical)
		multitranscript_directmatch a $o $isoformfiles
	} else {
		# use approximate matching (where applicable)
		# for multi-exon transcripts with the same introns will be merged
		# for mono-exon transcripts overlapping will be merged
		multitranscript_approxmatch a $o $isoformfiles $match
	}
	close $o
	foreach file $isoformfiles {
		catch {close $a(f,$file)}
	}
	file rename -force $compar_file.temp $compar_file

}
