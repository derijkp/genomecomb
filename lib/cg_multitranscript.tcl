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
		lappend poss [lsearch $header category]
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
	lappend basefields category
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

proc cg_multitranscript {args} {
	set match {}
	set exact 0
	set skipempty 1
	cg_options multitranscript args {
		-exact {
			set exact 1
			set match {}
		}
		-match {
			set match $value
		}
		-skipempty {
			set skipempty $value
		}
	} compar_file 2
	set isoformfiles $args
	if {$skipempty} {
		set todo {}
		foreach isoformfile $isoformfiles {
			if {[file size $isoformfile] == 0} continue
			lappend todo $isoformfile
		}
		set isoformfiles $todo
		if {![llength $isoformfiles]} {
			error "multitranscript error: all given isoformfiles are empty"
		}
	}

	analysisinfo_combine $compar_file $isoformfiles
	foreach file $isoformfiles {
		catch {close $a(f,$file)}
	}
	set header [multitranscript_open $isoformfiles a]
	# open result file
	catch {close $o}
	set o [wgzopen $compar_file.temp[gzext $compar_file]]
	if {$a(comment,0) ne ""} {puts -nonewline $o $a(comment,0)}
	puts $o [join $header \t]
	# process by groups of overlapping transcripts
	while 1 {
		# find earliest transcript to start region
		set curids {}
		set fnum -1
		foreach file $isoformfiles {
			incr fnum
			if {$a(status,$fnum) == -1} continue
			lappend curids [lrange $a(curid,$fnum) 0 5]
		}
		if {![llength $curids]} break
		set curid [lindex [bsort $curids] 0]
		foreach {curchr curbegin curend} $curid break
		# load all transcripts (from all files) that overlap
		# loop until no files have overlapping transcripts
		# seta will contain all known, setma all novel models, key is strand+junctions
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
					if {[string index $exonStarts end] eq ","} {
						set exonStarts [string range $exonStarts 0 end-1]
						lset line 4 $exonStarts
					}
					if {[string index $exonEnds end] eq ","} {
						set exonEnds [string range $exonEnds 0 end-1]
						lset line 5 $exonEnds
					}
					if {$chromosome ne $curchr} break
					if {$begin >= $curend} break
					incr num
					if {$end > $curend} {set curend $end}
					set common [list_sub $a(curline,$fnum) $a(common,$fnum)]
					set data [list_sub $a(curline,$fnum) $a(data,$fnum)]
					lappend line $common $data $fnum
					set cat [lindex $line 11]
					if {$exact} {
						set amatch 0
					} elseif {$match eq {}} {
						if {$cat eq "known"} {set amatch 0} else {set amatch 1}
					} elseif {[regexp $match $transcript]} {
						set amatch 1
					} else {
						set amatch 0
					}
					if {!$amatch} {
						set id [list $begin $end $strand $exonStarts $exonEnds]
						lappend seta($id) $line
					} else {
						set id [list $strand]
						set starts [split $exonStarts ,]
						set ends [split $exonEnds ,]
						if {[llength $starts] == 1} {
							set id single
						} else {
							foreach s [lrange $starts 1 end] \
								e [lrange $ends 0 end-1] {
								lappend id $e $s
							}
						}
						lappend setma($id) $line
					}
					set a(status,$fnum) [gets $a(f,$fnum) a(curline,$fnum)]
					set a(curline,$fnum) [split $a(curline,$fnum) \t]
					set temp [list_sub $a(curline,$fnum) $a(id,$fnum)]
					set ctemp [lrange $temp 0 5]
					if {[lindex [bsort [list $curid $ctemp]] 1] ne $ctemp} {
						error "file $file not sorted correctly; should be sorted on: chromosome begin end strand exonStarts exonEnds"
					}
					set a(curid,$fnum) $temp
					if {$a(status,$fnum) == -1} break
				}
			}
			if {!$num} break
		}
		if {[info exists setma(single)]} {
			# hande single exon transcripts
			# join $setma(single) \n
			# check for matches to known
			set knowns [array names seta]
			if {[llength $knowns]} {
				# if it matches known transcript -> add there
				set list {}
				foreach line [bsort $setma(single)] {
					foreach {chromosome begin end strand exonStarts exonEnds} $line break
					set matched 0
					foreach kline $knowns {
						foreach {kbegin kend kstrand kexonStarts kexonEnds} $kline break
						if {[llength $kexonStarts] > 1} continue
						if {$begin >= $kbegin && $end <= $kend} {
							lappend seta($kline) $line
							set matched 1
							break
						}
					}
					if {!$matched} {
						lappend list $line
					}
				}
			} else {
				set list [bsort $setma(single)]
			}
			# merge
			if {[llength $list] == 0} {
				#everything from setma(single) moved to known -> do nothing
			} elseif {[llength $list] == 1} {
				set line [lindex $list 0]
				foreach {chromosome begin end strand exonStarts exonEnds} $line break
				set id [list $begin $end $strand $exonStarts $exonEnds]
				lappend seta($id) $line
			} else {
				# merge overlapping (take min begin and max end)
				set line [lindex $list 0]
				set result [list $line]
				foreach {curchr curbegin curend} $line break
				set todo [lrange $list 1 end]
				lappend todo {}
				foreach line $todo {
					foreach {chr begin end} $line break
					if {$begin >= $curend || ![llength $line]} {
						set id [list $curbegin $curend $strand $curbegin $curend]
						foreach l $result {
							lset l 1 $curbegin
							lset l 2 $curend
							lset l 4 $curbegin
							lset l 5 $curend
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
				# calculate min begin and max end for merged transcript
				set begin [lmath_min [list_subindex $list 1]]
				set end [lmath_max [list_subindex $list 2]]
				foreach {strand starts ends} [list_sub [lindex $list 0] {3 4 5}] break
				set starts [join [lreplace [split $starts ,] 0 0 $begin] ,]
				set ends [join [lreplace [split $ends ,] end end $end] ,]
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
			set ts $seta($id)
			unset -nocomplain va
			foreach line $ts {
				foreach {common data fnum} [lrange $line end-2 end] break
				set va($fnum) $data
			}
			if {[llength $ts] > 1} {
				set cats [list_subindex $ts 11]
				set p [lsearch $cats known]
				if {$p == -1} {
					set p [lsearch $cats {}]
				}
				if {$p == -1} {
					set line [lindex $ts 0]
					set name [iso_name [lindex $line 0] $strand $starts $ends]
					lset line 8 $name
				} else {
					set line [lindex $ts $p]
				}
			} else {
				set line [lindex $ts 0]
				set cat [lindex $line 11]
				if {$cat ni {known {}}} {
					set name [iso_name [lindex $line 0] $strand $starts $ends]
					lset line 8 $name
				}
			}
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

	close $o
	foreach file $isoformfiles {
		catch {close $a(f,$file)}
	}
	file rename -force $compar_file.temp[gzext $compar_file] $compar_file
}
