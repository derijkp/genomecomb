#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc tsv_open {f {keepheaderVar {}}} {
	if {$keepheaderVar ne ""} {
		upvar $keepheaderVar keepheader
	}
	set keepheader {}
	set keep 0
	set buffering [fconfigure $f -buffering]
	fconfigure $f -buffering line
	set split 1
	set line [gets $f]
	set fchar [string index $line 0]
	set fchar2 [string index $line 1]
	if {$fchar eq "#" && $fchar2 eq "#"} {
		set vcf 1
	} else {
		set vcf 0
	}
	while {![eof $f]} {
		if {![string length $line]} {
			lappend keepheader \#
			set line [gets $f]
			break
		}
		set fchar [string index $line 0]
		if {$fchar eq ">"} {
			break
		} elseif {$fchar ne "#"} {
			break
		} elseif {$vcf} {
			set fchar2 [string index $line 1]
			if {$fchar2 ne "#"} {
				set split 0
				break
			}
		}
		lappend keepheader $line
		set line [gets $f]
	}
	fconfigure $f -buffering $buffering
	set fchar [string index $line 0]
	if {[inlist {# >} $fchar]} {
		set keepheader [join $keepheader \n]\n
		if {!$split} {
			return [string range $line 1 end]
		} else {
			return [split [string range $line 1 end] \t]
		}
	} else {
		if {[llength $keepheader]} {
			set keepheader [join $keepheader \n]\n
		}
		return [split $line \t]
	}
}

proc tsv_next {f xpos next {shift 100000}} {
	# do a ~ binary search to get at next faster
	set start [tell $f]
	while 1 {
		seek $f $shift current
		gets $f
		set line [getnotempty $f]
		set x [lindex $line $xpos]
		if {![isdouble $x] || ($x >= $next)} {
			seek $f $start
			set shift [expr {$shift / 2}]
			if {$shift < 1000} break
		}
		set start [tell $f]
	}
	while {![eof $f]} {
		set fpos [tell $f]
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set x [lindex $line $xpos]
		if {$x >= $next} break
	}
	if {![eof $f]} {
		return $fpos
	} else {
		return $x
	}
}

proc tsv_nextline {f xpos next {shift 100000}} {
	# do a ~ binary search to get at next faster
	set start [tell $f]
	while 1 {
		seek $f $shift current
		gets $f
		set line [getnotempty $f]
		set x [lindex $line $xpos]
		if {![isdouble $x] || ($x >= $next)} {
			seek $f $start
			set shift [expr {$shift / 2}]
			if {$shift < 1000} break
		}
		set start [tell $f]
	}
	while {![eof $f]} {
		set fpos [tell $f]
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set x [lindex $line $xpos]
		if {$x >= $next} break
	}
	return $line
}

proc tsv_index {xfield file} {
	set indexname [gzroot $file].${xfield}_index
	if {[inlist {.rz} [file extension $file]]} {
		set tempfile [scratchfile]
		exec razip -d -c $file > $tempfile
		set f [open $tempfile]
	} elseif {[inlist {.bgz} [file extension $file]]} {
		set tempfile [scratchfile]
		exec bgzip -d -c $file > $tempfile
		set f [open $tempfile]
	} elseif {[inlist {.gz} [file extension $file]]} {
		set tempfile [scratchfile]
		exec gunzip -c $file > $tempfile
		set f [open $tempfile]
	} else {
		set f [open $file]
	}
	set header [tsv_open $f]
	set xpos [lsearch $header $xfield]
	if {$xpos == -1} {error "field $xfield not present in file $file"}
	set fstart [tell $f]
	set fpos $fstart
	set line [split [gets $f] \t]
	set xmin [lindex $line $xpos]
	set findex [expr {$xmin-$xmin%10000}]
	set prev $findex
	set next [expr {$prev + 10000}]
	set index [list $fpos]
	catch {progress start [file size $file] "Making index"}
	while {![eof $f]} {
		set fpos [tsv_next $f $xpos $next]
		if {[eof $f]} {
			set xmax $fpos
			break
		}
		lappend index $fpos
		incr prev 10000
		incr next 10000
		catch {progress set [tell $f]}
		if {![expr $next%1000000]} {putslog $next}
	}
	catch {progress stop}
	close $f
	set size [file size $file]
	set f [open $file]
	seek $f [expr {$size-5000}] start
	gets $f
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set temp [lindex $line $xpos]
		if {[isint $temp]} {
			set xmax $temp
		}
	}
	close $f
	set o [open $indexname.temp w]
	puts $o 10000
	puts $o $findex
	puts $o $xmin
	puts $o $xmax
	puts $o [join $index \n]
	close $o
	file rename -force $indexname.temp $indexname
	if {[info exists tempfile]} {
		file delete $tempfile
	}
}

proc cg_tsv_index {args} {
	global scriptname action
	if {[llength $args] < 2} {
		puts stderr "format is: $scriptname $action field tsvfile ..."
		exit 1
	}
	set field [list_shift args]
	foreach tsvfile $args {
		putslog "Indexing $tsvfile"
		tsv_index $field $tsvfile
	}
}

proc tsv_index_header {file} {
	global cache
	set file [file_absolute $file]
	return [get cache(tsv_index,$file,header)]
}

proc tsv_index_open {file field {uncompress 0}} {
	global cache
	set file [file_absolute $file]
	if {[info exists cache(tsv_index,$file,$field,step)]} return
	set ext [file extension $file]
	if {$ext eq ".gz" || $ext eq ".rz" || $ext eq ".bgz"} {set uncompress 1}
	set root [gzroot $file]
	set workfile $file
	set uncompressed 0
	set remove 0
	if {[inlist {.rz .gz .bgz} $ext]} {
		if {$uncompress} {
			set workfile [scratchfile]
			puts stderr "temporarily uncompressing $file to $workfile"
			file delete -force $workfile.temp
			file delete -force $workfile
			gunzip $file $workfile.temp
			file rename -force $workfile.temp $workfile
			set uncompressed 1
			set remove 1
		}
	} else {
		set uncompressed 1
	}
	set indexname $root.${field}_index
	if {![file exists $indexname]} {
		tsv_index $field $workfile
		if {$workfile ne $file} {
			file rename -force [gzroot $workfile].${field}_index $indexname
		}
	}
	set o [open $indexname]
	set cache(tsv_index,$file,$field,step) [gets $o]
	set cache(tsv_index,$file,$field,findex) [gets $o]
	set xmin [gets $o]
	set cache(tsv_index,$file,$field,xmin) $xmin
	set cache(tsv_index,$file,$field,xmax) [gets $o]
	set cache(tsv_index,$file,$field,fx) [expr {$xmin-$xmin%10000}]
	set cache(tsv_index,$file,$field,index) [split [string trim [read $o]] \n]
	set cache(tsv_index,$file,$field,workfile) $workfile
	set cache(tsv_index,$file,$field,uncompressed) $uncompressed
	set cache(tsv_index,$file,$field,remove) $remove
	close $o
	set f [gzopen $workfile]
	set cache(tsv_index,$file,header) [tsv_open $f]
	if {$uncompressed} {
		set cache(tsv_index,$file,$field,channel) $f
	} else {
		catch {close $f}
	}
}

proc tsv_index_close {file field} {
	global cache
	set file [file_absolute $file]
	if {[info exists cache(tsv_index,$file,$field,workfile)] && [get cache(tsv_index,$file,$field,remove) 0] && ($cache(tsv_index,$file,$field,workfile) ne $file)} {
		close $cache(tsv_index,$file,$field,channel)
		# puts "remove $cache(tsv_index,$file,$field,workfile)"
		file delete $cache(tsv_index,$file,$field,workfile)
	}
	set indexname [gzroot $file].${field}_index
	unset -nocomplain cache(tsv_index,$file,$field,step)
	unset -nocomplain cache(tsv_index,$file,$field,findex)
	unset -nocomplain cache(tsv_index,$file,$field,xmin)
	unset -nocomplain cache(tsv_index,$file,$field,xmax)
	unset -nocomplain cache(tsv_index,$file,$field,fx)
	unset -nocomplain cache(tsv_index,$file,$field,index)
	unset -nocomplain cache(tsv_index,$file,$field,workfile)
	unset -nocomplain cache(tsv_index,$file,$field,uncompressed)
	unset -nocomplain cache(tsv_index,$file,$field,remove)
	unset -nocomplain cache(tsv_index,$file,$field,channel)
	unset -nocomplain cache(tsv_index,$file,header)
}

proc tsv_index_apprstop {file field} {
	global cache
	set uncompressed $cache(tsv_index,$file,$field,uncompressed)
	if {!$uncompressed} {
		catch {close $f}
	}
}

proc tsv_index_apprgoto {file field pos} {
	global cache
	set file [file_absolute $file]
	set index $cache(tsv_index,$file,$field,index)
	set step $cache(tsv_index,$file,$field,step)
	set start [expr {round($pos)-round($pos)%$step}]
	set indexpos [expr {($start-$cache(tsv_index,$file,$field,findex))/$step}]
	if {$indexpos < 0} {set indexpos 0}
	if {$indexpos >= [llength $index]} {set indexpos end}
	set fpos [expr {round([lindex $index $indexpos])}]
	set uncompressed $cache(tsv_index,$file,$field,uncompressed)
	if {$uncompressed} {
		set f $cache(tsv_index,$file,$field,channel)
		seek $f $fpos
	} else {
		set f [gzopen $file $fpos]
	}
	return $f
}

proc tsv_index_get {file field pos} {
	global cache
	set file [file_absolute $file]
	set header $cache(tsv_index,$file,header)
	set fieldpos [lsearch $header $field]
	set uncompressed [get cache(tsv_index,$file,$field,uncompressed) 0]
	set f [tsv_index_apprgoto $file $field $pos]
	if {$uncompressed} {
		set line [tsv_nextline $f $fieldpos $pos]
	} else {
		set line {}
		while 1 {
			# read in mem in chunks
			# do a ~ binary search to get at target faster
			set chunk [read $f 20480]
			append chunk [gets $f]
			if {![string length $chunk]} break
			set table [split $chunk \n]
			set temp [lindex $table end $fieldpos]
			if {$temp < $pos} continue
			if {$temp == $pos} break
			set ipos [binsearch $table $fieldpos $pos]
			set line [lindex $table $ipos]
			break
		}
	}
	tsv_index_apprstop $file $field
	set temp [lindex $line $fieldpos]
	if {$temp != $pos} {error "$pos not found in $file,$field"}
	return $line
}

# -1 if comp1 < comp2
# -2 if comp1 > comp2
# 0 if comp1 = comp2
# > 0 if comp1 ~ comp2
proc tsv_align_match {comp1 comp2} {
	foreach e1 $comp1 e2 $comp2 {
		if {$e1 ne $e2} {
			if {$e2 eq ""} {
				return -1
			} elseif {$e1 eq ""} {
				return -2
			}
			if {[ssort -natural [list $e1 $e2]] eq [list $e1 $e2]} {return -1} else {return -2}
		}
	}
	return 0
}

proc tsv_basicfields {header {num 6} {giveerror 1}} {
	set poss [list_cor $header {chromosome begin end type ref alt}]
	set nfposs [list_find $poss -1]
	foreach nfpos $nfposs {
		switch $nfpos {
			0 {
				foreach name {chrom chr chr1 genoName tName contig} {
					set v [lsearch $header $name]
					if {$v != -1} break
				}
			}
			1 {
				foreach name {start end1 chromStart genoStart tStart txStart} {
					set v [lsearch $header $name]
					if {$v != -1} break
				}
			}
			2 {
				foreach name {start2 chromEnd genoEnd tEnd txEnd} {
					set v [lsearch $header $name]
					if {$v != -1} break
				}
			}
			4 {
				set v [lsearch $header reference]
			}
			5 {
				set v [lsearch $header alternative]
			}
			default {
				continue
			}
		}
		lset poss $nfpos $v
	}
	incr num -1
	set poss [lrange $poss 0 $num]
	set pos [lsearch $poss -1]
	if {$pos == 1} {
		set pos [lsearch $header pos]
		if {$poss == -1} {set pos [lsearch $header offset]}
		lset poss 1 $pos
		lset poss 2 $pos
		set pos [lsearch $poss -1]
	}
	if {$giveerror && ($pos != -1)} {
		set notfound [list_sub {chromosome begin end type ref alt} [list_find $poss -1]]
		error "header error: fields (or alternatives) not found: $notfound"
	}
	return $poss
}

proc tabix {file chromosome begin end} {
	# tabix thinks a variant x-$begin overlaps the interval $begin-$end
	# incr begin to stop these from (incorrectly) appearing
	incr begin
	set temp [split [exec tabix $file $chromosome:$begin-$end] \n]
}

proc cg_maketabix {args} {
	foreach file $args {
		set ext [file extension $file]
		if {$ext ne ".gz"} {
			cg_bgzip $file
			set file [gzroot $file].gz
		}
		if {[file exists $file.tbi]} {
			putslog "Skipping $file: $file.tbi exists"
			continue
		}
		putslog "making tabix for $file"
		set f [gzopen $file]
		set header [tsv_open $f comment]
		catch {close $f}
		set skip [llength [split $comment \n]]
		foreach {chrompos beginpos endpos} [lmath_calc [tsv_basicfields $header 3] + 1] break
		exec tabix -s $chrompos -b $beginpos -e $endpos -0 -S $skip $file
	}
}
