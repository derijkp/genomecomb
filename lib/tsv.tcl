proc tsv_select_idtopos {header id fields} {
	set poss [list_cor $header $fields]
	if {[inlist $poss -1]} {error "sample $id not present"}
	return [lmath_calc $poss + 1]
}

proc tsv_select_sm {header ids} {
	set id1 [list_pop ids]
	foreach {a11 a21 ref} [tsv_select_idtopos $header $id1 [list alleleSeq1-$id1 alleleSeq2-$id1 reference]] break
	set temp [list "((\$$a11 != \"-\")||(\$$a21 != \"-\")) && (\$$a11 != \$$ref || \$$a21 != \$$ref)"]
	foreach id $ids {
		foreach {a12 a22} [tsv_select_idtopos $header $id1 [list alleleSeq1-$id alleleSeq2-$id]] break
		lappend temp "((\$$a11 == \$$a12 && \$$a21 == \$$a22) || (\$$a11 == \$$a22 && \$$a21 == \$$a12))"
	}
	set temp "([join $temp " && "])"
}

proc tsv_select_df {header ids} {
	set temp1 {}
	set temp2 {}
	foreach id $ids {
		foreach {a1 a2 ref} [tsv_select_idtopos $header $id [list alleleSeq1-$id alleleSeq2-$id reference]] break
		lappend temp1 "((\$$a1 != \"-\" && \$$a1 != \$$ref) || (\$$a2 != \"-\" && \$$a2 != \$$ref))"
		lappend temp2 "(\$$a1 == \$$ref && \$$a2 == \$$ref)"
	}
	set temp "(([join $temp1 " || " ]) && ([join $temp2 " || "]))"
}

proc tsv_select_mm {header ids} {
	set temp1 {}
	set temp2 {}
	set list {}
	foreach id $ids {
		foreach {a1 a2 ref} [tsv_select_idtopos $header $id [list alleleSeq1-$id alleleSeq2-$id reference]] break
		lappend list [list $a1 $a2]
	}
	while {[llength $list]} {
		foreach {a1 a2} [list_pop list] break
		lappend temp1 "((\$$a1 != \"-\" && \$$a1 != \$$ref) || (\$$a2 != \"-\" && \$$a2 != \$$ref))"
		list_foreach {a12 a22} $list {
			lappend temp2 "\$$a1 != \$$a12 || \$$a2 != \$$a22"
		}
	}
	set temp "(([join $temp1 " && " ]) && ([join $temp2 " || "]))"
}

proc tsv_select_un {header ids} {
	set temp1 {}
	set temp2 {}
	foreach id $ids {
		foreach {a1 a2 ref} [tsv_select_idtopos $header $id [list alleleSeq1-$id alleleSeq2-$id reference]] break
		lappend temp1 "(\$$a1 != \"-\" && (\$$a1 != \$$ref || \$$a2 != \$$ref))"
		lappend temp2 "\$$a1 == \"-\""
	}
	set temp "(([join $temp1 " || " ]) && ([join $temp2 " || "]))"
}

proc tsv_select {query {qfields {}} {sortfields {}} {f stdin} {out stdout}} {
	fconfigure $f -buffering none
	fconfigure $out -buffering none
	set header [tsv_open $f]
	set awk ""
	set sort ""
	set cut ""
	if {[llength $sortfields]} {
		set poss [list_cor $qfields $sortfields]
		if {[lsearch $poss -1] != -1} {
			set poss [list_cor $header $sortfields]
			if {[lsearch $poss -1] != -1} {error "fields [join [list_sub $sortfields [list_find $poss -1]] ,] not found"}
			if {$qfields ne ""} {
				set qposs [list_cor $header $qfields]
				set qposs [lmath_calc $qposs + 1]
				set cut "cut -d \\t -f [join $qposs ,]"
			}
		}
		set poss [lmath_calc $poss + 1]
		set keys {}
		foreach pos $poss {
			lappend keys $pos,$pos
		}
		set sort "gnusort8 -t \\t -V -s -k[join $keys " -k"]"
	}
	if {$query ne ""} {
		set indices [list_unmerge [regexp -all -indices -inline {[$]([a-zA-z0-9_(),-]+)} $query]]
		set indices [list_reverse $indices]
		list_foreach {start end} $indices {
			set field [string range $query [expr {$start+1}] $end]
			if {[regexp {^(.*)\((.*)\)$} $field temp func args]} {
				switch $func {
					sm {
						set ids [split $args ,]
						set temp [tsv_select_sm $header $ids]
						set query [string_replace $query $start $end $temp]
					}
					df {
						set ids [split $args ,]
						set temp [tsv_select_df $header $ids]
						set query [string_replace $query $start $end $temp]
					}
					mm {
						set ids [split $args ,]
						set temp [tsv_select_mm $header $ids]
						set query [string_replace $query $start $end $temp]
					}
					un {
						set ids [split $args ,]
						set temp [tsv_select_un $header $ids]
						set query [string_replace $query $start $end $temp]
					}
					default {
						error "Unkown function $func"
					}
				}
			} else {
				set pos [lsearch $header $field]
				if {$pos == -1} {error "field \"$field\" not present"}
				incr pos
				set query [string_replace $query $start $end \$$pos]
			}
		}
		set awk {BEGIN {FS="\t" ; OFS="\t"} }
		append awk $query
		if {($qfields ne "") && ($cut eq "")} {
			set qposs [list_cor $header $qfields]
		} else {
			set qposs [list_cor $header $header]
		}
		set qposs [lmath_calc $qposs + 1]
		append awk " \{print $[join $qposs ,$]\}"
	} elseif {($qfields ne "") && ($cut eq "")} {
		set qposs [list_cor $header $qfields]
		set qposs [lmath_calc $qposs + 1]
		set cut "cut -d \\t -f [join $qposs ,]"
	}
	set pipe {}
	if {$awk ne ""} {
		lappend pipe [list awk $awk]
	}
	if {$sort ne ""} {
		lappend pipe $sort
	}
	if {$cut ne ""} {
		lappend pipe $cut
	}
	# putslog pipe:[join $pipe " | "]
	if {$qfields ne ""} {puts $out [join $qfields \t]} else {puts $out [join $header \t]}
	if {![llength $pipe]} {
		fcopy $f $out
	} else {
		set o [open "|\ [join $pipe " | "]\ >@\ $out" w]
		fcopy $f $o
		close $o
	}
}

if 0 {
	lappend auto_path ~/dev/completegenomics/lib
	package require Tclx
	signal -restart error SIGINT
	package require Extral
	set f [open GS102/ASM/var-GS000000078-ASM.tsv]
	set query "\$begin < 2000"
	set qfields "chromosome begin end"
	set sortfields "haplotype"
	set sortfields "chromosome begin"

cg select -q '' < GS102/ASM/var-GS000000078-ASM.tsv | less
cg select -f 'haplotype chromosome begin' < GS102/ASM/var-GS000000078-ASM.tsv | less
cg select -q '$begin < 2000' < GS102/ASM/var-GS000000078-ASM.tsv | less

cg select -q '$begin < 2000' < GS102/ASM/var-GS000000078-ASM.tsv > /tmp/test1
cg select -q '$begin < 2000' -f 'chromosome begin end' -s 'chromosome begin' < GS102/ASM/var-GS000000078-ASM.tsv > /tmp/test2
cg select -q '$begin < 2000' -f 'chromosome begin end' -s 'haplotype' < GS102/ASM/var-GS000000078-ASM.tsv > /tmp/test3

}

proc tsv_sort {filename fields} {
	set f [open $filename]
	set line [gets $f]
	close $f
	if {[string index $line 0] eq "#"} {set line [string range $line 1 end]}
	set header [split $line \t]
	set poss [list_cor $header $fields]
	if {[lsearch $poss -1] != -1} {error "fields [join [list_sub $fields [list_find $poss -1]] ,] not found"}
	set poss [lmath_calc $poss + 1]
	puts [join $header \t]
	set command "tail +2 [list $filename] | gnusort8 -t \\t -V -s -k[join $poss " -k"] >@ stdout"
	eval exec $command
}

proc file_rootgz {filename} {
	if {[file extension $filename] eq ".gz"} {
		return [file root [file root $filename]]
	} elseif {[file extension $filename] eq ".rz"} {
		return [file root [file root $filename]]
	} else {
		file root [file root $filename]
	}
}

proc tsv_open {f} {
	set keep 0
	while {![eof $f]} {
		set line [gets $f]
		if {![string length $line]} continue
		if {[string index $line 0] ne "#"} break
		set keep [tell $f]
		set header $line
	}
	if {[string index $line 0] eq ">"} {
		return [split [string range $line 1 end] \t]
	}
	if {![info exists header]} {
		return [split $line \t]
	} elseif {[string index $header 0] eq "#"} {
		seek $f $keep
		return [split [string range $header 1 end] \t]
	} else {
		seek $f $keep
		return [split $header \t]
	}
}

if 0 {
	lappend auto_path ~/dev/completegenomics/lib
	package require Tclx
	signal -restart error SIGINT
	package require Extral
	cd /complgen/compar

	set filename /data/db/ucsc_ori/_data_db_ucsc-exapted_repeats.tsv
	set f [open $filename]
	set fields {chrom chromStart chromEnd}

	# set f [open /complgen/compar/78vs79_compar-filter-sc.tsv]
	set f [open /complgen/compar/78vs79_compar_pvt.tsv]
	set query {compar df sample "|79 78,79" refcons "" ns "" lowscore "" trf "" str "" rp "" sd "" sc "" dbsnp "" loc EXON}

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

proc tsv_index {file xfield} {
	if {[inlist {.rz} [file extension $file]]} {
		set indexname [file root $file].${xfield}_index
		set tempfile [tempfile]
		exec razip -d -c $file > $tempfile
		set f [open $tempfile]
	} else {
		set indexname $file.${xfield}_index
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
	catch {Classy::Progress start [file size $file] "Making index"}
	while {![eof $f]} {
		set fpos [tsv_next $f $xpos $next]
		if {[eof $f]} {
			set xmax $fpos
			break
		}
		lappend index $fpos
		incr prev 10000
		incr next 10000
		catch {Classy::Progress set [tell $f]}
		if {![expr $next%1000000]} {putslog $next}
	}
	catch {Classy::Progress stop}
	close $f
	set size [file size $file]
	set f [open $file]
	seek $f [expr {$size-5000}] start
	gets $f
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set xmax [lindex $line $xpos]
	}
	close $f
	set o [open $indexname w]
	puts $o 10000
	puts $o $findex
	puts $o $xmin
	puts $o $xmax
	puts $o [join $index \n]
	close $o
}

proc tsv_index_open {file field {uncompress 0}} {
	global cache
	set file [file normalize $file]
	if {[info exists cache(tsv_index,$file,$field,step)]} return
	set ext [file extension $file]
	if {$ext eq ".gz"} {set uncompress 1}
	set root [rzroot $file]
	set workfile $file
	set uncompressed 0
	set remove 0
	if {[inlist {.rz .gz} $ext]} {
		if {$uncompress} {
			set workfile $root
			puts "temporarily uncompressing $file"
			catch {exec gunzip -c $file > $workfile}
			set uncompressed 1
			set remove 1
		}
	} else {
		set uncompressed 1
	}
	set indexname $root.${field}_index
	if {![file exists $indexname]} {
		tsv_index $file $field
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
	close $o	lappend auto_path ~/dev/completegenomics/lib
	lappend auto_path /complgen/bin/complgen/apps/cg/lib
	package require Extral
	package require Tclx
	signal -restart error SIGINT
set compar_file /complgen/multicompar/compar-X.tsv

	set f [rzopen $workfile]
	set cache(tsv_index,$file,header) [tsv_open $f]
	if {$uncompressed} {
		set cache(tsv_index,$file,$field,channel) $f
	} else {
		catch {close $f}
	}
}

proc tsv_index_close {file field} {
	global cache
	set file [file normalize $file]
	if {$cache(tsv_index,$file,$field,remove) && ($cache(tsv_index,$file,$field,workfile) ne $file)} {
		close $cache(tsv_index,$file,$field,channel)
		# puts "remove $cache(tsv_index,$file,$field,workfile)"
		file remove $cache(tsv_index,$file,$field,workfile)
	}
	set indexname [rzroot $file].${field}_index
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
	set file [file normalize $file]
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
		set f [rzopen $file $fpos]
	}
	return $f
}

proc tsv_index_get {file field pos} {
	global cache
	set file [file normalize $file]
	set header $cache(tsv_index,$file,header)
	set fieldpos [lsearch $header $field]
	set uncompressed $cache(tsv_index,$file,$field,uncompressed)
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
