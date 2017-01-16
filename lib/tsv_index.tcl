#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc tsv_index {xfield file} {
	set indexname [gzroot $file].${xfield}_index
	set f [gzopen $file]
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
	gzclose $f
	set size [file size $file]
	set f [gzopen $file]
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
	gzclose $f
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
		error "format is: $scriptname $action field tsvfile ..."
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
	if {$ext eq ".gz" || $ext eq ".rz" || $ext eq ".lz4" || $ext eq ".bgz"} {set uncompress 1}
	set root [gzroot $file]
	set workfile $file
	set uncompressed 0
	set remove 0
	if {[inlist {.rz .gz .bgz} $ext]} {
		if {$uncompress} {
			set workfile [scratchfile]
			putslog "temporarily uncompressing $file to $workfile"
			file delete -force $workfile.temp
			file delete -force $workfile
			gunzip $file $workfile.temp
			file rename -force $workfile.temp $workfile
			set uncompressed 1
			set remove 1
		}
	} elseif {[inlist {.lz4} $ext]} {
		if {$uncompress} {
			set workfile [scratchfile]
			putslog "temporarily uncompressing $file to $workfile"
			file delete -force $workfile.temp
			file delete -force $workfile
			exec lz4c -q -d $file $workfile.temp
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
		catch {gzclose $f}
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
		catch {gzclose $f}
		set skip [llength [split $comment \n]]
		foreach {chrompos beginpos endpos} [lmath_calc [tsv_basicfields $header 3] + 1] break
		exec tabix -s $chrompos -b $beginpos -e $endpos -0 -S $skip $file
	}
}
