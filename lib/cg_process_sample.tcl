#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc process_compare_checkfield {file field} {
	set f [open $file]
	set header [tsv_open $f]
	close $f
	inlist $header $field
}

proc process_indexcompress {file} {
	set ext [file extension $file]
	if {$ext eq ".gz"} {
		gunzip $file
		set file [file root $file]
	}
	set f [gzopen $file]
	set header [tsv_open $f]
	foreach field {offset end1} {
		set fpos [lsearch $header $field]
		if {$fpos != -1} break
	}
	if {$fpos == -1} {error "no column offset, end1 in file $file"}
	if {([gets $f] eq "") && [eof $f]} return
	close $f
	if {![file exists $file.${field}_index]} {
		tsv_index $field $file
	}
	putslog "Compressing $file"
	if {![inlist {.rz .bgz} $ext]} {
		exec bgzip -c $file > $file.gz.temp
		file rename -force $file.gz.temp $file.gz
		file delete $file
	}
}

proc cg_process_indexcompress {args} {
	global scriptname action
	if {[llength $args] != 1} {
		puts stderr "format is: $scriptname $action file"
		puts stderr " - makes index, and compresses to bgzip"
		exit 1
	}
	foreach {file} $args break
	process_indexcompress $file
}

proc cg_process_sample {args} {
	cgmakelib_process_sample
	global sample
	if {([llength $args] < 2) || ([llength $args] > 3)} {
		errorformat process_sample ?force?
		exit 1
	}
	set force 0
	foreach {dir destdir force} $args break
	set dir [file normalize $dir]
	set destdir [file normalize $destdir]
	file mkdir $destdir
	catch {file delete $destdir/oricg}
	if {[file exists $dir]} {
		mklink $dir $destdir/oricg
	}
	set sample [file tail $destdir]
	catch {file delete $destdir/cg_process_sample.finished}
	cgmake --force $force $destdir/cg_process_sample-$sample.finished
}
