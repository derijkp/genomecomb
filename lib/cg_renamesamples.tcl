#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc file_rename {file newfile} {
	if {$newfile eq $file} return
	puts "$file -> $newfile"
	file rename $file $newfile
}

proc renamesamples_newfilename {file changes} {
	set dir [file dir $file]
	set tail [file tail $file]
	if {![catch {dict get $changes $tail} new]} {
		return $dir/$new
	}
	set newfile $file
	dict for {name newname} $changes {
		if {[string match *-$name $file]} {
			set newfile [string range $file 0 [expr {[string length $file] - [string length $name] - 1}]]$newname
			break
		}
	}
	return $newfile
}

proc renamesamples_file {file changes} {
	set gzroot [gzroot $file]
	if {$gzroot ne $file} {
		set gzext [file extension $file]
	} else {
		set gzext ""
	}
	set ext [file extension $gzroot]
	set basefile [file root $gzroot]
	set newbasefile [renamesamples_newfilename $basefile $changes]
	if {[inlist {.tsv .sft .tab} $ext]} {
		set f [gzopen $file]
		set header [tsv_open $f comment]
		set newheader {}
		set changed 0
		foreach field $header {
			dict for {name newname} $changes {
				if {[string match *-$name $field]} {
					set field [string range $field 0 [expr {[string length $field] - [string length $name] - 1}]]$newname
					set changed 1
					break
				}
			}
			lappend newheader $field
		}
		if {$changed} {
			set o [open $newbasefile$ext.temp w]
			puts -nonewline $o $comment
			puts $o [join $newheader \t]
			fcopy $f $o
			close $o
			gzclose $f
			if {$gzext ne ""} {compress $newbasefile$ext.temp $gzext}
			file delete $file
			file_rename $newbasefile$ext.temp$gzext $newbasefile$ext$gzext
			puts "Adapted $newbasefile$ext$gzext"
		} else {
			gzclose $f
			file_rename $file $newbasefile$ext$gzext
		}
	} else {
		file_rename $file $newbasefile$ext$gzext
	}
}

proc renamesamples {dir changes} {
	foreach file [glob -nocomplain $dir/*] {
		if {[file isdir $file]} {
			renamesamples $file $changes
			set newfile [renamesamples_newfilename $file $changes]
			file_rename $file $newfile
		} else {
			renamesamples_file $file $changes
		}
	}
}

proc cg_renamesamples {dir args} {
	if {[llength $args] == 1} {
		set changes [lindex $args 0]
	} else {
		set changes $args
	}
	if {![file isdir $dir]} {
		renamesamples_file $dir $changes
	} else {
		cd [file dir $dir]
		set dir [file tail $dir]
		file delete -force $dir.temp
		exec cp -al $dir $dir.temp
		renamesamples $dir.temp $changes
		file delete -force $dir.old
		file rename $dir $dir.old
		file rename $dir.temp $dir
	}
}
