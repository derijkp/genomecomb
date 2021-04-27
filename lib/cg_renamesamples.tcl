#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc file_rename {file newfile} {
	set time [file_mtime $file]
	if {$newfile eq $file} return
	puts "$file -> $newfile"
	file rename $file $newfile
	file mtime $newfile $time
}

proc renamesamples_newfilename {file changes} {
	set dir [file dir $file]
	set tail [file tail $file]
	if {![catch {dict get $changes $tail} new]} {
		return $dir/$new
	}
	set newfile $file
	dict for {name newname} $changes {
		if {[string match *-$name $tail]} {
			set newfile $dir/[string range $tail 0 [expr {[string length $tail] - [string length $name] - 1}]]$newname
			break
		}
		set pos [string first -$name. $tail]
		if {$pos != -1} {
			incr pos -1
			set newfile $dir/[string range $tail 0 $pos]-$newname.[string range $tail [expr {$pos + [string length $name] + 3}] end]
			break
		}
	}
	return $newfile
}

proc renamesamples_file {file changes {relink 0}} {
	set gzroot [gzroot $file]
	set gzext [gzext $file]
	set ext [file extension $gzroot]
	set basefile [file root $gzroot]
	set newbasefile [renamesamples_newfilename $basefile $changes]
	if {$newbasefile eq $basefile} {
		# try first . (if more than 2 exts)
		set tbasefile [file root $basefile]
		while {$tbasefile ne $basefile} {
			set basefile $tbasefile
			set tbasefile [file root $basefile]
		}
		set newbasefile [renamesamples_newfilename $basefile $changes]
		set ext [string range [gzroot $file] [string length $basefile] end]
	}
	set newfile $newbasefile$ext$gzext
	if {![catch {file link $file} link]} {
		set newlink [renamesamples_newfilename $link $changes]
		if {$relink && $link ne $newlink && ![regexp ^/ $link]} {
			puts "relinking $file to $newfile"
			file delete $file
			mklink $newlink $newfile
		} elseif {$file ne $newfile} {
			file_rename $file $newfile
		}
	} elseif {[inlist {.tsv .sft .tab} $ext]} {
		puts "converting $file to $newfile"
		set time [file_mtime $file]
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
		set pos [lsearch $header sample]
		if {$pos != -1} {
			set changed 1
		}
		if {$changed} {
			set tempfile [filetemp_ext $newfile 0]
			set o [wgzopen $tempfile]
			puts -nonewline $o $comment
			puts $o [join $newheader \t]
			if {$pos == -1} {
				fcopy $f $o
			} else {
				while {[gets $f line] != -1} {
					set split [split $line \t]
					set name [lindex $split $pos]
					set nsplit [split $name -]
					set sample [lindex $nsplit end]
					if {[dict exists $changes $sample]} {
						set pre [join [lrange $nsplit 0 end-1] -]
						if {$pre ne ""} {append pre "-"}
						lset split $pos $pre[dict get $changes $sample]
					}
					puts $o [join $split \t]
				}
			}
			close $o
			gzclose $f
			file delete $file
			file_rename $tempfile $newfile
			file mtime $newfile $time
			puts "Adapted $newfile"
		} elseif {$file ne $newfile} {
			gzclose $f
			file_rename $file $newfile
		}
	} elseif {$file ne $newfile} {
		file_rename $file $newfile
	}
}

proc renamesamples {dir changes {relink 0}} {
	putslog "converting $dir"
	foreach file [glob -nocomplain $dir/*] {
		if {[file isdir $file] && [catch {file link $file} link]} {
			renamesamples $file $changes
		}
		renamesamples_file $file $changes $relink
	}
	set newdir [renamesamples_newfilename $dir $changes]
	if {$newdir ne $dir} {
		if {[file exists $newdir]} {error "cannot rename $dir to $newdir: already exists"}
		file_rename $dir $newdir
	}
}

proc cg_renamesamples {dir args} {
	if {[llength $args] == 1} {
		set changes [lindex $args 0]
	} else {
		set changes $args
	}
	if {[expr {[llength $changes]%2}]} {error "renamesamples requires an even number of arguments after the file/dirname (oldname newname)"}
	set dir [file_absolute $dir]
	if {![file isdir $dir]} {
		renamesamples_file $dir $changes
	} else {
		putslog "convert root $dir"
		cd [file dir $dir]
		set dir [file tail $dir]
		file delete -force $dir.temp
		hardlink $dir $dir.temp
		renamesamples $dir.temp $changes
		file delete -force $dir.old
		file rename $dir $dir.old
		file rename $dir.temp $dir
	}
}
