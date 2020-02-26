proc analysisinfo_file {file} {
	return [gzroot $file].analysisinfo
}

proc analysisinfo_write {dep target args} {
	global env
	if {[file_root $dep] eq "-" || [file_root $target] eq "-"} {
		return
	}
	file mkdir [file dir $target]
	set dep [gzroot $dep]
	set target [gzroot $target]
	set depanalysisinfo [analysisinfo_file $dep]
	set targetanalysisinfo [analysisinfo_file $target]
	if {![llength $args]} {
		if {[file exists $depanalysisinfo] && $dep ne $target} {
			file copy -force $depanalysisinfo $targetanalysisinfo
		}
		return
	}
	if {[file exists $depanalysisinfo]} {
		if {$dep eq $target} {
			file rename -force -- $depanalysisinfo $depanalysisinfo.old
			set f [open $depanalysisinfo.old]
		} else {
			set f [open $depanalysisinfo]
		}
		set fields [tsv_open $f]
		foreach {field value} $args {
			if {$field eq "sample"} {
				set fields [list sample {*}$fields]
			} else {
				lappend fields $field
			}
		}
		set o [open $targetanalysisinfo w]
		puts $o [join $fields \t]
		while {[gets $f line] != -1} {
			set values [split $line \t]
			foreach {field value} $args {
				if {$field eq "sample"} {
					set line $value\t$line
				} else {
					append line \t $value
				}
			}
			puts $o $line
			
		}
		close $o
		close $f
		return
	} elseif {[info exists ::analysisinfo] && ![expr {[llength $::analysisinfo]%2}]} {
		set fields [list_unmerge $::analysisinfo 1 values]
	} else {
		set fields {}
		set values {}
	}
	foreach {field value} $args {
		if {$field eq "sample"} {
			set fields [list sample {*}$fields]
			set values [list $value {*}$values]
		} else {
			lappend fields $field
			lappend values $value
		}
	}
	set o [open $targetanalysisinfo w]
	puts $o [join $fields \t]
	puts $o [join $values \t]
	close $o
}

proc analysisinfo_copy {src dest {changes {}}} {
	if {![llength $changes]} {
		file copy $src $dest
	} else {
		set c [string trim [file_read $src]]
		set c [split $c \n]
		if {[llength $c] > 2} {
			error "can only change 1 line analysisinfo file: $src contains [expr {[llength $c]-1}]"
		}
		set header [split [lindex $c 0] \t]
		set data [split [lindex $c 1] \t]
		foreach {key value} $changes {
			set pos [lsearch $header $key]
			if {$pos == -1} {
				error "can not change analysisinfo file $src: field $key not present"
			}
			lset data $pos $value
		}
		file_write $dest [join $header \t]\n[join $data \t]\n
	}
}

proc result_rename {src target} {
	file rename -force -- $src $target
	catch {file rename -force -- [analysisinfo_file $src] [analysisinfo_file $target]}
	catch {file rename -force -- $src.[indexext $dest] $dest.[indexext $dest]}
}
