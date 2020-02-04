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

proc result_rename {src target} {
	file rename -force -- $src $target
	catch {file rename -force -- [analysisinfo_file $src] [analysisinfo_file $target]}
}
