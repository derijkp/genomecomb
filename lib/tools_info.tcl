proc analysisinfo_write {dep target args} {
	file mkdir [file dir $target]
	set dep [gzroot $dep]
	set target [gzroot $target]
	if {![llength $args]} {
		if {[file exists $dep.analysisinfo]} {
			file copy -force $dep.analysisinfo $target.analysisinfo
		}
		return
	}
	if {[file exists $dep.analysisinfo]} {
		set f [open $dep.analysisinfo]
		set fields [tsv_open $f]
		foreach {field value} $args {
			if {$field eq "sample"} {
				set fields [list sample {*}$fields]
			} else {
				lappend fields $field
			}
		}
		set o [open $target.analysisinfo w]
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
	} elseif {[info exists ::analysisinfo]} {
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
	set o [open $target.analysisinfo w]
	puts $o [join $fields \t]
	puts $o [join $values \t]
	close $o
}