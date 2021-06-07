proc analysisinfo_file {file} {
	return [gzroot $file].analysisinfo
}

proc analysisinfo_write {dep target args} {
	global env
	if {[ispipe $dep] || [ispipe $target]} {
		return
	}
	file mkdir [file dir $target]
	set dep [gzroot $dep]
	set target [gzroot $target]
	set depanalysisinfo [analysisinfo_file $dep]
	set targetanalysisinfo [analysisinfo_file $target]
	if {![llength $args]} {
		if {[file exists $depanalysisinfo] && $dep ne $target} {
			file_copy $depanalysisinfo $targetanalysisinfo
		}
		return
	}
	if {[file exists $depanalysisinfo] && [file size $depanalysisinfo] != 0} {
		if {$dep eq $target} {
			file rename -force -- $depanalysisinfo $depanalysisinfo.old
			set f [open $depanalysisinfo.old]
		} else {
			set f [open $depanalysisinfo]
		}
		set fields [tsv_open $f]
		set addvalues {}
		unset -nocomplain samplevalue
		foreach {field value} $args {
			if {$field in $fields} continue
			if {$field eq "sample"} {
				set fields [list sample {*}$fields]
				set samplevalue $value
			} else {
				lappend fields $field
				lappend addvalues $value
			}
		}
		set addvalues [join $addvalues \t]
		if {[llength $addvalues]} {set addvalues \t$addvalues}
		set o [open $targetanalysisinfo w]
		puts $o [join $fields \t]
		while {[gets $f line] != -1} {
			if {[info exists samplevalue]} {
				set line $samplevalue\t$line$addvalues
			} else {
				append line $addvalues
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
		file_copy $src $dest
	} else {
		set c [string trim [file_read $src]]
		set c [split $c \n]
		if {[llength $c] > 2} {
			error "can only change 1 line analysisinfo file: $src contains [expr {[llength $c]-1}]"
		}
		set header [split [lindex $c 0] \t]
		set data [split [lindex $c 1] \t]
		if {[llength $data] < [llength $header]} {
			lappend data {*}[list_fill [expr {[llength $header] - [llength $data]}] {}]
		}
		foreach {key value} $changes {
			set pos [lsearch $header $key]
			if {$pos == -1} continue
			lset data $pos $value
		}
		file_write $dest [join $header \t]\n[join $data \t]\n
	}
}

proc result_rename {src target} {
	if {[ispipe $src]} return
	file rename -force -- $src $target
	catch {file rename -force -- [analysisinfo_file $src] [analysisinfo_file $target]}
	set indexext [indexext $target]
	if {[file exists $src.$indexext]} {file rename -force -- $src.$indexext $target.$indexext}
	if {[file exists $src.tbi]} {file rename -force -- $src.tbi $target.tbi}
	if {[file exists $src.zsti]} {file rename -force -- $src.zsti $target.zsti}
}
