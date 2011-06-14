proc cg_makepvt {args} {
	set pos 0
	set sumfields {}
	set sorted 1
	set tempfile [tempfile]
	foreach {key value} $args {
		switch -- $key {
			-sumfields {
				set sumfields $value
			}
			-sorted {
				set sorted $value
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	set fields {}
	if {([llength $args] < 2) || ([llength $args] > 3)} {
		errorformat makepvt
		exit 1
	}
	foreach {file resultfile fields} $args break
	if {$fields eq ""} {
		set fields [list_remove [cg select -h $file] chromosome begin end]
	}
	set ufields $fields
	lappend ufields {numbases=$end - $begin} {numlines=1}
	cg select -f $ufields $file $tempfile
	if {$sorted} {
		set tempfile2 [tempfile]
		cg select -s $fields $tempfile $tempfile2
		set tempfile $tempfile2
	}
	lappend sumfields numbases numlines
	cg groupby -sorted $sorted -sumfields $sumfields $fields $tempfile $resultfile
}
