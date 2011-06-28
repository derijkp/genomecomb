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
	set h [cg select -h $file]
	if {$fields eq ""} {
		set fields [list_remove $h chromosome begin end]
	}
	set ufields $fields
	if {[lsearch [list_find $h {chromosome begin end}] -1] != -1} {
		set region 1
		lappend ufields {numbases=$end - $begin} {numlines=1}
	} else {
		set region 0
	}
	cg select -f $ufields $file $tempfile
	if {$sorted} {
		set tempfile2 [tempfile]
		cg select -s $fields $tempfile $tempfile2
		set tempfile $tempfile2
	}
	if {$region} {
		lappend sumfields numbases numlines
	}
	cg groupby -sorted $sorted -sumfields $sumfields $fields $tempfile $resultfile
}
