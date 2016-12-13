proc cg_makepvt {args} {
	set sumfields {}
	set sorted 1
	set tempfile [tempfile]
	set fields {}
	cg_options makepvt args {
		-sumfields {
			set sumfields $value
		}
		-sorted {
			set sorted $value
		}
	} {file resultfile fields} 2 3
	set h [cg select -h $file]
	if {$fields eq ""} {
		set fields [list_remove $h chromosome begin end]
	}
	set ufields $fields
	if {[lsearch [tsv_basicfields $h 3] -1] == -1} {
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
