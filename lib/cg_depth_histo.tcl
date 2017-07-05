proc depth_histo {bamfile regionfile {max 1000}} {
	global cache
	# catch {close $f1}
	# catch {close $f2}
	if {$regionfile ne ""} {
		set f2 [gzopen $regionfile]
		set poss2 [open_region $f2]
		gzclose $f2
	} else {
		set poss2 {0 0 0}
	}
	# puts [list depth_histo $regionfile {*}$poss2 $near]
	exec samtools depth $bamfile | depth_histo $regionfile {*}$poss2 $max >@ stdout 2>@ stderr
}

proc cg_depth_histo {args} {
	set max 1000
	cg_options depth_histo args {
		-max {set max $value}
	} {bamfile regionfile} 2 2 {
		make histogram of depth in a bamfile, seperating targeted by a regionfile
	}
	depth_histo $bamfile $regionfile $max
}
