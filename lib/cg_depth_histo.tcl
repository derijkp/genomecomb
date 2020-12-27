proc depth_histo {bamfile regionfile {max 1000} {Q 0} {q 0} args} {
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
	if {$max > 1000000} {
		set d $max
	} else {
		set d 1000000
	}
	# puts [list samtools depth -d$d -Q $Q -q $q $bamfile | depth_histo $regionfile {*}$poss2 $max]
	exec samtools depth -d$d -Q $Q -q $q $bamfile | depth_histo $regionfile {*}$poss2 $max >@ stdout 2>@ stderr
}

proc cg_depth_histo {args} {
	set regionfile {}
	set max 1000
	set opts {}
	set Q 0
	set q 0
	cg_options depth_histo args {
		-max {set max $value}
		-q {set q $value}
		-Q {set Q $value}
	} {bamfile regionfile} 1 2 {
		make histogram of depth in a bamfile, seperating targeted by a regionfile
	}
	depth_histo $bamfile $regionfile $max $Q $q
}
