proc regjoin {regfile1 regfile2} {
	global cache
	# catch {close $f1}
	# catch {close $f2}
	set f1 [open $regfile1]
	set poss1 [open_region $f1]
	close $f1
	if {$regfile2 ne ""} {
		set f2 [open $regfile2]
		set poss2 [open_region $f2]
		close $f2
	} else {
		set poss2 {0 0 0}
	}
	# puts [list reg_join $regfile1 {*}$poss1 $regfile2 {*}$poss2]
	exec reg_join $regfile1 {*}$poss1 $regfile2 {*}$poss2 >@ stdout 2>@ stderr
}

proc cg_regjoin {args} {
	if {([llength $args] != 1) && ([llength $args] != 2)} {
		errorformat regjoin
		exit 1
	}
	foreach {region_file1 region_file2} $args break
	regjoin $region_file1 $region_file2
}
