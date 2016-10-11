proc regsubtract {regfile1 regfile2} {
	global cache
	# catch {close $f1}
	# catch {close $f2}
	set f1 [open $regfile1]
	set poss1 [open_region $f1]
	set f2 [open $regfile2]
	set poss2 [open_region $f2]
	close $f1; close $f2
	# puts [list reg_subtract $regfile1 {*}$poss1 $regfile2 {*}$poss2]
	exec reg_subtract $regfile1 {*}$poss1 $regfile2 {*}$poss2 >@ stdout 2>@ stderr
}

proc cg_regsubtract {args} {
	if {[llength $args] != 2} {
		errorformat regsubtract
		exit 1
	}
	foreach {region_file1 region_file2} $args break
	regsubtract $region_file1 $region_file2
}
