proc regsubtract {regfile1 regfile2} {
	global cache
	# catch {close $f1}
	# catch {close $f2}
	set f1 [gzopen $regfile1]
	set poss1 [open_region $f1]
	set f2 [gzopen $regfile2]
	set poss2 [open_region $f2]
	gzclose $f1; gzclose $f2
	set temp1 [gztemp $regfile1]
	set temp2 [gztemp $regfile2]
	# puts [list reg_subtract $temp1 {*}$poss1 $temp2 {*}$poss2]
	exec reg_subtract $temp1 {*}$poss1 $temp2 {*}$poss2 >@ stdout 2>@ stderr
	gzrmtemp $temp1 ; gzrmtemp $temp2
}

proc cg_regsubtract {args} {
	if {[llength $args] != 2} {
		errorformat regsubtract
	}
	foreach {region_file1 region_file2} $args break
	regsubtract $region_file1 $region_file2
}
