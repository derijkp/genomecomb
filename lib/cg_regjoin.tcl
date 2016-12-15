proc regjoin {regfile1 regfile2} {
	global cache
	# catch {close $f1}
	# catch {close $f2}
	if {$regfile2 ne ""} {
		set f2 [gzopen $regfile2]
		set poss2 [open_region $f2]
		gzclose $f2
	} else {
		set poss2 {0 0 0}
	}
	if {$regfile1 ne ""} {
		set f1 [gzopen $regfile1]
		set poss1 [open_region $f1]
		gzclose $f1
		set temp1 [gztemp $regfile1]
		set temp2 [gztemp $regfile2]
		# puts [list reg_join $regfile1 {*}$poss1 $regfile2 {*}$poss2]
		exec reg_join $temp1 {*}$poss1 $temp2 {*}$poss2 >@ stdout 2>@ stderr
#		gzrmtemp $temp1 $temp2
	} else {
		set poss1 [open_region stdin]
		set temp2 [gztemp $regfile2]
		set o [open "| [list reg_join - {*}$poss1 $temp2 {*}$poss2] >@ stdout 2>@ stderr" w]
		fcopy stdin $o
		close $o
		gzrmtemp $temp2
	}
}

proc cg_regjoin {args} {
	if {[llength $args] > 2} {
		errorformat regjoin
	}
	foreach {region_file1 region_file2} {{} {}} break
	foreach {region_file1 region_file2} $args break
	regjoin $region_file1 $region_file2
}
