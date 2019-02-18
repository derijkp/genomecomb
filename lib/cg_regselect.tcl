proc regselect {regfile1 regfile2 {near -1} {o stdout}} {
	global cache
	# catch {close $f1}
	# catch {close $f2}
	set f1 [gzopen $regfile1]
	set poss1 [open_region $f1]
	gzclose $f1
	if {$regfile2 ne ""} {
		set f2 [gzopen $regfile2]
		set poss2 [open_region $f2]
		gzclose $f2
	} else {
		set poss2 {0 0 0}
	}
	# puts [list reg_select {*}$poss1 $regfile2 {*}$poss2 $near]
	set cat [gzcat $regfile1]
	exec {*}$cat $regfile1 | reg_select {*}$poss1 $regfile2 {*}$poss2 $near >@ $o 2>@ stderr
}

proc cg_regselect {args} {
	set o stdout
	set near -1
	cg_options reports_vars args {
		-o {
			set o [wgzopen $value]
		}
		-near {
			set near $value
		}
	} {region_file1 region_file2 near} 2 3
	if {$near eq ""} {set near -1}
	regselect $region_file1 $region_file2 $near $o
}
