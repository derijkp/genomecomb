proc regselect {regfile1 regfile2 {near -1} {output {}}} {
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
	# puts [list reg_select {*}$poss1 $regfile2 {*}$poss2 $near]
	if {$regfile1 eq "-"} {
		set header [tsv_open stdin]
		set poss1 [tsv_basicfields $header 3]
		if {$output eq ""} {
			set o [open [list | reg_select {*}$poss1 $regfile2 {*}$poss2 $near >@ stdout 2>@ stderr] w]
		} else {
			set o [wgzopen $output -1 {} [list | reg_select {*}$poss1 $regfile2 {*}$poss2 $near]]
		}
		puts $o [join $header \t]
		fcopy stdin $o
		close $o
	} else {
		set f1 [gzopen $regfile1]
		set poss1 [open_region $f1]
		gzclose $f1
		if {$output eq ""} {
			set o stdout
		} else {
			set o [wgzopen $output]
		}
		set cat [gzcat $regfile1]
		exec {*}$cat $regfile1 | reg_select {*}$poss1 $regfile2 {*}$poss2 $near >@ $o 2>@ stderr
	}
}

proc cg_regselect {args} {
	set output {}
	set near -1
	if {[llength $args] == 1} {
		set args [list - [lindex $args 0]]
	}
	cg_options regselect args {
		-o {
			set output $value
		}
		-near {
			set near $value
		}
	} {file region_select_file near} 2 3
	if {$near eq ""} {set near -1}
	regselect $file $region_select_file $near $output
}
