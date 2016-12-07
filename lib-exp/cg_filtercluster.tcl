# filter out variants that are closer together than dist
proc cg_filtercluster {args} {
	if {[llength $args] != 3} {
		error "format is: cg filtercluster file dist'
	}
	foreach {file dist destfile} $args break

	set f [gzopen $file]
	set header [tsv_open $f]
	set poss [tsv_basicfields $header 3]
	set o [open $destfile w]
	puts $o [join $header \t]
	set pline [split [gets $f] \t]
	foreach {pchr ps pe} [list_sub $pline $poss] break
	set write 1
	while {![eof $f]} {
		set line [split [gets $f] \t]
		foreach {chr s e} [list_sub $line $poss] break
		if {($chr ne $pchr) || ($s > [expr {$pe+$dist}])} {
			if {$write} {puts $o [join $pline \t]}
			set write 1
		} else {
			set write 0
		}
		set pline $line
		set pchr $chr
		set ps $s
		set pe $e
	}
	if {$write} {puts $o [join $pline \t]}
	close $o
	gzclose $f

}
