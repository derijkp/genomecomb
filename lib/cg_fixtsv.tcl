proc cg_fixtsv {file outfile} {
	set f [gzopen $file]
	set header [tsv_open $f comment]
	# set poss [tsv_basicfields $header 4 0]
	set o [open $outfile.temp w]
	puts -nonewline $o $comment
	puts $o [join $header \t]
	set tline [gets $f]
	set line [split $tline \t]
	# foreach {pchr pbegin pend ptype} [list_sub $line $poss] break
	# set prev [list_sub $line $poss]
	set len [llength $header]
	set llen [llength $line]
	set linenr 0
	while 1 {
		if {$llen != $len} {
			puts "line $linenr is of wrong length: [llength $line] iso $len\t$line, fixing"
			if {$llen < $len} {
				lappend line {*}[list_fill [expr {$len-$llen}] {}]
			}
			puts $o [join $line \t]
		} else {
			puts $o $tline
		}
		if {[eof $f]} break
		set tline [gets $f]
		set line [split $tline \t]
		set llen [llength $line]
		if {!$llen && [eof $f]} break
		incr linenr
	}
	close $f
	close $o
	if {[file exists $outfile]} {catch {file rename $outfile $outfile.old}}
	file rename -force $outfile.temp $outfile
}
