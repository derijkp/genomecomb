proc cg_fixtsv {args} {
	cg_options fixtsv args {
	} {file outfile} 2 2
	set f [gzopen $file]
	set header [tsv_open $f comment]
	# set poss [tsv_basicfields $header 4 0]
	set o [open $outfile.temp w]
	puts -nonewline $o $comment
	set rheader [list_remdup $header]
	if {[llength $rheader] < [llength $header]} {
		set cor [list_cor $header $rheader]
		set duplicates [list_sub $header -exclude $cor] 
		puts "header shows duplicate fields ([join $duplicates ", "]), removing"
		puts $o [join [list_sub $header $cor] \t]
	} else {
		set cor {}
		puts $o [join $header \t]
	}
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
			} else {
				set line [lrange $line 0 [expr {$len-1}]]
			}
			if {[llength $cor]} {set line [list_sub $line $cor]}
			puts $o [join $line \t]
		} else {
			if {[llength $cor]} {
				set line [list_sub $line $cor]
				puts $o [join $line \t]
			} else {
				puts $o $tline
			}
		}
		if {[eof $f]} break
		set tline [gets $f]
		set line [split $tline \t]
		set llen [llength $line]
		if {!$llen && [eof $f]} break
		incr linenr
	}
	gzclose $f
	close $o
	if {[file exists $outfile]} {catch {file rename -- $outfile $outfile.old}}
	file rename -force -- $outfile.temp $outfile
}
