proc cg_checktsv {file} {
	set f [gzopen $file]
	set header [tsv_open $f]
	set poss [tsv_basicfields $header 4 0]
	set line [split [gets $f] \t]
	# foreach {pchr pbegin pend ptype} [list_sub $line $poss] break
	set prev [list_sub $line $poss]
	set len [llength $header]
	set linenr 0
	while 1 {
		if {[llength $line] != $len} {
			puts "line $linenr is of wrong length: [llength $line] iso $len\t$line"
		}
		if {[eof $f]} break
		set line [split [gets $f] \t]
		set llen [llength $line]
		if {!$llen && [eof $f]} break
		incr linenr
		set cur [list_sub $line $poss]
		if {[list $prev $cur] ne [ssort -natural [list $prev $cur]]} {
			puts "line $linenr is sorted wrong:\t$line"
		}
		set prev $cur
	}
	close $f
}
