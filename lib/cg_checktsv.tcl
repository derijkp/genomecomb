proc cg_checktsv {args} {
	set checksort 1
	cg_options checktsv args {
		-checksort {set checksort $value}
	} {file}
	set f [gzopen $file]
	set header [tsv_open $f]
	set poss [tsv_basicfields $header 4 0]
	set line [split [gets $f] \t]
	# foreach {pchr pbegin pend ptype} [list_sub $line $poss] break
	set rheader [list_remdup $header]
	if {[llength $rheader] < [llength $header]} {
		puts "header has duplicate fields: [list_sub $header -exclude [list_cor $header $rheader]]"
	}
	set prev [list_sub $line $poss]
	set len [llength $header]
	set linenr 0
	set error 0
	while 1 {
		if {[llength $line] != $len} {
			puts stderr "line $linenr is of wrong length: [llength $line] iso $len\t$line"
			set error 1
		}
		if {[eof $f]} break
		set line [split [gets $f] \t]
		set llen [llength $line]
		if {!$llen && [eof $f]} break
		incr linenr
		if {$checksort} {
			set cur [list_sub $line $poss]
			if {[list $prev $cur] ne [bsort -sortchromosome [list $prev $cur]]} {
				puts stderr "line $linenr is sorted wrong:\t$line (prev = $prev)"
				set error 1
			}
			set prev $cur
		}
	}
	close $f
	if {$error} {error "checking tsv file $file concluded with errors"}
}
