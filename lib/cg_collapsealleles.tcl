proc cg_collapsealleles {args} {
	if {([llength $args] > 2)} {
		errorformat collapsealleles
		exit 1
	}
	foreach {file outfile} {{} {}} break
	foreach {file outfile} $args break
	if {$file eq ""} {
		set f stdin
	} else {
		set f [gzopen $file]
	}
	if {$outfile eq ""} {
		set o stdout
	} else {
		set o [open $outfile w]
	}
	set header [tsv_open $f comment]
	set poss [tsv_basicfields $header]
	set apos [lindex $poss 5]
	set rpos [lindex $poss 4]
	set poss [lrange $poss 0 4]
	set samples [samples $header]
	if {[llength $samples]} {
		set sposs {}
		foreach sample $samples {
			lappend sposs {*}[list_cor $header [list sequenced-$sample zyg-$sample]]
		}
	} else {
		set sposs [list_cor $header [list sequenced zyg]]
	}
	if {[string length $comment]} {
		puts [string trim $comment]
	}
	puts [join $header \t]
	set prevline [split [gets $f] \t]
	set prevloc [list_sub $prevline $poss]
	set cur [list $prevline]
	while {1} {
		set line [split [gets $f] \t]
		set loc [list_sub $line $poss]
		if {$prevloc ne $loc} {
			if {[llength $cur] > 1} {
				set len [llength [lindex $cur 0]]
				set resultline {}
				for {set i 0} {$i < $len} {incr i} {
					set v [list_subindex $cur $i]
					set temp [list_remove $v ?]
					if {[llength $temp]} {set v $temp}
					set cv [list_remdup $v]
					if {[llength $cv] > 1} {
						lappend resultline [join $v ,] 
					} else {
						lappend resultline [lindex $cv 0]
					}
				}
				set ref [lindex $resultline $rpos]
				foreach {seqpos zygpos} $sposs {
					if {$seqpos != -1} {
						set seq [lindex $resultline $seqpos]
						if {[inlist $seq v]} {set seq v} elseif {[inlist $seq u]} {set seq u} else {set seq r}
						lset resultline $seqpos $seq
					}
					if {$zygpos != -1} {
						set zyg [list_remdup [lindex $resultline $zygpos]]
						if {[inlist $zyg u]} {
							set zyg u
						} elseif {[inlist $zyg m]} {
							set zyg m
						} elseif {[inlist $zyg t]} {
							set zyg t
						} elseif {[inlist $zyg c]} {
							set zyg c
						} elseif {[inlist $zyg o]} {
							set zyg c
						} else {
							set zyg r
						}
						lset resultline $zygpos $zyg
					}
				}
				puts [join $resultline \t]
			} else {
				puts [join [lindex $cur 0] \t]
			}
			set cur {}
			set prevloc $loc
		}
		if {![llength $line] && [eof $f]} break
		lappend cur $line
	}
	close $f
}
