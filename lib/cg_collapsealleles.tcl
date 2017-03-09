proc cg_collapsealleles {args} {
	foreach {file outfile} {{} {}} break
	set duplicates k
	array set duplicatesa {merge m keep k first f error e}
	cg_options collapsealleles args {
		-duplicates {
			if {[lindex $value 0] eq "max"} {
				set duplicates +
				set duplicatesmax [lindex $value 1]
			} elseif {[lindex $value 0] eq "min"} {
				set duplicates -
				set duplicatesmax [lindex $value 1]
			} elseif {![info exists duplicatesa($value)]} {
				error "wrong value for -duplicates: $value"
			} else {
				set duplicates $duplicatesa($value)
			}
		}
	} {file outfile} 0 2 {
		Convert a split variant file to unsplit by collapsing alleles of the same variant into one line.
	}
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
	set poss [tsv_basicfields $header 6 0]
	set apos [lindex $poss 5]
	set rpos [lindex $poss 4]
	set poss [lrange $poss 0 4]
	if {[info exists duplicatesmax]} {set mpos [lsearch $header $duplicatesmax]}
	unset -nocomplain typea
	set pos 0
	foreach field $header {
		if {[regexp {^sequenced(-|$)} $field]} {
			set typea($pos) s
		} elseif {[regexp {^zyg(-|$)} $field]} {
			set typea($pos) z
		} elseif {[regexp {^genotypes(-|$)} $field]} {
			set typea($pos) g
		} else {
			set typea($pos) n
		}
		incr pos
	}
	if {$apos != -1 && $duplicates eq "k"} {
		set typea($apos) a
	}
	if {[string length $comment]} {
		regsub "#split\t1" $comment "#split\t0" comment
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
			set alts [list_subindex $cur $apos]
			if {[llength $alts] > 1} {
				set altsd [list_remdup $alts]
				if {[llength $altsd] < [llength $alts]} {
					if {$duplicates eq "f"} {
						set cur [list_sub $cur [list_cor $alts $altsd]]
					} elseif {$duplicates eq "+"} {
						set list [list_subindex $cur $mpos]
						set slist [lsort -dict -decreasing $list]
						set cor [list_cor $list $slist]
						set cur [list_sub $cur $cor]
						set alts [list_sub $alts $cor]
						set cur [list_sub $cur [list_cor $alts $altsd]]
					} elseif {$duplicates eq "-"} {
						set list [list_subindex $cur $mpos]
						set slist [lsort -dict $list]
						set cor [list_cor $list $slist]
						set cur [list_sub $cur $cor]
						set alts [list_sub $alts $cor]
						set cur [list_sub $cur [list_cor $alts $altsd]]
					} elseif {$duplicates eq "e"} {
						error "duplicates found:\n[join $cur \n]"
					}
				}
			}
			if {[llength $cur] == 1} {
				puts [join [lindex $cur 0] \t]
			} else {
				set len [llength [lindex $cur 0]]
				set resultline {}
				for {set i 0} {$i < $len} {incr i} {
					set list [list_subindex $cur $i]
					set type $typea($i)
					if {$type eq "a"} {
						set list [join $list ,] 
					} elseif {$type eq "n"} {
						set temp [list_remove $list ?]
						if {[llength $temp]} {set list $temp}
						set clist [list_remdup $list]
						if {[llength $clist] > 1} {
							set list [join $list ,] 
						} else {
							set list [lindex $clist 0]
						}
					} elseif {$type eq "s"} {
						if {[inlist $seq v]} {set list v} elseif {[inlist $seq u]} {set list u} else {set list r}
					} elseif {$type eq "z"} {
						if {[inlist $list u]} {
							set list u
						} elseif {[inlist $list m]} {
							set list m
						} elseif {[inlist $list t]} {
							set list t
						} elseif {[inlist $list c]} {
							set list c
						} elseif {[inlist $list o]} {
							set list c
						} else {
							set list r
						}
					} elseif {$type eq "g"} {
						set first [lindex $list 0]
						set connect ";"
						regexp {[,;]} $first connect
						set first [split $first ",;"]
						set glen [llength $first]
						set genotypes [list_fill $glen -1]
						set apos 1
						foreach pos [list_find $first 0] {
							lset genotypes $pos 0
						}
						foreach g $list {
							set g [split $g ",;"]
							foreach pos [list_find $g 1] {
								lset genotypes $pos $apos
							}
							incr apos
						}
						set list [join $genotypes $connect]
					}
					lappend resultline $list
				}
				puts [join $resultline \t]
			}
			set cur {}
			set prevloc $loc
		}
		if {![llength $line] && [eof $f]} break
		lappend cur $line
	}
	if {$f ne "stdin"} {gzclose $f}
	if {$o ne "stdout"} {close $o}
}
