proc cg_bam_histo {args} {
	set pos 0
	set intervals {1 5 10 20 50 100 200 500 1000}
	cg_options bam_histo args {
		-n - -namecol {
			set namecol $value
		}
	} {regionfile bamfile intervals} 2 3
	set chrs [bam_chrs $bamfile]
	if {[regexp ^chr [lindex $chrs 0]]} {set pre chr} else {set pre {}}
	foreach bamchr $chrs {
		set bamchrsa($bamchr) 1
	}
	set f [gzopen $regionfile]
	set header [tsv_open $f]
	set poss [tsv_basicfields $header 3]
	if {[info exists namecol]} {
		set namepos [lsearch $header $namecol]
	} else {
		foreach namecol {name info} {
			set namepos [lsearch $header $namecol]
			if {$namepos != -1} break
		}
	}
	lappend poss $namepos
	set biv <[lindex $intervals 0]
	set iv $biv
	set header [list name]
	set p {}
	set tota($iv) 0
	foreach limit $intervals {
		set tota($limit) 0
	}
	set a($iv) 0
	foreach limit $intervals {
		set a($limit) 0
		lappend header r$p<$limit
		set p $limit
	}
	lappend header r${limit}<
	set totsize 0
	set totsum 0
	set totmin {}
	set totmax {}
	set size 0
	set sum 0
	set mins {}
	set maxs {}
	lappend header size avg min max
	puts [join $header \t]
	set line [getline $f]
	set prevname [lindex [list_sub $line $poss] 3]
	while {1} {
		foreach {chr begin end name} [list_sub $line $poss] break
		if {[eof $f] || $name ne $prevname} {
			incr tota($biv) $a($biv)
			set result [list $a($biv)]
			set iv $biv
			set a($iv) 0
			foreach limit $intervals {
				incr tota($limit) $a($limit)
				lappend result $a($limit)
				set a($limit) 0
			}
			set min [lmath_min $mins]
			set max [lmath_max $maxs]
			lappend result $size
			if {$size == 0} {
				lappend result ""
			} else {
				lappend result [format %.2f [expr {$sum/double($size)}]]
			}
			lappend result $min
			lappend result $max
			puts $prevname\t[join $result \t]
			set prevname $name
			incr totsize $size
			incr totsum $sum
			if {$totmin eq "" || $min < $totmin} {set totmin $min}
			if {$totmax eq "" || $max > $totmax} {set totmax $max}
			set size 0
			set sum 0
			set mins {}
			set maxs {}
			if {[eof $f]} break
		}
		incr end -1
		set chr [chr_clip $chr]
		putslog $chr:$begin-$end
		if {![info exists bamchrsa($pre$chr)]} {
			if {$chr eq "M"} {set chr MT} elseif {$chr eq "MT"} {set chr M}
		}
		set temp [exec samtools depth -d1000000 -q 1 -r $pre$chr:[expr {$begin+1}]-[expr {$end+1}] $bamfile]
		set data [split [string trim $temp] \n]
		set data [list_subindex $data 2]
		set explen [expr {$end-$begin+1}]
		set missing [expr {$explen-[llength $data]}]
		if {$missing} {
			lappend data {*}[list_fill $missing 0]
		}
		foreach v $data {
			set iv $biv
			foreach limit $intervals {
				if {$v < $limit} break
				set iv $limit
			}
			incr a($iv)
		}
		incr size [llength $data]
		set sum [expr {$sum + round([lmath_sum $data])}]
		lappend mins [lmath_min $data]
		lappend maxs [lmath_max $data]
		set line [getline $f]
	}
	gzclose $f
	puts ----------
	set result [list $tota($biv)]
	foreach limit $intervals {
		lappend result $tota($limit)
	}
	set tot [lmath_sum $result]
	if {$tot > 0} {
		set presult [list [format %.2f [expr {100*$tota($biv)/$tot}]]]
		foreach limit $intervals {
			lappend presult [format %.2f [expr {100*$tota($limit)/$tot}]]
		}
		lappend result $totsize [format %.2f [expr {$totsum/double($totsize)}]] $totmin $totmax
	} else {
		set presult [list [format %.2f 0]]
		foreach limit $intervals {
			lappend presult [format %.2f 0]
		}
		lappend result $totsize [format %.2f 0] $totmin $totmax
	}
	puts Total\t[join $result \t]
	puts Totalpercent\t[join $presult \t]
}
