proc cg_bam_histo {args} {
	set pos 0
	set namecol name
	foreach {key value} $args {
		switch -- $key {
			-n - --namecol {
				set namecol $value
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	if {[llength $args] != 3} {
		exiterror "wrong # args: should be \"cg bam_histo regionfile bamfile intervals\""
	}
	foreach {regionfile bamfile intervals} $args break
	set f [gzopen $regionfile]
	set header [tsv_open $f]
	set poss [tsv_basicfields $header 3]
	set namepos [lsearch $header $namecol]
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
	set prevname [lindex $line 3]
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
		putslog $chr:$begin-$end
		set data [split [string trim [exec samtools depth -q 1 -r $chr:[expr {$begin+1}]-[expr {$end+1}] $bamfile]] \n]
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
	set presult [list [format %.2f [expr {100*$tota($biv)/$tot}]]]
	foreach limit $intervals {
		lappend presult [format %.2f [expr {100*$tota($limit)/$tot}]]
	}
	lappend result $totsize [format %.2f [expr {$totsum/double($totsize)}]] $totmin $totmax
	puts Total\t[join $result \t]
	puts Totalpercent\t[join $presult \t]
}
