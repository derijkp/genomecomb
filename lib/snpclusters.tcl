#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg_clusters_distgraph {} {
	set f stdin
	set o stdout

	set poss [open_region $f]
	set num 0
	puts $o "begin\tdistance"
	foreach {pchr pbegin pend} [get_region $f $poss] break
	set pnchr [chr2num $pchr]
			while {![eof $f]} {
				incr num
				if {![expr $num%100000]} {putslog $num}
				set line [get_region $f $poss]
				foreach {pchr pbegin pend} $line break
				set pnchr [chr2num $pchr]
				if {$pnchr == 6} break
			}
			# seek $f 97574037
	while {![eof $f]} {
		incr num
		if {![expr $num%100000]} {putsprogress $num}
		set line [get_region $f $poss]
		if {![isint [lindex $line 2]]} continue
		foreach {chr begin end} $line break
		set nchr [chr2num $chr]
		if {$nchr == $pnchr} {
			set dist [expr {$begin - $pbegin}]
			puts $o $pbegin\t$dist
		} else {
			error "only one chromsome at the time"
		}
		foreach {pchr pnchr pbegin pend} [list $chr $nchr $begin $end] break
	}
	close $o
	close $f
}

proc clusters {} {
	set f stdin
	set o stdout
	set poss [open_region $f header]
	set fstart [tell $f]
	#
	puts $o "chromosome\tbegin\tend\tsize"
	set checknum 20
	set schecknum [expr {$checknum-1}]
	set maxdiff 150
	set breakdiff 10000
	#
	foreach {chr begin end} [get_region $f $poss] break
	foreach {pchr pbegin pend} [list $chr $begin $end] break
	set keep [list [list $pbegin 5000 0]]
	set checking 0
	set numtotal 0
	set numsmall 0
	set match 2
	set mismatch -5
	set maxp 0
	#
	set num 0
	while {![eof $f]} {
		set pos 1
		while {![eof $f]} {
			incr pos
			incr num
			if {![expr $num%100000]} {putsprogress "${num} ($chr:$begin)"}
			foreach {chr begin end} [get_region $f $poss] break
			if {$chr != $pchr} break
			set diff [expr {$begin-$pbegin}]
			if {$diff < $maxdiff} {incr numsmall}
			set percent [expr {double($numsmall)/$checknum}]
			lappend keep [list $begin $diff $percent]
			if {$pos >= $checknum} break
			foreach {pchr pbegin pend} [list $chr $begin $end] break
		}
		if {$pos < $checknum} continue
		set checking 0
		set score 0
		while {![eof $f]} {
			incr num
			if {![expr $num%100000]} {putsprogress "${num} ($chr:$begin)"}
			foreach {chr begin end} [get_region $f $poss] break
			if {$chr == $pchr} {
				set diff [expr {$begin-$pbegin}]
				if {$diff < $maxdiff} {incr numsmall}
				set rdiff [lindex $keep end-$schecknum 1]
				if {$rdiff < $maxdiff} {incr numsmall -1}
				set percent [expr {double($numsmall)/$checknum}]
				lappend keep [list $begin $diff $percent $score]
			} else {
				set diff 10000000
				set percent 0
				set score -10
			}
			# puts $o [join [list $begin $diff $percent $score] \t]
			if {$checking} {
				incr rend
				if {$diff < $maxdiff} {incr score $match} else {incr score $mismatch}
				if {$score >= $maxscore} {
					set maxscore $score
					set maxpos $rend
					set limit [expr {round($maxscore/4)}] 
				}
				if {($score <= $limit) || ($diff > $breakdiff)} {
					set rend $maxpos
					set size [expr {$rend-$rstart+1}]
					if {$size > 20} {
						# puts $o -----end-----
						if 0 {
							puts $o \n\n
							puts $o "[llength $keep]\t$chr:[lindex $keep 0 0]-[lindex $keep end 0] ($num)"
							puts $o ----------
							puts $o [join [lrange $keep 0 [expr {$rstart-1}]] \n]
							puts $o *****start*****
							puts $o [join [lrange $keep $rstart $rend] \n]
							puts $o *****end*****
							puts $o [join [lrange $keep [expr {$rend+1}] end] \n]
						}
						set begin [lindex $keep $rstart 0]
						set end [expr {[lindex $keep $rend 0]+1}]
						puts $o $chr\t$begin\t$end\t[expr {$end - $begin + 1}]
						# putslog $chr\t$begin\t$end\t[expr {$end - $begin + 1}]
					}
					set checking 0
					set keep [lrange $keep end-$checknum end]
					set score 0
				}
			} elseif {$percent >= 0.9} {
				set score 0
				set limit [expr {(round($checknum*(1-$percent))+2)*+$mismatch}]
				set pos [expr {[llength $keep]-1}]
				set maxscore $score
				set maxpos $pos
				while {$pos >= 0} {
					set diff [lindex $keep $pos 1]
					if {$diff > $breakdiff} break
					if {$diff < $maxdiff} {incr score $match} else {incr score $mismatch}
					if {$score >= $maxscore} {
						set maxscore $score
						set maxpos $pos
						if {$score > 20} {
							set limit [expr {round($maxscore/4)}]
						}
					}
					if {$score <= $limit} break
					incr pos -1
				}
				set rstart $maxpos
				set checking 1
				set rend [expr {[llength $keep]-1}]
				set score [expr {round(10*$percent)}]
				set limit 0
				set maxscore $score
				# puts $o -----start-----
			} elseif {$percent < 0.2} {
				set keep [lrange $keep end-$checknum end]
			}
			if {$chr != $pchr} {break}
			foreach {pchr pbegin pend} [list $chr $begin $end] break
		}
		flush $o
		foreach {pchr pbegin pend} [list $chr $begin $end] break
		set keep [list [list $begin 10000 0]]
		set checking 0
		set numtotal 0
		set numsmall 0
		set maxp 0
	}
}

proc cg_clusterregions {args} {
	global scriptname action
	if {[llength $args] != 0} {
		puts stderr "format is: $scriptname $action"
		puts stderr " - outputs regions with clusters of variations"
		puts stderr " - input is an annotated variations file"
		exit 1
	}
	clusters
}
