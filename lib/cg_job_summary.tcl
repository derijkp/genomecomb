proc log_addtime {infoaVar sample timehours starttime endtime cores} {
	upvar $infoaVar infoa
	if {![info exists infoa(time,$sample)]} {
		set infoa(time,$sample) [expr {$timehours * $cores}]
	} else {
		set infoa(time,$sample) [expr {$infoa(time,$sample) + $timehours * $cores}]
	}
	if {![info exists infoa(starttime,$sample)]} {
		set infoa(starttime,$sample) $starttime
	} elseif {[time_gt $infoa(starttime,$sample) $starttime]} {
		set infoa(starttime,$sample) $starttime
	}
	if {![info exists infoa(endtime,$sample)]} {
		set infoa(endtime,$sample) $endtime
	} elseif {[time_gt $endtime $infoa(endtime,$sample)]} {
		set infoa(endtime,$sample) $endtime
	}
	if {![info exists infoa(longest,$sample)]} {
		set infoa(longest,$sample) $timehours
	} elseif {$timehours > $infoa(longest,$sample)} {
		set infoa(longest,$sample) $timehours
	}
}

proc format4f {value} {
	if {$value eq ""} {return ""}
	format %.4f $value
}

proc cg_job_summary args {
	# init needs to be here, because otherwise variables like cgpub(pid) will 
	# be reset on first use of functions in jobs.tcl
	set cleanup success
	set force 0
	set removeold 0
	set clustercores 144
	set listanalyses {
		longshot longshot-
		meth_nanopolish meth_nanopolish
		clair3 clair3-
	}
	cg_options job_summary args {
		-clustercores {
			set clustercores $value
		}
		-listanalyses {
			set listanalyses $value
		}
	} logfile
	set logfile [file_absolute $logfile]
	catch {close $f} ; set f [open $logfile]	
	set header [tsv_open $f]
	if {![inlist $header cores]} {
		puts stderr "incomplete because logfile has no cores field: $logfile"
		set cores {}
	}
	set incomplete 0
	unset -nocomplain infoa
	unset -nocomplain samplesa
	while 1 {
		if {[gets $f line] == -1} break
		foreach $header [split $line \t] break
		# putsvars {*}$header
		if {$job eq "total"} continue
		set target [lindex $targets 0]
		set starget [file split $target]
		set pos [lsearch $starget samples]
		if {$pos > -1} {
			incr pos
			set sample [lindex $starget $pos]
		} else {
			set sample {}
		}
		if {![info exists samplesa($sample)]} {
			set samplesa($sample) 1
			set infoa(incomplete,$sample) 0
			set infoa(analyses,$sample) {}
		}
		set analysis {}
		foreach {aname pattern} $listanalyses {
			if {[regexp $pattern $job]} {
				set analysis $aname
				if {![info exists infoa(a,$aname,$sample)]} {
					set infoa(a,$aname,$sample) 1
					lappend infoa(analyses,$sample) $aname
				}
				break
			}
		}
		if {$status eq "skipped"} {
			set incomplete 1
			incr infoa(incomplete,$sample)
		}
		if {$starttime eq "" || $endtime eq ""} {
			continue
		}
		# lappend infoa(target,$sample) $target
		if {![isint $cores]} {
			set cores 1
			set incomplete 1
			incr infoa(incomplete,$sample)
		}
		set timehours [timebetween_inhours $starttime $endtime]
		log_addtime infoa s,$sample $timehours $starttime $endtime $cores
		log_addtime infoa a,$analysis,$sample $timehours $starttime $endtime $cores
	}
	# parray infoa
	close $f
	#
	# output
	puts [join {type sample analyses incomplete coretime clustertime walltime starttime endtime longest} \t]
	set samples [lsort -dict [array names samplesa]]
	set totaltime 0.0
	set totalanalyses start
	set starttimes {}
	set endtimes {}
	set longests {}
	set totalmb 0
	foreach sample $samples {
		# get timing
		if {[get infoa(incomplete,$sample) 1] > 0} {set incomplete incompletetime} else {set incomplete {}}
		set time [format4f [get infoa(time,s,$sample) 0]]
		set clustertime [format4f [expr {[get infoa(time,s,$sample) 0]/$clustercores}]]
		set analyses $infoa(analyses,$sample)
		if {$totalanalyses eq "start"} {
			set totalanalyses $analyses
		} elseif {$analyses ne $totalanalyses} {
			list_addnew totalanalyses mixed {*}$analyses
		}
		set starttime [get infoa(starttime,s,$sample) {}]
		set endtime [get infoa(endtime,s,$sample) {}]
		lappend starttimes $starttime
		lappend endtimes $endtime
		set walltime [format4f [timebetween_inhours $starttime $endtime]]
		set longest [format4f [get infoa(longest,s,$sample) {}]]
		lappend longests $longest
		puts [join [list $type $sample $analyses $incomplete $time $clustertime $walltime $starttime $endtime $longest] \t]
		set totaltime [expr {$totaltime + [get infoa(time,s,$sample) 0]}]
		if {[info exists infoa(starttime,s,$sample)] && [time_gt $starttime $infoa(starttime,s,$sample)]} {
			set starttime $infoa(starttime,s,$sample)
		}
		if {[info exists infoa(endtime,s,$sample)] && [time_gt $infoa(endtime,s,$sample) $endtime]} {
			set endtime $infoa(endtime,s,$sample)
		}
		# analyses
		if {$sample ne ""} {
			foreach analysis $infoa(analyses,$sample) {
				if {![info exists infoa(time,a,$analysis,$sample)]} continue
				set type $analysis
				set time [format4f $infoa(time,a,$analysis,$sample)]
				set clustertime [format4f [expr {$time/$clustercores}]]
				set starttime $infoa(starttime,a,$analysis,$sample)
				set endtime $infoa(endtime,a,$analysis,$sample)
				set walltime [format4f [timebetween_inhours $starttime $endtime]]
				set longest [format4f $infoa(longest,a,$analysis,$sample)]
				puts [join [list $type $sample {} $incomplete $time $clustertime $walltime $starttime $endtime $longest] \t]
			}
		}
	}
	set type total
	set time [format4f $totaltime]
	set clustertime [format4f [expr {$totaltime/$clustercores}]]
	if {$incomplete > 0} {set incomplete incompletetime} else {set incomplete {}}
	set starttime [lindex [lsort -dict $starttimes] 0]
	set endtime [lindex [lsort -dict $endtimes] end]
	set longest [lindex [lsort -dict $longests] end]
	set walltime [format4f [timebetween_inhours $starttime $endtime]]
	set sample [llength $samples]
	if {[lsearch $totalanalyses mixed]} {set totalanalyses [list_union mixed $totalanalyses]}
	puts [join [list $type $sample $totalanalyses $incomplete $time $clustertime $walltime $starttime $endtime $longest] \t]
}
