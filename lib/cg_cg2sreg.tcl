proc cg_cg2sreg {args} {
	global scriptname action
	set sorted 0
	set pos 0
	foreach {key value} $args {
		switch -- $key {
			-sorted {
				set sorted [true $value]
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	if {[llength $args] == 2} {
		foreach {file outfile} $args break
	} else {
		errorformat cg2sreg
	}
	putslog "Extract $outfile from $file"
	if {$sorted} {
		cg select -q {$varType != "no-call" && $varType != "no-ref"} -f "chromosome begin end" $file $outfile.temp
	} else {
		cg select -q {$varType != "no-call" && $varType != "no-ref"} -f "chromosome begin end" -s "chromosome begin end" $file $outfile.temp
	}
	cg regjoin $outfile.temp > $outfile.temp2
	file rename -force $outfile.temp2 $outfile
	file delete $outfile.temp
}
