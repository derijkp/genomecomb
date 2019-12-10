proc cg_cg2sreg {args} {
	global scriptname action
	set sorted 0
	cg_options cg2sreg args {
		-sorted {
			set sorted [true $value]
		}
	} {file outfile}
	putslog "Extract $outfile from $file"
	if {$sorted} {
		cg select -q {$varType != "no-call" && $varType != "no-ref"} -f "chromosome begin end" $file $outfile.temp
	} else {
		cg select -q {$varType != "no-call" && $varType != "no-ref"} -f "chromosome begin end" -s "chromosome begin end" $file $outfile.temp
	}
	cg regjoin $outfile.temp > $outfile.temp2
	file rename -force -- $outfile.temp2 $outfile
	file delete $outfile.temp
}
