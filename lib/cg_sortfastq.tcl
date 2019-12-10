proc cg_sortfastq {args} {
	cg_options sortfastq args {
	} {infile outfile} 0 2
	if {![info exists infile]} {
		exec paste - - - - | gnusort8 -T [scratchdir] -t \t -V -s -k 1,1 | tr {\t} {\n} <@ stdin >@ stdout
	} elseif {![info exists outfile]} {
		exec {*}[gzcat $infile] $infile | paste - - - - | gnusort8 -T [scratchdir] -t \t -V -s -k 1,1 | tr {\t} {\n} >@ stdout
	} else {
		exec {*}[gzcat $infile] $infile | paste - - - - | gnusort8 -T [scratchdir] -t \t -V -s -k 1,1 | tr {\t} {\n} {*}[compresspipe $outfile] > $outfile.temp
		file rename -force -- $outfile.temp $outfile
	}
}
