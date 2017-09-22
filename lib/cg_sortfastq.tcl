proc cg_sortfastq {args} {
	cg_options sortfastq args {
	} {infile outfile} 0 2
	if {![info exists infile]} {
		exec paste - - - - | gnusort8 -T [scratchdir] -t \t -V -s -k 1,1 | tr {\t} {\n} <@ stdin >@ stdout
	} elseif {![info exists outfile]} {
		exec {*}[gzcat $infile] $infile | paste - - - - | gnusort8 -T [scratchdir] -t \t -V -s -k 1,1 | tr {\t} {\n} >@ stdout
	} else {
		exec {*}[gzcat $infile] $infile | paste - - - - | gnusort8 -T [scratchdir] -t \t -V -s -k 1,1 | tr {\t} {\n} {*}[compresspipe $outfile] > $outfile.temp
		file rename -force $outfile.temp $outfile
	}
}

proc cg_fastq2tsv {args} {
	cg_options fastq2tsv args {
	} {infile outfile} 0 2
	set header [join {id sequence temp quality} \t]
	if {[info exists outfile]} {
		set compresspipe [compresspipe $outfile]
		if {$compresspipe ne ""} {
			set o [open [list {*}$compresspipe > $outfile.temp] w]
		} else {
			set o [open $outfile.temp w]
		}
	} else {
		set o stdout
	}
	if {![info exists infile]} {
		puts $header
		exec paste - - - - <@ stdin >@ stdout
	} elseif {![info exists outfile]} {
		puts $header
		exec {*}[gzcat $infile] $infile | paste - - - - >@ stdout
	} else {
		puts $o $header
		exec {*}[gzcat $infile] $infile | paste - - - - >@ $o
		close $o
		file rename -force $outfile.temp $outfile
	}
}