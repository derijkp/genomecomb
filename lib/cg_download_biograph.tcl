proc cg_download_biograph {file diseaseid {genefile {}}} {
	unset -nocomplain a
	if {$genefile ne ""} {
		set f [gzopen $genefile]
		set header [tsv_open $f]
		set temp [split [string trim [read $f]] \n]
		gzclose $f
		set genes [list_subindex $temp 0]
		foreach gene $genes {
			set a($gene) 1
		}
	} 
	exec wget --quiet -O $file.temp "http://biograph.be/concept/tsv/$diseaseid?filter_directness=Known+and+inferred&filter_type=Gene"
	# set count [lindex [exec wc -l $file.temp] 0]
	set f [open $file.temp]
	set o [open $file.temp2 w]
	puts $o [join {gene	rank score} \t]
	set currank 0
	while 1 {
		if {[gets $f line] == -1} break
		set line [split $line \t]
		foreach {temp id gene type known score} $line break
		set gene [lindex $gene 0]
		if {[info exists a($gene)]} {
			set rank 0
			unset a($gene)
		} elseif {$known eq "Known"} {
			set rank 1
		} else {
			incr currank
			set rank $currank
		}
		puts $o [join [list $gene $rank [format %.4g $score]] \t]
	}
	foreach gene [array names a] {
		puts $o [join [list $gene 0 -] \t]
	}
	close $o
	close $f
	if {[file extension $file] eq ",lz4"} {
		cg_lz4 $file.temp2
		file rename -force $file.temp2.lz4 $file
	} else {
		file rename -force $file.temp2 $file
	}
	file delete $file.temp
}
# cg download_biograph biograph_ad.tsv C0002395 
