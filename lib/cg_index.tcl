proc cg_index {args} {
	if {[lindex $args 0] eq "-cols"} {
		set cols 1
		set args [lrange $args 1 end]
	} else {
		set cols 0
	}
	if {[llength $args] == 0} {
		errorformat index
		exit 1
	}
	set file [lindex $args 0]
	set time [file mtime $file]
	set indexdir [gzroot $file].index
	set ext [file extension $file]
	if {[inlist {.rz .bgz .gz} $ext]} {set compressed 1} else {set compressed 0}
	file mkdir $indexdir
	set indexfile $indexdir/lines.bcol
	bcol_indexlines $file $indexfile
	if {![file exists $indexdir/info.tsv] || [file mtime $indexdir/info.tsv] < $time} {
		catch {file delete $indexdir/info.tsv}
		set f [gzopen $file]
		set header [tsv_open $f]
		catch {close $f}
		set bcol [bcol_open $indexfile]
		set size [bcol_size $bcol]
		bcol_close $bcol
		set f [open $indexdir/info.tsv.temp w]
		puts $f key\tvalue
		puts $f file\t$file
		puts $f lineindexfile\t[file tail $indexfile]
		puts $f header\t$header
		puts $f size\t$size
		close $f
		file rename $indexdir/info.tsv.temp $indexdir/info.tsv
	}
	if {$cols} {
		file mkdir $indexdir/cols
		set f [gzopen $file]
		set header [tsv_open $f]
		set os {}
		foreach field $header {
			lappend os [open $indexdir/cols/$field.col w]
		}
		while {![eof $f]} {
			set line [split [gets $f] \t]
			if {![llength $line] && [eof $f]} break
			foreach value $line o $os {
				puts $o $value
			}
		}
		close $f
		foreach o $os {close $o}
		foreach field $header {
			exec gnusort8 -N $indexdir/cols/$field.col | uniq -c | gnusort8 -n > $indexdir/cols/$field.col.histo
		}
	}
	return $indexfile
}
