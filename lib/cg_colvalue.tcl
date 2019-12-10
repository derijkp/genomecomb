proc cg_colvalue {args} {
	set file ""
	set outfile ""
	set idfields {}
	set keyfields {}
	set valuefield {}
	cg_options colvalue args {
		-idfields {
			set idfields $value
		}
		-keyfields {
			set keyfields $value
		}
		-valuefield {
			set valuefield $value
		}
	} {file outfile} 0 2
	if {$file eq ""} {
		set f stdin
		set header [tsv_open $f comment]
		set tempfile [tempfile]
		set o [open $tempfile w]
		puts $o [join $header \t]
		fcopy $f $o
		close $o
	} else {
		set tempfile $file
		set f [gzopen $file]
		set header [tsv_open $f comment]
		gzclose $f
	}
	if {$outfile eq ""} {
		set o stdout
	} else {
		set tempoutfile [filetemp $outfile]
		set o [open $tempoutfile w]
	}
	if {![llength $keyfields]} {
		foreach keyfields {key parameter} {
			set keyposs [lsearch $header $keyfields]
			if {$keyposs != -1} break
		}
		if {$keyposs == -1} {error "keyfield \"key\" or \"parameter\" not found"}
	} else {
		set keyposs [list_cor $header $keyfields]
		if {[inlist $keyposs -1]} {
			error "Some keyfield(s) not present: [list_sub $keyfields [list_find $keyposs -1]]"
		}
	}
	if {![llength $valuefield]} {
		foreach valuefield {value} {
			set valuepos [lsearch $header $valuefield]
			if {$valuepos != -1} break
		}
		if {$valuepos == -1} {error "valuefield \"value\" not found"}
	} else {
		set valuepos [lsearch $header $valuefield]
		if {$valuepos == -1} {
			error "valuefield not present: $valuefield"
		}
	}
	if {![llength $idfields]} {
		set idfields [list_lremove [list_remove $header $valuefield] $keyfields]
	}
	set idposs [list_cor $header $idfields]
	set fields {}
	foreach line [lrange [split [cg select -g [list_merge $keyfields {}] $tempfile] \n] 1 end] {
		lappend fields [join [lrange $line 0 end-1] -]
	}
	#
	# make output
	set f [gzopen $tempfile]
	set header [tsv_open $f]
	puts $o [join [list {*}$idfields {*}$fields] \t]
	set previd {}
	unset -nocomplain a
	while 1 {
		set line [split [gets $f] \t]
		set id [list_sub $line $idposs]
		if {$id ne $previd && $previd ne "" || [eof $f]} {
			set result $previd
			foreach field $fields {
				lappend result [get a($field) ""]
			}
			puts $o [join $result \t]
			unset -nocomplain a
			if {[eof $f]} break
		}
		if {![llength $line]} continue
		set key [join [list_sub $line $keyposs] -]
		set value [lindex $line $valuepos]
		set a($key) $value
		set previd $id
	}
	#
	# finish
	if {$o ne "stdout"} {
		close $o
		file rename -force -- $tempoutfile $outfile
	}
	gzclose $f
}

if 0 {
	set tsvfile tests/data/reg4.tsv
	set tsvfile tests/data/testvars.tsv
	set o stdout
}
