proc cg_index {args} {
	set updated 0
	set pos 0
	set len [llength $args]
	set cols 0
	set colinfo 0
	set dbstring {}
	set refdir {}
	set verbose 0
	while {$pos < $len} {
		set opt [lindex $args $pos]
		switch $opt {
			-cols {
				set cols 1
				incr pos
			}
			-colinfo {
				set colinfo 1
				incr pos
			}
			-db {
				incr pos
				set dbstring [lindex $args $pos]
				incr pos
			}
			-refdir {
				incr pos
				set refdir [lindex $args $pos]
				incr pos
			}
			-v {
				set verbose 1
				incr pos
			}
			-- break
			default {
				if {[string index $opt 0] ne "-"} break
				puts stderr "ERROR: Unkown option $opt: should be one of -cols, -db, -colinfo, -refdir, -v"
				exit 1
			}
		}
	}
	set args [lrange $args $pos end]
	if {[llength $args] != 1} {
		errorformat index
		exit 1
	}
	set file [lindex $args 0]
	set time [file mtime $file]
	set indexdir [gzroot $file].index
	set ext [file extension $file]
	if {[inlist {.rz .lz4 .bgz .gz} $ext]} {set compressed 1} else {set compressed 0}
	if {[file exists $indexdir/info.tsv]} {
		set info [infofile_read $indexdir/info.tsv]
	} else {
		set info {}
	}
	file mkdir $indexdir
	set indexfile $indexdir/lines.bcol
	if {$refdir ne ""} {
		set updated 1
		dict set info refdir $refdir
	}
	if {![file exists $indexdir/info.tsv] || [file mtime $indexdir/info.tsv] < $time || ![file exists $indexfile]} {
		if {$verbose} {
			putslog "Creating lineindex"
		}
		bcol_indexlines $file $indexfile $colinfo
		catch {file delete $indexdir/info.tsv}
		set f [gzopen $file]
		set header [tsv_open $f]
		catch {close $f}
		set bcol [bcol_open $indexfile]
		set size [bcol_size $bcol]
		bcol_close $bcol
		dict set info file $file
		dict set info lineindexfile [file tail $indexfile]
		dict set info header $header
		dict set info size $size
		set updated 1
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
	if {[dict exists $info sqlbackend_db]} {
		set db [dict get $info sqlbackend_db]
		set table [dict get $info sqlbackend_table]
		set user [dict get $info user]
		set pw [dict get $info pw]
	}
	if {$dbstring ne ""} {
		if {[llength $dbstring] < 2} {
			error "option -db must have format: database table ?user? ?pw?"
		}
		foreach {db table user pw} $dbstring break
	}
	if {[get db ""] ne "" || [dict exists $info sqlbackend_table]} {
		set updated 1
		dict set info sqlbackend_db $db 
		dict set info sqlbackend_table $table
		dict set info user $user
		dict set info pw $pw
		set error [catch {
			set dbtime [cg_monetdb_sql $db [subst {select "value" from "_genomecomb_info" where "table" = '$table' and "key" = 'time'}]]
		}]
		if {$error || ![isint $dbtime] || $dbtime < $time} {
			if {$verbose} {
				putslog "Loading into database $db ($table)"
			}
			cg_tomonetdb $db $table $file
			if {$user ne ""} {
				set o [open $indexdir/.monetdb w]
				puts $o "user=$user"
				puts $o "password=$pw"
				puts $o "language=sql"
				puts $o "save_history=true"
				close $o
			} else {
				file delete $indexdir/.monetdb
			}
		}
	}
	if {$updated} {
		infofile_write $indexdir/info.tsv $info
	}
	return $indexfile
}
