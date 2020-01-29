proc cg_indexcol {args} {
	set samplingnum 0
	cg_options tsv2bed args {
		-sample {
			set samplingnum $value
		}
	} {file fieldname} 2
	set todofields [list $fieldname {*}$args]
	progress start [llength $todofields] "Indexing columns"
	foreach field $todofields {
		progress message "Indexing column $field"
		set histofile [indexdir_filewrite $file colinfo/$field.colinfo]
		if {$samplingnum == 0} {
			set fields [list "v=\[join \[split \$$field ,\\\\;\] \\\\n\]"]
			set filesize [file size $file]
			set step [expr {$filesize/10}]
			if {$step > 100000} {set step 100000} elseif {$step < 1} {set step 10000}
			progress start 2 "Indexing column" "Indexing column"
			progress start $filesize "Indexing column $field of file $file" "Indexing column $field"
			bgcg bgcg_progress bgexechandle \
				select -v -$step -f $fields -sh $histofile.temph $file $histofile.temp
			progress stop
			progress next "Sorting (no progress shown)"
			exec gnusort8 -N $histofile.temp | uniq > $histofile.temp2
			progress next "Sorting finished"
			progress stop
			file delete $histofile.temp $histofile.temph
			file rename -- $histofile.temp2 $histofile
		} else {
			progress start $samplingnum "Sampling data"  "Sampling $file to show example values for fields, please be patient"
			catch {update idletasks}
			if {![gziscompressed $file]} {
				set filesize [file size $file]
				set skip [expr {int($filesize/double($samplingnum))}]
			} else {
				set skip 0
			}
			catch {close $f} ; set f [gzopen $file]
			set header [tsv_open $f]
			set pos [lsearch $header $field]
			set num 0
			set row 0
			set break 0
			unset -nocomplain a
			set haslists 0
			while {![eof $f]} {
				set line [split [gets $f] \t]
				set valuelist [split [lindex $line $pos] {,;}]
				if {[llength $valuelist] > 1} {set haslists 1}
				foreach v $valuelist {
					set a($v) 1
				}
				incr num
				if {$num >= $samplingnum} {
					set break 1
					break
				}
				if {$skip} {
					seek $f $skip current
					gets $f
				}
				progress next
			}
			gzclose $f
			set result [array names a]
			set result [bsort $result]
			set result [tsv_defaultvalues $field $result]
			if {$haslists} {lappend result {present! lists}}
			lappend result {...}
			file_write $histofile.temp [join $result \n]
			file rename -force -- $histofile.temp $histofile
			progress stop
		}
		progress next
	}
	progress stop
}

proc cg_index {args} {
	set updated 0
	set cols 0
	set colinfo 0
	set dbstring {}
	set refdir {}
	cg_options index args {
		-cols {
			if {$value ni {0 1}} {set value 1; incr pos -1}
			set cols $value
		}
		-colinfo {
			if {$value ni {0 1}} {set value 1; incr pos -1}
			set colinfo $value
		}
		-db {
			set dbstring $value
		}
		-refdir {
			set refdir $value
		}
	} {file} 1 1
	set compressed [gziscompressed $file]
	set indexdir [indexdir $file]
	set infofile [indexdir_file $file info.tsv ok]
	if {$ok} {
		set info [infofile_read $infofile]
	} else {
		set info {}
	}
	if {$refdir ne ""} {
		set updated 1
		dict set info refdir $refdir
	}
	set indexfile [indexdir_file $file lines.bcol ok]
	set infofile [indexdir_file $file info.tsv infook]
	if {!$ok} {
		putslog "Creating lineindex"
		bcol_indexlines $file $indexfile $colinfo
		catch {file delete $indexdir/info.tsv}
		set f [gzopen $file]
		set header [tsv_open $f]
		catch {gzclose $f}
		set bcol [bcol_open $indexfile]
		set size [bcol_size $bcol]
		bcol_close $bcol
		dict set info file $file
		dict set info lineindexfile [file tail $indexfile]
		dict set info header $header
		dict set info size $size
		set updated 1
	} elseif {!$infook} {
		set f [gzopen $file]
		set header [tsv_open $f]
		catch {gzclose $f}
		set bcol [bcol_open $indexfile]
		set size [bcol_size $bcol]
		dict set info file $file
		dict set info lineindexfile [file tail $indexfile]
		dict set info header $header
		dict set info size $size
		set updated 1
	}
	if {$cols} {
		file mkdir $indexdir/cols
		set colsdir [indexdir_file $file cols ok]
		if {!$ok} {
			set f [gzopen $file]
			set header [tsv_open $f]
			set os {}
			foreach field $header {
				lappend os [open $colsdir/$field.col w]
			}
			while {![eof $f]} {
				set line [split [gets $f] \t]
				if {![llength $line] && [eof $f]} break
				foreach value $line o $os {
					puts $o $value
				}
			}
			gzclose $f
			foreach o $os {close $o}
			foreach field $header {
				exec gnusort8 -N $colsdir/$field.col | uniq -c | gnusort8 -n > $colsdir/$field.col.histo
			}
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
			putslog "Loading into database $db ($table)"
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
		infofile_write $infofile $info
	}
	return $indexfile
}
