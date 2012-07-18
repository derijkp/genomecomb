Class subclass table_tsv

table_tsv method init args {
	private $object tdata cache
	if {[llength $args]} {
		$object open {*}$args
	}
	set tdata(numcache) 0
}

table_tsv method destroy {} {
	$object close
}

if 0 {
	set file /complgen/projects/kr1270/compar/annot_sv-kr1270.tsv.index/query_results.txt
	set field rowid
	set binfile /complgen/projects/kr1270/compar/annot_sv-kr1270.tsv.index/query_results.bincol
}

proc _table_tsv_calcline {object line cor} {
	if {[llength $cor] == 1} {
		return [list_sub $line [lindex $cor 0]]
	} elseif {[llength $cor]} {
		# $object.calcline is a function defined by method fields
		# when there are calculated fields
		return [$object.calcline $cor $line]
	} else {
		return $line
	}
}

proc _table_tsv_get {object row col args} {
	private $object tdata cache
	if {![info exists tdata(f)]} return
	if {$row == 0} {
		return [lindex [get tdata(qfields) {}] $col]
	}
	incr row -1
	if {![info exists cache($row)]} {
		set limit 50
		set lineindex $tdata(lineindex)
		if {$tdata(numcache) > 10000} {
			unset -nocomplain cache
			set tdata(numcache) 0
		}
		set cor $tdata(fieldscor)
		if {$tdata(query) eq ""} {
			set filepos [bcol_get $lineindex $row $row]
			if {$tdata(compressed)} {
				set f [gzopen $tdata(file) $filepos]
			} else {
				set f $tdata(f)
				seek $f $filepos
			}
			set r $row
			for {set i 0} {$i < $limit} {incr i} {
				if {[eof $f]} break
				set line [split [gets $f] \t]
				set cache($r) [_table_tsv_calcline $object $line $cor]
				incr r
			}
			if {$tdata(compressed)} {
				catch {close $f}
			}
		} else {
			set qrows [bcol_get $tdata(query_results) $row [expr {$row+$limit}]]
			set r $row
			foreach qrow $qrows {
				set filepos [bcol_get $lineindex $qrow $qrow]
				if {$tdata(compressed)} {
					set f [gzopen $tdata(file) $filepos]
				} else {
					set f $tdata(f)
					seek $f $filepos
				}
				set line [split [gets $f] \t]
				set cache($r) [_table_tsv_calcline $object $line $cor]
				if {$tdata(compressed)} {
					catch {close $f}
				}
				incr r
			}
		}
		incr tdata(numcache) $limit
	}
	return [lindex $cache($row) $col]
}

table_tsv method qlen {} {
	private $object tdata
	return $tdata(len)	
}

table_tsv method table {args} {
	private $object tdata
	return $tdata(table)	
}

table_tsv method query {query} {
	private $object tdata
	set tdata(query) $query
	if {[get tdata(query_results) ""] ne ""} {
		bcol_close $tdata(query_results)
	}
	if {$tdata(query) eq ""} {
		set tdata(query_results) {}
		set tdata(len) $tdata(tlen)
		$object reset
		Extral::event generate querychanged $object
		return
	}
	Classy::Progress start 2
	Classy::Progress message "Running query, please be patient (no progress shown)"
	putslog "Doing query $query"
	regsub -all \n $query { } query
	exec cg select -q $query -f {rowid=$ROW} $tdata(file) $tdata(indexdir)/query_results.tsv
	Classy::Progress next "Converting results"
	putslog "Converting results"
	set f [open $tdata(indexdir)/query_results.tsv]
	gets $f
	set o [open $tdata(indexdir)/query_results.bcol.bin.temp w]
	fconfigure $o -encoding binary -translation binary
	set len 0
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line] && [eof $f]} break
		incr len
		puts -nonewline $o [binary format i $line]
	}
	close $o
	set o [open $tdata(indexdir)/query_results.bcol.temp w]
	puts $o "# binary column"
	puts $o "# type iu"
	# puts $o "# len $len"
	puts $o "# [list query $query]"
#	puts $o [join {begin end type} \t]
	puts $o [join {begin type offset} \t]
	puts $o 0\tiu\t0
	puts $o [expr {$len-1}]\tend\t0
	close $o
	file rename -force $tdata(indexdir)/query_results.bcol.bin.temp $tdata(indexdir)/query_results.bcol.bin
	file rename -force $tdata(indexdir)/query_results.bcol.temp $tdata(indexdir)/query_results.bcol
	set tdata(len) $len
	$object reset
	set tdata(query_results) [bcol_open $tdata(indexdir)/query_results.bcol]
	putslog "query done"
	Classy::Progress stop
	Extral::event generate querychanged $object
}

#table_tsv method subsql {args} {
#	private $object tdata
#	set query [lindex $args 0]
#	set qfields $tdata(fields)
#	lappend qfields {*}[lrange $args 1 end]
#	set sql [monetdb_makesql $tdata(table) $tdata(tfields) $tdata(query) qfields {} 0 {} {}]
#	set fsql "with \"temp\" as ($sql) $query"
#	$object sql $fsql
#}
#
#table_tsv method sql {sql} {
#	private $object tdata
#	exec mclient -d $tdata(database) -f tab -s $sql
#}

table_tsv method info {key} {
	private $object tdata
	return $tdata($key)
}

table_tsv method tfields {} {
	private $object tdata
	return $tdata(tfields)
}

table_tsv method fields {args} {
	private $object tdata
	if {[llength $args]} {
		set tdata(fields) [lindex $args 0]
		if {[inlist {* {}} $tdata(fields)]}	{
			set tdata(qfields) $tdata(tfields)
			set tdata(fieldscor) {}
		} else {
			set header $tdata(tfields)
			set tdata(qfields) [tsv_select_expandfields $header $tdata(fields) tdata(qcode)]
			if {$tdata(qfields) eq $header} {
				set tdata(fieldscor) {}
			} else {
				set cor [list_cor $header $tdata(qfields)]
				set poss [list_find $cor -1]
				if {![llength $poss]} {
					set tdata(fieldscor) [list $cor]
				} else {
					set neededfields {}
					set num 0
					# neededfields will be the same for all procs, so first gather all we need using todo list
					# then make the procs from todo list later
					set todo {}
					foreach pos $poss {
						set el [lindex $tdata(qcode) $pos]
						set code [tsv_select_expandcode $header [lindex $el 1] neededfields]
						lappend todo $object.make_col$num $code
						incr num
					}
					# make procs for each calculated field
					set tclcode ""
					foreach {name code} $todo {
						append tclcode [subst -nocommands {
							proc $name {$neededfields} {
								if {[catch {expr {$code}} e]} {
									switch \$e {
										{domain error: argument not in valid range} {return NaN}
										{divide by zero} {return NaN}
									}
								}
								return \$e
							}
						}]
						lappend procs $name
					}
					# make main proc
					set neededcols [list_cor $header $neededfields]
					append tclcode "\nproc $object.calcline {cor line} \{\n"
					append tclcode "set result \[list_sub \$line [list $cor]\]\n"
					append tclcode "set neededfields \[list_sub \$line [list $neededcols]\]\n"
					foreach {name code} $todo pos $poss {
						append tclcode "lset result $pos \[$name {*}\$neededfields\]\n"
					}
					append tclcode "return \$result\n\}\n"
					uplevel #0 $tclcode
					set tdata(fieldscor) [list $cor $object.calcline]
				}
			}
		}
		$object reset
		Extral::event generate querychanged $object
	} else {
		return $tdata(fields)
	}
}

table_tsv method link {tktable} {
	private $object tdata
	set tdata(tktable)	$tktable
	$tktable configure -variabletype tktable -usecommand 1 -command "_table_tsv_get [list $object] %r %c"
	[winfo parent $tdata(tktable)].buttons.query configure -textvariable [privatevar $object tdata(query)]
}

table_tsv method index {file} {
	set time [file mtime $file]
	set indexdir [gzroot $file].index
	file mkdir $indexdir
	set indexfile $indexdir/lines.bcol
	set ext [file extension $file]
	if {[inlist {.rz .bgz .gz} $ext]} {set compressed 1} else {set compressed 0}
	cg_index $file
	set result [split [string trim [file_read $indexdir/info.tsv]] \n\t]
	dict set result lineindex [bcol_open $indexfile]
	dict set result file [file normalize $file]
	dict set result compressed $compressed
	dict set result lineindexfile [file normalize $indexfile]
	dict set result file $file
	lappend result indexdir $indexdir
	return $result
}

table_tsv method close {} {
	private $object tdata cache
	if {[get tdata(query_results) ""] ne ""} {
		bcol_close $tdata(query_results)
		set tdata(query_results) {}
		set tdata(query) {}
	}
	if {[get tdata(lineindex) ""] ne ""} {
		bcol_close $tdata(lineindex)
		set tdata(lineindex) {}
	}
	catch {close $tdata(f)}
	unset -nocomplain cache
	unset -nocomplain tdata
	set tdata(numcache) 0
}

table_tsv method open {file} {
	private $object tdata cache
	$object close
	set file [file normalize $file]
	set info [$object index $file]
	set file [dict get $info file]
	set database [file dir $file]
	set table $file
	setprivate $object database $database
	setprivate $object table $table
	set tdata(indexdir) [dict get $info indexdir]
	set tdata(lineindex) [dict get $info lineindex]
	set tdata(file) $file
	set tdata(database) $database
	set tdata(table) $table
	set tdata(tfields) [dict get $info header]
	set tdata(compressed) [dict get $info compressed]
	set tdata(tlen) [dict get $info size]
	set tdata(fields) $tdata(tfields)
	set tdata(qfields) $tdata(tfields)
	set tdata(fieldscor) {}
	set tdata(f) [open $file]
	if {![file exists $tdata(indexdir)/query_results.bcol]} {
		set tdata(query) {}
		set tdata(query_results) {}
		set tdata(len) $tdata(tlen)
	} else {
		set f [open $tdata(indexdir)/query_results.bcol]
		set header [tsv_open $f comment]
		close $f
		regsub -all {# } $comment { } comment
		set query [dict get $comment query]
		if {[file mtime $tdata(indexdir)/query_results.bcol] < [file mtime $file]} {
			puts "cache is older than file, redo query"
			$object query $query
		} else {
			set tdata(query) $query
			set tdata(query_results) [bcol_open $tdata(indexdir)/query_results.bcol]
			set tdata(len) [bcol_size $tdata(query_results)]
		}
	}
	$object reset
	Extral::event generate querychanged $object
	return $file
}

table_tsv method reset {args} {
	private $object tdata cache
	unset -nocomplain cache
	set tdata(numcache) 0
#	catch {
#		$tdata(tktable) configure -variabletype tktable -usecommand 1 -command "_table_tsv_get [list $object] %r %c"
#	}
}

table_tsv method save {file} {
	private $object tdata cache
	if {$tdata(query) eq ""} {
		file copy $tdata(file) $file
		return
	}
	set len $tdata(len)
	set row 0
	set lineindex $tdata(lineindex)
	set cor $tdata(fieldscor)
	set o [open $file w]
	puts $o "# [list orifile $tdata(file)]"
	puts $o "# [list query $tdata(query)]"
	puts $o [join $tdata(qfields) \t]
	while {$row < $len} {
		set qrows [bcol_get $tdata(query_results) $row [expr {$row+19}]]
		foreach qrow $qrows {
			set filepos [bcol_get $lineindex $qrow $qrow]
			if {$tdata(compressed)} {
				set f [gzopen $tdata(file) $filepos]
			} else {
				set f $tdata(f)
				seek $f $filepos
			}
			if {[llength $cor]} {
				puts $o [list_sub [split [gets $f] \t] $cor]
			} else {
				puts $o [gets $f]
			}
			if {$tdata(compressed)} {
				catch {close $f}
			}
			incr r
		}
		incr row 20
	}
	close $o
}
