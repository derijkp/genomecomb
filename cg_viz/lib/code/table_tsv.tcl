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
	if {![info exists tdata(sqlbackend_db)]} {
		exec cg select -q $query -f {rowid=$ROW} $tdata(file) $tdata(indexdir)/query_results.tsv
	} else {
		set qfields ROW
		set sql [monetdb_makesql $tdata(sqlbackend_table) $tdata(tfields) $query qfields rowid 0 $tdata(monetfieldtrans)]
		file_write $tdata(indexdir)/query_results.tsv ROW\n
		exec mclient -d $tdata(sqlbackend_db) -f tab -s $sql >> $tdata(indexdir)/query_results.tsv
		# exec cg mselect -q $query -f ROW $tdata(sqlbackend_db) $tdata(sqlbackend_table) $tdata(indexdir)/query_results.tsv
	}
	Classy::Progress next "Converting results"
	putslog "Converting results"
	set f [open $tdata(indexdir)/query_results.tsv]
	set header [tsv_open $f]
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

table_tsv method subsql {args} {
putsvars args
	private $object tdata
	if {![info exists tdata(sqlbackend_db)]} {
		error "no fast access without db backend"
	}
	set query [lindex $args 0]
	set qfields $tdata(fields)
#	set qfields $tdata(monetfields)
	lappend qfields {*}[lrange $args 1 end]
	set supersql [monetdb_makesql $tdata(sqlbackend_table) $tdata(tfields) $tdata(query) qfields {} 0 $tdata(monetfieldtrans)]
	set fsql "with \"temp\" as ($supersql) $query"
	$object sql $fsql
}

table_tsv method sql {sql {file {}}} {
	private $object tdata
	if {![info exists tdata(sqlbackend_db)]} {
		error "no fast access without db backend"
	}
	if {![info exists tdata(sqlbackend_db)]} {return {}}
	exec mclient -d $tdata(sqlbackend_db) -f tab -s $sql
}

table_tsv method info {key} {
	private $object tdata
	return $tdata($key)
}

table_tsv method tfields {} {
	private $object tdata
	return $tdata(tfields)
}

table_tsv method qfields {} {
	private $object tdata
	return $tdata(qfields)
}

table_tsv method values {field {max 5000}} {
	private $object tdata
	set histofile $tdata(file).index/colinfo/$field.colinfo
	if (![file exists $histofile]) {
		private $object values
		if {![info exists values($field)]} {
			set lineindex $tdata(lineindex)
			set step [expr {int([dict get $lineindex max]/double($max))}]
			if {$step < 1} {set step 1}
			catch {close $f} ; set f [gzopen $tdata(file)]
			set header [tsv_open $f]
			set pos [lsearch $header $field]
			set num 0
			set row 0
			set break 0
			unset -nocomplain a
			while {![eof $f]} {
				set line [split [gets $f] \t]
				set v [lindex $line $pos]
				if {![info exists a($v)]} {
					set a($v) 1
				} else {
					incr a($v)
				}
				incr num
				if {$num >= $max} {
					set break 1
					break
				}
				if {$step > 1} {
					incr row $step
					seek $f [bcol_get $lineindex $row $row]
				}
			}
			close $f
			set result {}
			foreach v [array names a] {
				lappend result [list $v $a($v)]
			}
			set result [lsort -dict -index 1 -decreasing $result]
			if {$break} {lappend result {sampled incomplete}}
			set values($field) $result
		}
		return $values($field)
	}
	set result {}
	set f [open $histofile]
	while {![eof $f]} {
		set line [gets $f]
		if {![llength $line]} continue
		if {[string index $line 0] eq "#"}		{
			if {[regexp {^# *([^: ]+): *(.*)$} $line temp key value]} {
				set a($key) $value
			}
		} else {
			lappend result [split $line \t]
		}
	}
	close $f
	if {[lindex $result end] eq "..."} {
		set a(incomplete) 1
		list_pop result
	} else {
		set a(incomplete) 0
	}
	set result [lsort -dict -index 1 -decreasing $result]
	if {[isdouble [get a(min) {}]] && [isdouble [get a(max) {}]]} {
		regsub {\.0+$} $a(min) {} min
		regsub {\.0+$} $a(max) {} max
		lappend result [list $max maxnum] [list $min minnum]
	}
	if {$a(incomplete)} {
		lappend result {... incomplete}
	}
	return $result
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

table_tsv method link {tktable button} {
	private $object tdata
	set tdata(tktable)	$tktable
	$tktable configure -variabletype tktable -usecommand 1 -command "_table_tsv_get [list $object] %r %c"
	$button configure -textvariable [privatevar $object tdata(query)]
}

table_tsv method index {file} {
	set time [file mtime $file]
	set indexdir [gzroot $file].index
	file mkdir $indexdir
	set indexfile $indexdir/lines.bcol
	set ext [file extension $file]
	if {[inlist {.rz .bgz .gz} $ext]} {set compressed 1} else {set compressed 0}
	cg_index $file
	set result [infofile_read $indexdir/info.tsv]
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

table_tsv method open {file parent} {
	private $object tdata cache
	$object close
	set file [file normalize $file]
	set info [$object index $file]
	unset -nocomplain tdata
	# tdata(query) is linked to entry, so clear it explicitely
	set tdata(query) ""
	array set tdata $info
	set tdata(file) $file
	set tdata(tfields) [dict get $info header]
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
	if {[info exists tdata(sqlbackend_db)] && [info exists tdata(refdir)]} {
		set tdata(monetfields) [split [cg_monetdb_fields $tdata(sqlbackend_db)  $tdata(sqlbackend_table)] \n]
		set tdata(monetfieldtrans) [cg_monetdb_sql $tdata(sqlbackend_db) [subst {select "value" from "genomecomb_info" where "table" = '$tdata(sqlbackend_table)' and "key" = 'fieldtrans'}]]
		catch {$object.disp1 destroy}
		display_chr new $object.disp1 $parent.canvas.data $object $tdata(refdir)
		Classy::todo $object.disp1 redraw
	}
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
	set cor $tdata(fieldscor)
	if {$tdata(query) eq ""} {
		if {![llength $cor]} {
			file copy $tdata(file) $file
		} else {
			cg select -f $tdata(fields) $tdata(file) $file
		}
		return
	}
	set len $tdata(len)
	set row 0
	set lineindex $tdata(lineindex)
	set o [open $file w]
	puts $o "# [list orifile $tdata(file)]"
	puts $o "# [list query $tdata(query)]"
	puts $o [join $tdata(qfields) \t]
	while {$row < $len} {
		set end [expr {$row+19}]
		if {$end >= $len} {set end [expr {$len-1}]}
		set qrows [bcol_get $tdata(query_results) $row $end]
		foreach qrow $qrows {
			set filepos [bcol_get $lineindex $qrow $qrow]
			if {$tdata(compressed)} {
				set f [gzopen $tdata(file) $filepos]
			} else {
				set f $tdata(f)
				seek $f $filepos
			}
			if {[llength $cor]} {
				set line [split [gets $f] \t]
				puts $o [join [_table_tsv_calcline $object $line $cor] \t]
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
