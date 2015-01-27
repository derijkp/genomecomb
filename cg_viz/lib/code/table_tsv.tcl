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

proc _table_tsv_calcline {object line cor row} {
	if {[llength $cor] == 1} {
		return [list_sub $line [lindex $cor 0]]
	} elseif {[llength $cor]} {
		# $object.calcline is a function defined by method fields
		# when there are calculated fields
		return [$object.calcline $cor $line $row]
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
				set cache($r) [_table_tsv_calcline $object $line $cor $r]
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
				set cache($r) [_table_tsv_calcline $object $line $cor $qrow]
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

table_tsv method query {args} {
	private $object tdata
	if {![llength $args]} {
		return $tdata(query)
	}
	set query [lindex $args 0]
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
	progress start {70 30}
	putslog "Doing query $query"
	regsub -all \n $query { } query
	if {![info exists tdata(sqlbackend_db)]} {
		if {[file exists $tdata(indexdir)/info.tsv]} {
			set info [infofile_read $tdata(indexdir)/info.tsv]
			catch {set numlines [dict get $info size]}
		}
		set fieldopt [list {rowid=$ROW}]
		if {$query ne ""} {
			set header $tdata(header)
			set fields $tdata(fields)
			tsv_select_expandcode $header $query neededfields
			set neededfields [list_lremove [list_remdup $neededfields] $header]
			foreach field $neededfields {
				set pos [lsearch -regexp $fields " *-?$field *="]
				if {$pos == -1} continue
				set temp [lindex $fields $pos]
				if {[string index $temp 0] != "-"} {set temp -$temp}
				lappend fieldopt $temp
			}
		}
		if {![info exists numlines]} {
			progress message "Running query, please be patient (no progress shown)"
			exec cg select -q $query -f $fieldopt $tdata(file) $tdata(indexdir)/query_results.tsv
		} else {
			set step [expr {$numlines/10}]
			if {$step > 50000} {set step 50000} elseif {$step < 1} {set step 1}
			progress start $numlines "Running query" "Running query"
			set ::bgerror {}
			Extral::bgexec -progresscommand [list $object queryprogress] -no_error_redir -channelvar [privatevar $object bgexechandle] \
				cg select -v $step -q $query -f $fieldopt $tdata(file) $tdata(indexdir)/query_results.tsv 2>@1
			if {$::bgerror ne ""} {error $::bgerror}
			progress stop
		}
	} else {
		progress message "Running query, please be patient (no progress shown)"
		set qfields ROW
		set sql [monetdb_makesql $tdata(sqlbackend_table) $tdata(tfields) $query qfields rowid 0 $tdata(monetfieldtrans)]
		file_write $tdata(indexdir)/query_results.tsv ROW\n
		exec mclient -d $tdata(sqlbackend_db) -f tab -s $sql >> $tdata(indexdir)/query_results.tsv
		# exec cg mselect -q $query -f ROW $tdata(sqlbackend_db) $tdata(sqlbackend_table) $tdata(indexdir)/query_results.tsv
	}
	progress next "Converting results"
	putslog "Converting results"
	cg bcol make -t iu -co 0 $tdata(indexdir)/tempquery_results rowid < $tdata(indexdir)/query_results.tsv
	set c [file_read $tdata(indexdir)/tempquery_results.bcol]
	file_write $tdata(indexdir)/tempquery_results.bcol [string_change $c [list "\# default 0" "\# default 0\n\# [list query $query]"]]
	set len [expr {[lindex [split [string trim $c] \n] end 0]+1}]
	file rename -force $tdata(indexdir)/tempquery_results.bcol.bin $tdata(indexdir)/query_results.bcol.bin
	file rename -force $tdata(indexdir)/tempquery_results.bcol $tdata(indexdir)/query_results.bcol
	set tdata(len) $len
	$object reset
	set tdata(query_results) [bcol_open $tdata(indexdir)/query_results.bcol]
	putslog "query done"
	progress stop
	Extral::event generate querychanged $object
}

table_tsv method subsql {args} {
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

proc tsv_defaultvalues {field {result {}}} {
	if {[regexp ^sequenced- $field]} {
		set result [list_union {r v u} $result]
	} elseif {[regexp ^zyg- $field]} {
		set result [list_union {m t c o r u} $result]
	} elseif {[regexp _impact $field]} {
		annot_init
		set pre [list_subindex $::snp_annot_list 0]
		set result [list_union $pre $result]
	} elseif {$field eq "type"} {
		annot_init
		set result [list_union {snp del ins sub} $result]
	}
	return $result
}

table_tsv method values {field {samplevalues 0} {max 1000} {valuesonly 0}} {
	private $object tdata values
	switch $samplevalues {
		all {cg_indexcol $tdata(file) $field}
		allif0 {
			set histofile [indexdir_file $tdata(file) colinfo/$field.colinfo ok]
			if (!$ok) {cg_indexcol $tdata(file) $field}
		}
		sample - 1 {
			if {$tdata(size) < $max || $tdata(size) < 10000} {
				cg_indexcol $tdata(file) $field
			} else {
				cg_indexcol -sample $max $tdata(file) $field
			}
		}
	}
	set histofile [indexdir_file $tdata(file) colinfo/$field.colinfo ok]
	if (!$ok) {return {}}
	if {[file size $histofile] > 1000000} {
		set skip [expr {[file size $histofile]/$max}]
		set sample 1
	} else {
		set sample 0
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
			if {$sample} {
				seek $f $skip current
				gets $f
			}
		}
	}
	close $f
	if {[lindex $result end] eq "..."} {
		set a(incomplete) 1
		list_pop result
	} else {
		set a(incomplete) 0
	}
	if {$sample} {
		set a(incomplete) 1
	}
	set result [ssort -natural [list_subindex $result 0]]
	if {$valuesonly} {return $result}
	if {[isdouble [get a(min) {}]] && [isdouble [get a(max) {}]]} {
		regsub {\.0+$} $a(min) {} min
		regsub {\.0+$} $a(max) {} max
		lappend result [list $max maxnum] [list $min minnum]
	}
	if {$sample} {
		lappend result {{to large} incomplete}
	} elseif {$a(incomplete)} {
		lappend result {sampled incomplete}
	}
	return $result
}

table_tsv method fields {args} {
	private $object tdata
	if {[llength $args]} {
		set tdata(fields) [lindex $args 0]
		file_write $tdata(indexdir)/query_fields.txt $tdata(fields)
		if {[inlist {* {}} $tdata(fields)]}	{
			set tdata(qfields) $tdata(tfields)
			set tdata(fieldscor) {}
		} else {
			set header $tdata(tfields)
			set tdata(qfields) [tsv_select_expandfields $header $tdata(fields) tdata(qcode)]
			set pos 0
			set poss {}
			foreach el $tdata(qcode) {
				if {![isint $el]} {lappend poss $pos}
				incr pos
			}
			if {$tdata(qfields) eq $header && ![llength $poss]} {
				set tdata(fieldscor) {}
			} else {
				set cor [list_cor $header $tdata(qfields)]
				if {![llength $poss]} {
					set tdata(fieldscor) [list $cor]
				} else {
					set neededfields {}
					set num 0
					# neededfields will be the same for all procs, so first gather all we need using todo list
					# then make the procs from todo list later
					set todo {}
					set prequery {}
					foreach pos $poss {
						set el [lindex $tdata(qcode) $pos]
						set code [tsv_select_expandcode $header [lindex $el 1] neededfields prequery]
						lappend todo $object.make_col$num $code
						incr num
					}
					# make procs for each calculated field
					set tclcode ""
					foreach {name code} $todo {
						append tclcode [subst -nocommands {
							proc $name {$neededfields} {
								$prequery
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
					append tclcode "\nproc $object.calcline {cor line ROW} \{\n"
					append tclcode "set result \[list_sub \$line [list $cor]\]\n"
					append tclcode "set neededfields \[list_sub \$line [list $neededcols]\]\n"
					foreach {name code} $todo pos $poss {
						append tclcode "lset result $pos \[$name {*}\$neededfields\]\n"
					}
					set rowpos [lsearch $neededfields ROW]
					if {$rowpos != -1} {
						append tclcode "lset result $rowpos \$ROW\n"
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

table_tsv method index {file {colinfo 0}} {
	set time [file mtime $file]
	set indexfile [indexdir_file $file lines.bcol ok]
	set ext [file extension $file]
	if {[inlist {.rz .lz4 .bgz .gz} $ext]} {set compressed 1} else {set compressed 0}
	if {$colinfo} {
		cg_index -colinfo $file
	} else {
		cg_index $file
	}
	set infofile [indexdir_file $file info.tsv ok]
	if {$ok} {
		set result [infofile_read $infofile]
	} else {
		cg_index $file
		set infofile [indexdir_file $file info.tsv ok]
		set result [infofile_read $infofile]
	}
	set indexdir [indexdir $file]
	puts "Using indexdir $indexdir"
	dict set result indexdir $indexdir
	dict set result lineindex [bcol_open $indexfile]
	dict set result file [file_absolute $file]
	dict set result compressed $compressed
	dict set result lineindexfile [file_absolute $indexfile]
	dict set result file $file
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
	set file [file_absolute $file]
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
	tsv_select_sampleinfo_setfile $file
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
			if {[catch {
				$object query $query
			} errormsg]} {
				Classy::todo bgerror $errormsg
			}
		} else {
			set tdata(query) $query
			set tdata(query_results) [bcol_open $tdata(indexdir)/query_results.bcol]
			set tdata(len) [bcol_size $tdata(query_results)]
		}
	}
	$object reset
	if {[info exists tdata(sqlbackend_db)] && [info exists tdata(refdir)]} {
		set tdata(monetfields) [split [cg_monetdb_fields $tdata(sqlbackend_db)  $tdata(sqlbackend_table)] \n]
		set tdata(monetfieldtrans) [cg_monetdb_sql $tdata(sqlbackend_db) [subst {select "value" from "_genomecomb_info" where "table" = '$tdata(sqlbackend_table)' and "key" = 'fieldtrans'}]]
		catch {$object.disp1 destroy}
		display_chr new $object.disp1 $parent.canvas.data $object $tdata(refdir)
		Classy::todo $object.disp1 redraw
	}
	if {[file exists $tdata(indexdir)/query_fields.txt]} {
		if {[catch {
			$object fields [file_read $tdata(indexdir)/query_fields.txt]
		} errormsg]} {
			Classy::todo bgerror $errormsg
		}
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
				puts $o [join [_table_tsv_calcline $object $line $cor $row] \t]
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
