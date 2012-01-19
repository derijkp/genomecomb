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
				set cache($r) [split [gets $f] \t]
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
				set cache($r) [split [gets $f] \t]
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
	exec cg select -q $query -f {rowid=NR-1} $tdata(file) $tdata(indexdir)/query_results.tsv
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
	puts $o "# len $len"
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

table_tsv method fields {args} {
	private $object tdata
	if {[llength $args]} {
		set tdata(fields) [lindex $args 0]
		# set tdata(qfields) [mselect_expandfields $tdata(tfields) $tdata(fields) tdata(qcode)]
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
