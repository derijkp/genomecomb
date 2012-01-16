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

# lassign {/data/complgen/projects/kr1270/compar/annot_sv-kr1270.tsv.index/query_results.bcolinfo 1 4} bcolfile start end
proc bincol_get {bcolfile start end} {
	set f [open $bcolfile]
	set line [gets $f]
	if {$line ne "# binary column"} {error "$bcolfile is not a binary column file"}
	array set ainfo [string range [gets $f] 2 end]
	close $f
	if {$ainfo(type) ne "i"} {error "$bcolfile is not of the right type"}
	set max [expr {$ainfo(len)-1}]
	if {$start > $max} {
		return {}
	}
	if {$end > $max} {set end $max}
	set len [expr {$end-$start+1}]
	set pos [expr {4*$start}]
	set f [open [file root $bcolfile].bcolbin]
	fconfigure $f -encoding binary -translation binary
	seek $f $pos
	set b [read $f [expr {4*$len}]]
	binary scan $b i$len result
	close $f
	return $result
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
		set fi $tdata(fi)
		if {$tdata(numcache) > 10000} {
			unset -nocomplain cache
			set tdata(numcache) 0
		}
		if {$tdata(query) eq ""} {
			set ipos [expr {$row*8}]
			seek $fi $ipos
			unset -nocomplain filepos
			binary scan [read $fi 8] w filepos
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
			set indexfile $tdata(indexdir)/query_results.bcolinfo
			set qrows [bincol_get $indexfile $row [expr {$row+$limit}]]
			set r $row
			foreach qrow $qrows {
				set ipos [expr {$qrow*8}]
				seek $fi $ipos
				unset -nocomplain filepos
				binary scan [read $fi 8] w filepos
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

table_tsv method query {args} {
	private $object tdata
	if {[llength $args]} {
		set query [lindex $args 0]
		set tdata(query) $query
	} else {
		set query $tdata(query)
	}
	if {$tdata(query) eq ""} {
		Extral::event generate querychanged $object
		return
	}
	Classy::Progress start 2
	Classy::Progress message "Running query, please be patient (no progress shown)"
	putslog "Doing query $query"
	cg select -q $query -f {rowid=NR-1} $tdata(file) $tdata(indexdir)/query_results.tsv
	Classy::Progress next "Converting results"
	putslog "Converting results"
	set len [expr {[lindex [exec wc -l $tdata(indexdir)/query_results.tsv] 0] -1 }]
	set f [open $tdata(indexdir)/query_results.tsv]
	gets $f
	set o [open $tdata(indexdir)/query_results.bcolbin.temp w]
	fconfigure $o -encoding binary -translation binary
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line] && [eof $f]} break
		puts -nonewline $o [binary format i $line]
	}
	close $o
	set o [open $tdata(indexdir)/query_results.bcolinfo.temp w]
	puts $o "# binary column"
	puts $o "# [list query $query len $len type i]"
	puts $o [join {begin end type} \t]
	puts $o 0\t[expr {$len-1}]\ti
	close $o
	file rename -force $tdata(indexdir)/query_results.bcolbin.temp $tdata(indexdir)/query_results.bcolbin
	file rename -force $tdata(indexdir)/query_results.bcolinfo.temp $tdata(indexdir)/query_results.bcolinfo
	set tdata(len) $len
	$object reset
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
	}
	set tdata(qfields) [mselect_expandfields $tdata(tfields) $tdata(fields) tdata(qcode)]
	$object reset
	Extral::event generate querychanged $object
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
	set indexinfofile $indexdir/info.txt
	set indexlinesfile $indexdir/lines.binw
	set ext [file extension $file]
	if {[inlist {.rz .bgz .gz} $ext]} {set compressed 1} else {set compressed 0}
	if {![file exists $indexinfofile] || [file mtime $indexinfofile] < $time} {
		if {$compressed} {
			Classy::Progress start 1 "uncompressing $file for indexing, please be patient"
			Classy::Progress message "uncompressing $file for indexing, please be patient (no progress shown)"
			if {$ext eq ".rz"} {
				set tempfile $indexdir/[file root [file tail $file]]
				if {![file exists $tempfile]} {
					gunzip $file $tempfile
				}
			} else {
				gunzip $file
				set file [file root $file]
				set tempfile $file
			}
			Classy::Progress stop
		} else {
			set tempfile $file
		}
		set f [open $tempfile]
		set header [tsv_open $f comment]
		if {![file exists $indexlinesfile] || [file mtime $indexlinesfile] < $time} {
			Classy::Progress start [file size $tempfile] "Indexing $file, please be patient"
			Classy::Progress message "Indexing $file, please be patient"
			set o [open $indexlinesfile.temp w]
			fconfigure $o -encoding binary -translation binary
			while {![eof $f]} {
				set pos [tell $f]
				set line [gets $f]
				while {$line eq "" && ![eof $f]} {
					set pos [tell $f]
					set line [gets $f]
				}
				if {$line ne ""} {
					puts -nonewline $o [binary format w $pos]
				}
				Classy::Progress set $pos
			}
			close $o
			file rename -force $indexlinesfile.temp $indexlinesfile
			Classy::Progress stop
		}
		close $f
		if {$compressed} {
			if {$ext eq ".rz"} {
				file delete $tempfile
			} else {
				cg razip $file
				set file $file.rz
			}
		}
		set lines [expr {[file size $indexlinesfile]/8}]
		set o [open $indexinfofile.temp w]
		puts $o header\t$header
		puts $o comments\t$comment
		puts $o lines\t$lines
		close $o
		file rename -force $indexinfofile.temp $indexinfofile
	}
	set result [split [string trim [file_read $indexinfofile]] \t\n]
	lappend result file [file normalize $file]
	lappend result compressed $compressed
	lappend result indexlinesfile [file normalize $indexlinesfile]
	lappend result file $file
	lappend result indexdir $indexdir
	return $result
}

table_tsv method close {} {
	private $object tdata cache
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
	set tdata(file) $file
	set tdata(database) $database
	set tdata(table) $table
	set tdata(tfields) [dict get $info header]
	set tdata(compressed) [dict get $info compressed]
	set tdata(tlen) [dict get $info lines]
	set tdata(fields) $tdata(tfields)
	set tdata(qfields) $tdata(tfields)
	set tdata(len) $tdata(tlen)
	set tdata(query) {}
	set tdata(f) [open $file]
	set fi [open [dict get $info indexlinesfile]]
	fconfigure $fi -encoding binary -translation binary
	set tdata(fi) $fi
	$object reset
	Extral::event generate querychanged $object
	return $file
}

table_tsv method reset {args} {
	private $object tdata cache
	unset -nocomplain cache
	set tdata(numcache) 0
#	$tdata(tktable) configure -variabletype tktable -usecommand 1 -command "_table_tsv_get [list $object] %r %c"
}
