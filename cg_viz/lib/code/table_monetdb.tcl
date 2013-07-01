Class subclass table_monetdb

table_monetdb method init args {
	private $object tdata cache
	if {[llength $args]} {
		$object open {*}$args
	}
	set tdata(numcache) 0
}

proc _table_monetdb_get {object row col args} {
# putsvars object row col args
	private $object tdata cache
	if {$row == 0} {
		return [lindex [get tdata(qfields) {}] $col]
	}
	if {![info exists cache($row)]} {
		if {$tdata(numcache) > 10000} {
			unset -nocomplain cache
			set tdata(numcache) 0
		}
		set sortfields {}
		set inverse 0
		set qfields $tdata(fields)
		set limit 50
		set offset [expr {$row-1}]
		set end [expr {$row+$limit}]
		if {$tdata(query) eq ""} {
			set query [subst {"rowid" between $offset and $end}]
			set sql [monetdb_makesql $tdata(table) $tdata(tfields) $query qfields [list_union $sortfields rowid]]
			set c [split [$object sql $sql] \n]
		} else {
			set sql [subst {select * from "query" where "query_rowid" between $offset and $end}]
			set c [split [$object sql $sql] \n]
		}
		set r $row
		foreach line $c {
			set cache($r) [split $line \t]
			incr r
		}
		incr tdata(numcache) $limit
	}
	return [lindex $cache($row) $col]
}

table_monetdb method qlen {} {
	private $object tdata
	return $tdata(len)	
}

table_monetdb method table {args} {
	private $object tdata
	return $tdata(table)	
}

#table_monetdb method query {args} {
#	private $object tdata
#	if {[llength $args]} {
#		set tdata(query) [lindex $args 0]
#	}
#	unset -nocomplain cache
#	set tdata(numcache) 0
#	catch {$object sql [subst {drop view query}]}
#	set qfields $tdata(fields)
#	lappend qfields {query_rowid=row_number() over (order by rowid)}
#	set sql [monetdb_makesql $tdata(table) $tdata(tfields) $tdata(query) qfields {} 0]
#	set sql "create view query as $sql"
#	$object sql $sql
#	set tdata(len) [$object sql [subst {select count(*) from "query"}]]
#	$object reset
#	Extral::event generate querychanged $object
#}

table_monetdb method query {args} {
	private $object tdata
	if {[llength $args]} {
		set tdata(query) [lindex $args 0]
	}
	unset -nocomplain cache
	set tdata(numcache) 0
	catch {$object sql [subst {drop view query}]}
	set qfields $tdata(fields)
	lappend qfields {query_rowid=row_number() over (order by rowid)}
	set sql [monetdb_makesql $tdata(table) $tdata(tfields) $tdata(query) qfields rowid 0]
	catch {$object sql {drop view "query"}}
	set sql "create view \"query\" as $sql"
	$object sql $sql
	set tdata(len) [$object sql [subst {select count(*) from "query"}]]
	$object reset
	Extral::event generate querychanged $object
}

table_monetdb method subsql {args} {
	private $object tdata
	set query [lindex $args 0]
	set qfields $tdata(fields)
	lappend qfields {*}[lrange $args 1 end]
	set sql [monetdb_makesql $tdata(table) $tdata(tfields) $tdata(query) qfields {} 0]
	set fsql "with \"temp\" as ($sql) $query"
	$object sql $fsql
}

table_monetdb method sql {sql} {
	private $object tdata
	exec mclient -d $tdata(database) -f tab -s $sql
}

table_monetdb method info {key} {
	private $object tdata
	return $tdata($key)
}

table_monetdb method tfields {} {
	private $object tdata
	return $tdata(tfields)
}

table_monetdb method qfields {} {
	private $object tdata
	return [get tdata(qfields) ""]
}

table_monetdb method fields {args} {
	private $object tdata
	if {[llength $args]} {
		set tdata(fields) [lindex $args 0]
		if {[inlist {* {}} $tdata(fields)]}	{
			set tdata(qfields) $tdata(tfields)
			set tdata(fieldscor) {}
		} else {
			set tdata(qfields) [tsv_select_expandfields $tdata(tfields) $tdata(fields) tdata(qcode)]
			if {$tdata(qfields) eq $tdata(tfields)} {
				set tdata(fieldscor) {}
			} else {
				set tdata(fieldscor) [list_cor $tdata(tfields) $tdata(qfields)]
			}
		}
		$object reset
		Extral::event generate querychanged $object
	} else {
		return $tdata(fields)
	}
}

table_monetdb method link {tktable button} {
	private $object tdata
	set tdata(tktable)	$tktable
	$tktable configure -variabletype tktable -usecommand 1 -command "_table_monetdb_get [list $object] %r %c"
	$button configure -textvariable [privatevar $object tdata(query)]
}

table_monetdb method open {database table dbfarm} {
	private $object tdata cache
	unset -nocomplain cache
	set tdata(numcache) 0
	setprivate $object database $database
	setprivate $object table $table
	set tdata(database) $database
	set tdata(table) $table
	# table fields (fields in original table)
	set tdata(tfields) [cg_monetdb_fields $database $table]
	# fields (requested fields for query, may include calculated fields in the form field=...)
	set tdata(fields) $tdata(tfields)
	# query fields (requested fields for query)
	set tdata(qfields) $tdata(tfields)
	set tdata(tlen) [exec mclient -d $database -f tab -s [subst {select count(*) from \"$table\"}]]
	set tdata(len) $tdata(tlen)
	set tdata(query) {}
	$object reset
	Extral::event generate querychanged $object
}

table_monetdb method reset {args} {
	private $object tdata cache
	unset -nocomplain cache
	set tdata(numcache) 0
#	$tdata(tktable) configure -variabletype tktable -usecommand 1 -command "_table_monetdb_get [list $object] %r %c"
}
