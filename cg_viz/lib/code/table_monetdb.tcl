Class subclass table_monetdb

table_monetdb method init args {
	private $object tdata cache
	if {[llength $args]} {
		$object open {*}$args
	}
	set tdata(numcache) 0
}

proc _table_monetdb_get {object row col args} {
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
		set qfields rowid
		set limit 10
		set offset [expr {$row-1}]
		set sql [monetdb_makesql $tdata(table) $tdata(tfields) $tdata(query) qfields $sortfields $inverse $offset $limit]
		set ids [split [$object sql $sql] \n]
		set c [split [$object sql ""] \n]
		set qfields $tdata(fields)
		set limit 10
		set offset [expr {$row-1}]
		set qids "\"rowid\" in ([join $ids ,])"
		set sql [monetdb_makesql $tdata(table) $tdata(tfields) $qids qfields $sortfields 0 {} {}]
		set c [split [$object sql $sql] \n]
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

table_monetdb method query {args} {
	private $object tdata
	if {[llength $args]} {
		set tdata(query) [lindex $args 0]
	}
	unset -nocomplain cache
	set tdata(numcache) 0
	set qfields $tdata(fields)
	set sql [monetdb_makesql $tdata(table) $tdata(tfields) $tdata(query) qfields {} 0 {} {}]
	set sql "create view query_$tdata(table) as $sql"
	catch {$object sql [subst {drop view query_$tdata(table)}]}
	$object sql $sql
	set tdata(len) [$object sql [subst {select count(*) from "query_$tdata(table)"}]]
#	set qfields rowid
#	set sql [monetdb_makesql $tdata(table) $tdata(tfields) $tdata(query) qfields {} 0 {} {}]
#	regsub {\"rowid\"} $sql {count(*)} sql
#	set tdata(len) [$object sql $sql]
	$object reset
	Extral::event generate querychanged $object
}

table_monetdb method subsql {args} {
	private $object tdata
	set query [lindex $args 0]
	set qfields $tdata(fields)
	lappend qfields {*}[lrange $args 1 end]
	set sql [monetdb_makesql $tdata(table) $tdata(tfields) $tdata(query) qfields {} 0 {} {}]
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

table_monetdb method fields {args} {
	private $object tdata
	if {[llength $args]} {
		set tdata(fields) [lindex $args 0]
	}
	set tdata(qfields) [mselect_expandfields $tdata(tfields) $tdata(fields) tdata(qcode)]
	$object reset
	Extral::event generate querychanged $object
}

table_monetdb method link {tktable} {
	private $object tdata
	set tdata(tktable)	$tktable
	$tktable configure -variabletype tktable -usecommand 1 -command "_table_monetdb_get [list $object] %r %c"
	[winfo parent $tdata(tktable)].buttons.query configure -textvariable [privatevar $object tdata(query)]
}

table_monetdb method open {dbfarm database table} {
	private $object tdata cache
	unset -nocomplain cache
	set tdata(numcache) 0
	setprivate $object database $database
	setprivate $object table $table
	set tdata(database) $database
	set tdata(table) $table
	set tdata(tfields) [cg_monetdb fields $database $table]
	set tdata(tlen) [exec mclient -d $database -f tab -s [subst {select count(*) from \"$table\"}]]
	set tdata(fields) $tdata(tfields)
	set tdata(qfields) $tdata(tfields)
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
