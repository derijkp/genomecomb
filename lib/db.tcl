proc cg_tosqlite {args} {
	global scriptname action
	if {[llength $args] != 3} {
		puts stderr "format is: $scriptname $action dbfile table tsvfile"
		exit 1
	}
	foreach {dbfile table tsvfile} $args break
	set f [open $tsvfile]
	set header [split [gets $f] \t]
	set line [split [gets $f] \t]
	close $f
	set sql ""
	foreach field $header v $line {
		if {[isint $v]} {
			set type integer
		} elseif {[isdouble $v]} {
			set type real
		} else {
			set type text
		}
		lappend sql "\"$field\" $type"
	}
	set sql "create table $table ([join $sql ,])"
	catch {package require dbi}
	package require dbi_sqlite3
	dbi_sqlite3 db
	db create $dbfile
	db open $dbfile
	db exec $sql
	db import abort compar $tsvfile \t {}
	db close
}

