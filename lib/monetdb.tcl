proc monetdbinit {} {
	global monetdbdir
	if {![catch {set monetdbdir [file dir [file dir [exec which mclient]]]} e]} {
		return
	}
	foreach path [list \
		$::appdir/../extern/monetdb \
		$::appdir/../../extern/monetdb \
		/home/peter/apps/monetdb \
	] {
		if {[file exists $path/bin/mclient]} {
			set monetdbdir $path
			set ::env(PATH) ${::env(PATH)}:$path/bin
			break
		}
	}
	if {![info exists monetdbdir]} {
		error "could not find monetdb installation"
	}
}
monetdbinit

proc monetdb_getinfo {dbfarm} {
	set dbfarm [file normalize $dbfarm]
	set data [exec monetdbd get all $dbfarm]
	set result {}
	foreach line [split $data \n] {
		if {[regexp {^ *([^ ]+) +(.*) *$} $line temp key value]} {
			dict set result $key $value
		}
	}
	return $result
}

proc cg_monetdb {args} {
	# work in progress
	set cmd [lindex $args 0]
	set args [lrange $args 1 end]
	switch $cmd {
		sql {
			foreach {database sql} $args break
			set result [exec mclient -d $database -f tab -s $sql >@ stdout]
		}
		createdbfarm {
			foreach {dbfarm port} $args break
			# create
			exec monetdbd create $dbfarm
			exec monetdbd set port=$port $dbfarm
		}
		startdbfarm {
			foreach {dbfarm} $args break
			exec monetdbd start $dbfarm
		}
		stopdbfarm {
			foreach {dbfarm} $args break
			exec monetdbd start $dbfarm
		}
		createdb {
			foreach {dbfarm database} $args break
			set controlport [dict get [monetdb_getinfo $dbfarm] controlport]
			exec monetdb -p$controlport create $database
			exec monetdb -p$controlport start $database
			exec monetdb -p$controlport release $database
		}
		setuser {
			foreach {user password} $args break
			file_write ~/.monetdb "user=$user\npassword=$password\nlanguage=sql\nsave_history=true\n"
		}
		schema {
			foreach {database} $args break
			set sql {SELECT "current_schema"}
			set result [exec mclient -d $database -f tab -s $sql]
		}
		tables {
			foreach {database schema} $args break
			if {$schema eq ""} {
				set schema sys
			}
			set sql [format {
				SELECT t.name
				FROM "sys"."_tables" "t","sys"."schemas" "s"
			      	WHERE "t"."system" <> true
				AND "t"."schema_id" = "s"."id"
				AND "s"."name" = '%s'
			} $schema]
			set result [exec mclient -d $database -f tab -s $sql]
		}
		fields {
			foreach {database table schema} $args break
			if {$schema eq ""} {
				set schema sys
			}
			set sql [format {
				SELECT "c"."name"
				FROM "sys"."_columns" "c",
					"sys"."_tables" "t",
				 	"sys"."schemas" "s"
			 	WHERE "c"."table_id" = "t"."id"
				 AND '%s' = "t"."name"
				 AND "t"."schema_id" = "s"."id"
				 AND "s"."name" = '%s'
				 ORDER BY "number"
			} $table $schema]
			set result [exec mclient -d $database -f tab -s $sql]
		}
		tableinfo {
			foreach {database table schema} $args break
			if {$schema eq ""} {
				set schema sys
			}
			set sql [format {
				SELECT "c"."name",
				"c"."type",
				"c"."type_digits",
				"c"."type_scale",
				"c"."null",
				"c"."default",
				"c"."number"
				FROM "sys"."_columns" "c",
					"sys"."_tables" "t",
				 	"sys"."schemas" "s"
			 	WHERE "c"."table_id" = "t"."id"
				 AND '%s' = "t"."name"
				 AND "t"."schema_id" = "s"."id"
				 AND "s"."name" = '%s'
				 ORDER BY "number"
			} $table $schema]
			set result [exec mclient -d $database -f tab -s $sql]
		}
		default {
			error "unkown command \"$cmd\""
		}
	}
	return $result
}

proc cg_tomonetdb {args} {
	global scriptname action
	if {[llength $args] != 3} {
		puts stderr "format is: $scriptname $action db table tsvfile"
		exit 1
	}
	foreach {db table tsvfile} $args break
	set tsvfile [file normalize $tsvfile]
	set f [gzopen $tsvfile]
	set header [tsv_open $f]
	set o [open $tsvfile.temp w]
	fcopy $f $o
	close $o
	close $f
	array set ftype {begin integer end integer start integer type varchar(10)}
	set sql ""
	foreach field $header {
		set type [get ftype($field) clob]
		regsub -- {\.} $field {x} field
		regsub -- {_} $field {x} field
		regsub -- {-} $field {_} field
		regsub -- {\+} $field {x} field
		lappend sql "\"$field\" $type"
	}
	set sql "create table \"$table\" ([join $sql ,\n]);\n"
	set num [lindex [exec wc -l $tsvfile.temp] 0]
	exec mclient -d$db -s $sql
	exec mclient -d$db -s "copy $num records into \"$table\" from '$tsvfile.temp' delimiters '\t', '\n' null as '';"
	exec mclient -d$db -s "alter table \"$table\" add column \"rowid\" serial"
	file delete $tsvfile.temp
}

if 0 {
	set dbfarm ~/my-dbfarm
	set database test
	set port 50000

set db test
set table compar
set tsvfile /media/wd2t/complgen/projects/dlb1/dlb_compar.tsv
set tsvfile /complgen/projects/amg_md200/compar/annotamg_md200_compar.tsv

}
