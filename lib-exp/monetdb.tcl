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

proc cg_monetdb_tables {database {schema sys}} {
	set sql [format {
		SELECT t.name
		FROM "sys"."_tables" "t","sys"."schemas" "s"
	      	WHERE "t"."system" <> true
		AND "t"."schema_id" = "s"."id"
		AND "s"."name" = '%s'
	} $schema]
	set result [exec mclient -d $database -f tab -s $sql]
}

proc cg_monetdb_fields {database table {schema sys}} {
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

proc cg_monetdb_schema {database} {
	set sql {SELECT "current_schema"}
	set result [exec mclient -d $database -f tab -s $sql]
}

proc cg_monetdb_tableinfo {database table {schema sys}} {
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

proc cg_monetdb_initdb database {
	catch {exec mclient -d $database -s {
		-- case sensitive
		create function pcre_match(s string, pattern string)
		  returns BOOLEAN
		  external name pcre.match;
	}} e
	puts "adding pcre_match: $e"
	catch {exec mclient -d $database -s {
		-- case insensitive
		create function pcre_imatch(s string, pattern string)
		  returns BOOLEAN
		  external name pcre.imatch;
	}} e
	puts "adding pcre_imatch: $e"
	catch {exec mclient -d $database -s {
		create table genomecomb_info (
			"table" clob,
			"key" clob,
			"value" clob,
			primary key("table","value")
		)
	}} e
	puts "Creating table genomecomb_info: $e"
}

proc cg_monetdb_sql {database sql} {
	set result [exec mclient -d $database -f tab -s $sql]
}

proc cg_monetdb {args} {
	# work in progress
	set cmd [lindex $args 0]
	set args [lrange $args 1 end]
	set result {}
	switch $cmd {
		sql {
			foreach {database sql} $args break
			exec mclient -d $database -f tab -s $sql >@ stdout
			return
		}
		createdbfarm {
			foreach {dbfarm port} $args break
			# create
			exec monetdbd create $dbfarm
			if {[llength $args] > 1} {
				exec monetdbd set port=$port $dbfarm
			}
		}
		startdbfarm {
			foreach {dbfarm port} $args break
			set dbfarm [file normalize $dbfarm]
			if {[llength $args] > 1} {
				exec monetdbd set port=$port $dbfarm
			}
			set info [monetdb_getinfo $dbfarm]
			if {[dict get $info status] eq "no monetdbd is serving this dbfarm"} {
				puts "Starting $dbfarm"
				exec monetdbd start $dbfarm
			} else {
				puts "$dbfarm already running on port: [dict get $info port]"
			}
			set info [monetdb_getinfo $dbfarm]
			return $info
		}
		stopdbfarm {
			foreach {dbfarm} $args break
			exec monetdbd stop $dbfarm
		}
		port {
			foreach {dbfarm} $args break
			set result [dict get [monetdb_getinfo $dbfarm] port]
		}
		stopdb {
			set database [list_pop args]
			exec monetdb {*}$args stop $database
		}
		startdb {
			set database [list_pop args]
			exec monetdb {*}$args start $database
		}
		createdb {
			set database [list_pop args]
			# set controlport [dict get [monetdb_getinfo $dbfarm] port]
			# exec monetdb -p$controlport create $database
			## catch {exec monetdb -p$controlport release $database}
			#exec monetdb -p$controlport start $database
			exec monetdb {*}$args create $database
			exec monetdb {*}$args start $database
			cg_monetdb_initdb $database
			exec monetdb {*}$args release $database
		}
		initdb {
			set result [cg_monetdb_initdb {*}$args]
		}
		setuser {
			# set DOTMONETDBFILE to use other file than ~/.monetdb
			unset ::env(DOTMONETDBFILE)
			set tempfile [tempfile]
			foreach {user password} $args break
			file_write ~/.monetdb [subst {
				user=$user
				password=$password
				language=sql
				save_history=true
			}]
		}
		settempuser {
			# set DOTMONETDBFILE to use other file than ~/.monetdb
			set tempfile [tempfile]
			foreach {user password} $args break
			file_write $tempfile [subst {
				user=$user
				password=$password
				language=sql
				save_history=true
			}]
			file attributes $tempfile -permissions go-rwx
			set ::env(DOTMONETDBFILE) $tempfile
		}
		schema {
			set result [cg_monetdb_schema {*}$args]
		}
		tables {
			set result [cg_monetdb_tables {*}$args]
		}
		fields {
			set result [cg_monetdb_fields {*}$args]
		}
		tableinfo {
			set result [cg_monetdb_tableinfo {*}$args]
		}
		default {
			error "cg monetdb unkown sub command \"$cmd\""
		}
	}
	puts $result
}

proc monetdb_convfield field {
	set poss [regexp -all -indices -inline {[^a-zA-Z0-9-]+} $field]
	set newfield $field
	list_foreach {s e} [lreverse $poss] {
		incr e
		set newfield [string replace $newfield $s $e [string toupper [string index $newfield $e]]]
	}
	regsub -all -- - $newfield _ newfield
	return $newfield
}

proc cg_tomonetdb {args} {
	global scriptname action
	if {[llength $args] != 3} {
		puts stderr "format is: $scriptname $action db table tsvfile"
		exit 1
	}
	foreach {db table tsvfile} $args break

	set tsvfile [file normalize $tsvfile]
	set num [lindex [exec wc -l $tsvfile] 0]
	set f [gzopen $tsvfile]
	set header [tsv_open $f comment]
	set offset [expr {2+[llength [split $comment \n]]}]
	incr num [expr {1-$offset}]
	array set ftype {begin integer end integer start integer type varchar(10)}
	set sql ""
	set fieldtrans {}
	foreach field $header {
		set type [get ftype($field) clob]
		set newfield [monetdb_convfield $field]
		lappend fieldtrans $field $newfield
		lappend sql "\"$newfield\" $type"
	}
	cg_monetdb_sql $db "drop table \"$table\""
	cg_monetdb_sql $db "create table \"$table\" ([join $sql ,\n]);\n"
	# exec {*}[gzcat $tsvfile] $tsvfile | mclient -d$db -s "copy $num offset $offset records into \"$table\" from stdin delimiters '\t', '\n' null as '';"
	set o [open $tsvfile.temp w]
	fconfigure $f -encoding binary -translation binary
	fconfigure $o -encoding binary -translation binary
	fcopy $f $o
	flush $o
	close $f
	close $o
	cg_monetdb_sql $db "copy $num records into \"$table\" from '$tsvfile.temp' delimiters '\t', '\n' null as '';"
	# cg_monetdb_sql $db "select count(*) from \"$table\""
	cg_monetdb_sql $db [subst {alter table "$table" add column "rowid" serial}]
	foreach {key value} [list file $tsvfile time [file mtime $tsvfile] fieldtrans $fieldtrans] {
		if {[catch {
			cg_monetdb_sql $db [subst {insert into "genomecomb_info" ("table","key","value") values('$table','$key','$value')}]
		}]} {
			cg_monetdb_sql $db [subst {update "genomecomb_info" set "value" = '$value' where "table" = '$table' and "key" = '$key'}]
		}
	}
	# cg_monetdb_sql $db [subst {select * from "genomecomb_info" where "table" = '$table'}]
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
