proc monetdbinit {} {
	global monetdbdir
	if {[catch {
		set file [exec which mclient]
	}]} {
		foreach path [list \
			$::appdir/../extern/monetdb \
			$::appdir/../../extern/monetdb \
			/home/peter/apps/monetdb \
		] {
			set file $path/bin/mclient
			if {[file exists $file]} {
				break
			}
		}
	}
	while {1} {
		if {[catch {
			set link [file readlink $file]
		}]} break
		if {$link eq ""} break
		if {[file pathtype $link] eq "absolute"} {
			set file $link
		} else {
			set file [file_absolute [file dir $file]/$link]
		}
	}
	if {![file exists $file]} {
		error "could not find monetdb installation"
	}
	set monetdbdir [file dir [file dir $file]]
	set ::env(PATH) ${::env(PATH)}:$monetdbdir/bin
	if {[info exists ::env(LD_LIBRARY_PATH)]} {
		set ::env(LD_LIBRARY_PATH) ${::env(LD_LIBRARY_PATH)}:$monetdbdir/lib
	} else {
		set ::env(LD_LIBRARY_PATH) $monetdbdir/lib
	}
}
monetdbinit

proc monetdb_getinfo {dbfarm} {
	set dbfarm [file_absolute $dbfarm]
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
	catch {cg_monetdb_sql $database {
		-- case sensitive
		create function pcre_match(s string, pattern string)
		  returns BOOLEAN
		  external name pcre.match;
	}} e
	puts "Added pcre_match: $e"
	catch {cg_monetdb_sql $database {
		-- case insensitive
		create function pcre_imatch(s string, pattern string)
		  returns BOOLEAN
		  external name pcre.imatch;
	}} e
	puts "Added pcre_imatch: $e"
	catch {cg_monetdb_sql $database {
		create table _genomecomb_tableinfo (
			"table" clob not null primary key,
			"file" clob,
			"time" bigint,
			"fieldtrans" clob,
			"size" bigint,
			"type" varchar(10)
		)
	}} e
	puts "Created table _genomecomb_tableinfo: $e"
	catch {cg_monetdb_sql $database {
		create table _genomecomb_info (
			"table" clob,
			"key" clob,
			"value" clob,
			primary key("table","value")
		)
	}} e
	puts "Created table _genomecomb_info: $e"
}

proc cg_monetdb_sql {database sql} {
	set result [exec mclient -d $database -f tab -s $sql]
}

proc cg_monetdb_sql_list {database sql} {
	set result {}
	foreach line [split [exec mclient -d $database -f tab -s $sql] \n] {
		lappend result [split $line \t]
	}
	return $result
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
			set dbfarm [file_absolute $dbfarm]
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
		discover {
			return [exec monetdb discover]
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
		deletedb {
			set database [list_pop args]
			exec monetdb {*}$args stop $database
			exec monetdb {*}$args destroy -f $database
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
		set c [string index $newfield $e]
		set u [string toupper $c]
		if {$u eq $c} {set u x$u}
		set newfield [string replace $newfield $s $e $u]
	}
	regsub -all -- - $newfield _ newfield
	return $newfield
}

proc monetdb_backend {db file} {
	# create database if it does not exist yet
	if {[catch {
		cg_monetdb_sql $db {select 1}
	}]} {
		cg_monetdb createdb $db
	}
	# get table
	cg_monetdb_sql $db {select * from "_genomecomb_tableinfo"}
	set data [lindex [cg_monetdb_sql_list $db [subst {
		select "table","time" from "_genomecomb_tableinfo" where "file" = '$file'
	}]] 0]
	if {[llength $data]} {
		foreach {table time} $data break
	}
	set mtime [file mtime $file]
	if {[info exists table] && $time < $mtime} {
		putslog "$file is more recently updated than database $db, removing and reloading ..."
		catch {cg_monetdb_sql $db [subst {drop table "$table"}]} e
		puts $e
	}
	if {![info exists table] || $time < $mtime} {
		putslog "Inserting $file into database $db"
		set tableroot [file root [file tail $file]]
		regsub -all {[^A-Za-z0-9]} $tableroot {} tableroot
		set tables [cg_monetdb_tables $db]
		set table $tableroot
		set num 2
		while {$table in $tables} {
			set table ${tableroot}_$num
			incr num
		}
		cg_tomonetdb $db $table $file
		return $table
	} else {
		return $table
	}
}

proc cg_tomonetdb {args} {
	global scriptname action
	if {[llength $args] != 3} {
		error "format is: $scriptname $action db table tsvfile"
	}
	foreach {db table tsvfile} $args break

	set tsvfile [file_absolute $tsvfile]
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
	catch {cg_monetdb_sql $db "drop table \"$table\""}
	cg_monetdb_sql $db "create table \"$table\" ([join $sql ,\n]);\n"
	# exec {*}[gzcat $tsvfile] $tsvfile | mclient -d$db -s "copy $num offset $offset records into \"$table\" from stdin delimiters '\t', '\n' null as '';"
	set temptsvfile [file_tempwrite $tsvfile]
	set o [open $temptsvfile w]
	fconfigure $f -encoding binary -translation binary
	fconfigure $o -encoding binary -translation binary
	fcopy $f $o
	flush $o
	close $f
	close $o
	cg_monetdb_sql $db "copy $num records into \"$table\" from '$temptsvfile' delimiters '\t', '\n' null as '';"
	# cg_monetdb_sql $db "select count(*) from \"$table\""
	cg_monetdb_sql $db [subst {alter table "$table" add column "rowid" serial}]
	set time [file mtime $tsvfile]
	if {[catch {
		cg_monetdb_sql $db [subst {
			insert into "_genomecomb_tableinfo" ("table","file","time","fieldtrans")
			values('$table','$tsvfile','$time','$fieldtrans')
		}]
	}]} {
		cg_monetdb_sql $db [subst {
			update "_genomecomb_tableinfo"
			set "file" = '$tsvfile', "time" = '$tinme', "fieldtrans" = '$fieldtrans'
			where "table" = '$table'
		}]
	}
	file delete $temptsvfile
}

proc cg_genomecombinfo {cmd table args} {
	switch $cmd {
		get {
			set key [lindex $args 0]
			cg_monetdb_sql $db [subst {select "value" from "_genomecomb_info" where "table" = '$table' and "key" = '$key'}]
		}
		set {
			foreach {key value} $args {
				catch {
					cg_monetdb_sql $db [subst {delete from "_genomecomb_info" where "table" = '$table' and "key" = '$key'}]
				}
				cg_monetdb_sql $db [subst {insert into "_genomecomb_info" ("table","key","value") values('$table','$key','$value')}]
			}
		}
	}
}

#proc cg_tomonetdbnorm {args} {
#	global scriptname action
#	if {[llength $args] != 3} {
#		error "format is: $scriptname $action db table tsvfile"
#	}
#	foreach {db table tsvfile} $args break
#
#	unset -nocomplain fieldsa
#	unset -nocomplain posa
#	unset -nocomplain afieldsa
#	unset -nocomplain aposa
#	set aposa(post) {}
#	set aposa(pre) {}
#	set afieldsa(post) {}
#	set afieldsa(pre) {}
#	set samples {}
#	set samplecols {}
#	set fpos pre
#	set colpos 0
#	foreach col $header {
#		set pos [string first - $col]
#		if {$pos != -1} {
#			set field [string range $col 0 [expr {$pos-1}]]
#			list_addnew samplecols $field
#			incr pos
#			set sample [string range $col $pos end]
#			lappend samples $sample
#			lappend fieldsa($sample) $field
#			lappend posa($sample) $colpos
#			set fpos post
#		} else {
#			lappend afieldsa($fpos) $col
#			lappend aposa($fpos) $colpos
#		}
#		incr colpos
#	}
#	set samples [list_remdup $samples]
#	set samplecols [list_union [list_common {sequenced zyg quality totalScore1 totalScore2 alleleSeq1 alleleSeq2} $samplecols] $samplecols]
#	set todo {}
#	foreach sample $samples {
#		set cor [list_cor $fieldsa($sample) $samplecols]
#		set cor [list_sub $posa($sample) $cor]
#		set cor [list_change $cor {{} -1}]
#		lappend todo [list_concat $aposa(pre) $cor $aposa(post)]
#	}
#	set common [list_common $samplecols $afieldsa(pre)]
#	if {[llength $common]} {
#		foreach field $common {
#			set newfield ${field}_global
#			set num 1
#			while {[inlist $samplecols $newfield]} {
#				set newfield ${field}_global$num
#				incr num
#			}
#			set afieldsa(pre) [list_change $afieldsa(pre) [list $field $newfield]]
#		}
#	}
#	set common [list_common $samplecols $afieldsa(post)]
#	if {[llength $common]} {
#		foreach field $common {
#			set newfield ${field}_global
#			set num 1
#			while {[inlist $samplecols $newfield]} {
#				set newfield ${field}_global$num
#				incr num
#			}
#			set afieldsa(post) [list_change $afieldsa(post) [list $field $newfield]]
#		}
#	}
#	puts $o [join [list_concat sample $afieldsa(pre) $samplecols $afieldsa(post)] \t]
#	while {![eof $f]} {
#		set line [split [gets $f] \t]
#		if {![llength $line]} continue
#		foreach sample $samples cor $todo {
#			puts $o $sample\t[join [list_sub $line $cor] \t]
#		}
#	}
#
#
#
#
#
#	set tsvfile [file_absolute $tsvfile]
#	set num [lindex [exec wc -l $tsvfile] 0]
#	set f [gzopen $tsvfile]
#	set header [tsv_open $f comment]
#	set offset [expr {2+[llength [split $comment \n]]}]
#	incr num [expr {1-$offset}]
#	array set ftype {begin integer end integer start integer type varchar(10)}
#	set sql ""
#	set fieldtrans {}
#	foreach field $header {
#		set type [get ftype($field) clob]
#		set newfield [monetdb_convfield $field]
#		lappend fieldtrans $field $newfield
#		lappend sql "\"$newfield\" $type"
#	}
#	catch {cg_monetdb_sql $db "drop table \"$table\""}
#	cg_monetdb_sql $db "create table \"$table\" ([join $sql ,\n]);\n"
#	# exec {*}[gzcat $tsvfile] $tsvfile | mclient -d$db -s "copy $num offset $offset records into \"$table\" from stdin delimiters '\t', '\n' null as '';"
#	set o [open $tsvfile.temp w]
#	fconfigure $f -encoding binary -translation binary
#	fconfigure $o -encoding binary -translation binary
#	fcopy $f $o
#	flush $o
#	close $f
#	close $o
#	cg_monetdb_sql $db "copy $num records into \"$table\" from '$tsvfile.temp' delimiters '\t', '\n' null as '';"
#	# cg_monetdb_sql $db "select count(*) from \"$table\""
#	cg_monetdb_sql $db [subst {alter table "$table" add column "rowid" serial}]
#	foreach {key value} [list file $tsvfile time [file mtime $tsvfile] fieldtrans $fieldtrans] {
#		catch {
#			cg_monetdb_sql $db [subst {delete from "_genomecomb_info" where "table" = '$table' and "key" = '$key'}]
#		}
#		cg_monetdb_sql $db [subst {insert into "_genomecomb_info" ("table","key","value") values('$table','$key','$value')}]
#	}
#	# cg_monetdb_sql $db [subst {select * from "_genomecomb_info" where "table" = '$table'}]
#	file delete $tsvfile.temp
#}

if 0 {
	set dbfarm ~/my-dbfarm
	set database test
	set port 50000

set db test
set table compar
set tsvfile /media/wd2t/complgen/projects/dlb1/dlb_compar.tsv
set tsvfile /complgen/projects/amg_md200/compar/annotamg_md200_compar.tsv
}
