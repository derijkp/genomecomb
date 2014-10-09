proc multidb_monet_createtable {database table def} {
	set fields {}
	set def [split [string trim $def] \n]
	cg_monetdb_sql $database "create table \"$table\" (\n\t[join $def ,\n]\n)"
}

proc multidb_monet_create {database} {
	cg monetdb createdb $database
	# foreach table {geno analysis experiment var} {cg_monetdb_sql $database [subst {drop table "$table"}]}
	multidb_monet_createtable $database var {
		"chromosome" text
		"begin" integer
		"end" integer
		"type" char(3)
		"ref" text
		"alt" text
		"id" integer not null primary key
	}
	multidb_monet_createtable $database experiment {
		"id" text not null primary key
		"gentli_id" integer
		"technology" text
		"provider" text
		"time" time
	}
	multidb_monet_createtable $database analysis {
		"id" integer not null primary key
		"experiment" text
		"ngsseq" text
		"mapper" text
		"varcall" text
		"gentli_ngsseq" integer
		"reference" text
		"individual" text
		"sample" integer
		unique("experiment","ngsseq","mapper","varcall")
	}
	multidb_monet_createtable $database geno {
		"var" integer not null references "var"("id") on update cascade
		"analysis" integer not null references "analysis"("id") on update cascade
		"sequenced" char(1)
		"zyg" char(1)
		"alleleSeq1" text
		"alleleSeq2" text
		"quality" float
		"coverage" integer
		primary key ("var","analysis")
	}
	cg_monetdb_sql $database {create index "var_basicfields" on "var"("chromosome","begin","end","type")}
	cg_monetdb_sql $database {create index "geno_zygosity" on "geno"("zyg")}
	cg_monetdb_sql $database {create index "geno_sequenced" on "geno"("sequenced")}
}

proc multidb_monet_open {compar_dir database} {
	if {[catch {
		cg_monetdb_sql $database {select count(*) from "var"}
	} varcount]} {
		if {[regexp {please create it first} $varcount]} {
			multidb_monet_create $database
			set varcount 0
		} else {
			error $varcount
		}
	}
	if {![file exists $compar_dir/vars.tsv]} {
		file_write $compar_dir/vars.tsv [join {chromosome begin end type ref alt id} \t]\n
		set sql [subst {
			select "chromosome","begin","end","type","ref","alt","id" from "var"
		}]
		exec mclient -d $database -f tab -s $sql >> $compar_dir/vars.tsv
		file delete $compar_dir/vars.tsv.maxid
		file delete $compar_dir/vars.tsv.count
	}
	if {![file exists $compar_dir/vars.tsv.maxid]} {
		set maxid [cg_monetdb_sql $database {select max("id") from "var"}]
		if {$maxid eq ""} {set maxid 0}
		file_write $compar_dir/vars.tsv.maxid $maxid
	}
	if {![file exists $compar_dir/vars.tsv.count]} {
		set count [cg_monetdb_sql $database {select count("id") from "var"}]
		file_write $compar_dir/vars.tsv.count $count
	}
	if {![file exists $compar_dir/analysis.tsv]} {
		set fields [cg_monetdb_fields $database analysis]
		file_write $compar_dir/analysis.tsv [join $fields \t]\n
		set sql [subst {
			select "[join $fields \",\"]" from "analysis"
		}]
		exec mclient -d $database -f tab -s $sql >> $compar_dir/analysis.tsv
		file delete $compar_dir/analysis.tsv.maxid
		file delete $compar_dir/analysis.tsv.count
	}
	if {![file exists $compar_dir/analysis.tsv.maxid]} {
		set maxid [cg_monetdb_sql $database {select max("id") from "analysis"}]
		if {$maxid eq ""} {set maxid 0}
		file_write $compar_dir/analysis.tsv.maxid $maxid
	}
	if {![file exists $compar_dir/analysis.tsv.count]} {
		set count [cg_monetdb_sql $database {select count("id") from "analysis"}]
		file_write $compar_dir/analysis.tsv.count $count
	}
	cg_monetdb_fields $database geno
}

proc multidb_monet_import {compar_dir database} {
	set file $compar_dir/vars.tsv.insert
	set num [lindex [exec wc -l $file] 0]
	incr num -1
	set f [open $file]
	set header [tsv_open $f comment]
	close $f
	set offset [expr {2+[llength [split $comment \n]]}]
	cg_monetdb_sql $database "copy $num offset $offset records into \"var\" from '$compar_dir/vars.tsv.insert' delimiters '\t', '\n' null as '?';"
	file rename -force $compar_dir/vars.tsv.insert $compar_dir/vars.tsv.insert.done
#	cg_monetdb_sql $database "copy $num2 offset $offset2 records into \"geno\" from '$tempfile.sampledata.tsv' delimiters '\t', '\n' null as '?';"

}
