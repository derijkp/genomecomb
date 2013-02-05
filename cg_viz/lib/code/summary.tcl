table_tsv method summary {definition {file {}}} {
	private $object tdata
	if {![info exists tdata(sqlbackend_db)]} {
		error "No fast summaries without db backend"
	}

if 0 {
	set definition {
		{chromosome {}}
		{sample {} sequenced-*}
		{sequenced-* count}
	}
	set definition {
		{chromosome {} alt {A G T C}}
		{sequenced-testNA19240chr2122cg {v r}}
		{coverage-testNA19240chr2122cg count}
	}
}
	foreach {rowdef coldef celdef} $definition break
	set colnames {}
	set resultfields {}
	set groupby {}
	set where {}
	foreach {field values} $rowdef {
		set mfield [monetdb_convfield $field]
		lappend groupby \"$mfield\"
		lappend resultfields \"$mfield\"
		lappend colnames [list $field]
		if {[llength $values]} {
			lappend where "\"$mfield\" in \(\'[join $values \',\']\'\)"
		}
	}
	set fieldnames {{}}
	set fieldqueries {{}}
	foreach {field values} $coldef {
		set mfield [monetdb_convfield $field]
		# lappend groupby \"$mfield\"
		if {![llength $values]} {
			set values [$object subsql [subst {select distinct("$mfield") from "temp"}]]
		}
		set tempqueries {}
		set tempnames {}
		foreach fieldquery $fieldqueries fieldname $fieldnames {
			foreach value $values {
				lappend tempqueries [list_union $fieldquery [list "\"$mfield\" = \'$value\'"]]
				lappend tempnames [list_union $fieldname [list $value]]
			}
		}
		set fieldqueries $tempqueries
		set fieldnames $tempnames
	}
	foreach {cellfield aggregate} $celdef break
	set colnum 0
	set num 0
	foreach fieldquery $fieldqueries fieldname $fieldnames {
		set mfield [monetdb_convfield $cellfield]
		lappend resultfields "${aggregate}(case when [join $fieldquery " and "] then \"$mfield\" else null end) as \"col$num\""
		lappend colnames $fieldname
		incr num
	}
	set sql [subst {select [join $resultfields " , \n"] from "temp" }]
	if {[llength $where]} {append sql \n [subst {where [join $where " and "]}]}
	append sql \n [subst {group by [join $groupby " , "]}]
	append sql \n [subst {order by [join $groupby " , "]}]
	set data [$object subsql $sql]
	if {$file ne ""} {
		set f [open $file w]
		puts $f [join $colnames \t]
		puts $f $data
		close $f
	} else {
		set result [list $colnames]
		foreach line [lsort -dict [split $data \n]] {
			lappend result [split $line \t]
		}
		return $result
	}

}
