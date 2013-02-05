table_tsv method mfield {field extrafieldsVar colfieldVar} {
	upvar $extrafieldsVar extrafields
	upvar $colfieldVar colfield
	set posis [string first = $field]
	if {$posis == -1} {
		set mfield [monetdb_convfield $field]
		set colfield $field
	} else {
		set mfield [string range $field 0 [expr {$posis-1}]]
		lappend extrafields $field
		set colfield $mfield
	}
	return $mfield
}

table_tsv method summary {definition {file {}}} {
	private $object tdata
	if {![info exists tdata(sqlbackend_db)]} {
		error "No fast summaries without db backend"
	}
	foreach {rowdef coldef celdef} $definition break
	foreach {cellfield aggregate} $celdef break
	set colnames {}
	set resultfields {}
	set groupby {}
	set where {}
	set extrafields {}
	foreach {field values} $rowdef {
		set mfield [$object mfield $field extrafields colfield]
		lappend groupby \"$mfield\"
		lappend resultfields \"$mfield\"
		lappend colnames [list $colfield]
		if {[llength $values]} {
			lappend where "\"$mfield\" in \(\'[join $values \',\']\'\)"
		}
	}
	if {![llength $coldef]} {
		set fieldnames all
		set fieldqueries [list {}]
	} else {
		set fieldnames {{}}
		set fieldqueries {{}}
	}
	foreach {field values} $coldef {
		set mfield [$object mfield $field extrafields colfield]
		# lappend groupby \"$mfield\"
		if {![llength $values]} {
			set values [$object subsql [subst {select distinct("$mfield") from "temp"}] {*}$extrafields]
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
	set mfield [$object mfield $cellfield extrafields colfield]
	set colnum 0
	set num 0
	foreach fieldquery $fieldqueries fieldname $fieldnames {
		if {$fieldquery eq ""} {
			set temp \"$mfield\"
		} else {
			set temp "case when [join $fieldquery " and "] then \"$mfield\" else null end"
		}
		if {$aggregate eq "avg"} {
			lappend resultfields "round(avg($temp),2) as \"col$num\""
		} else {
			lappend resultfields "${aggregate}($temp) as \"col$num\""
		}
		lappend colnames $fieldname
		incr num
	}
	set sql [subst {select [join $resultfields " , \n"] from "temp" }]
	if {[llength $where]} {append sql \n [subst {where [join $where " and "]}]}
	if {[llength $groupby]} {
		append sql \n [subst {group by [join $groupby " , "]}]
		append sql \n [subst {order by [join $groupby " , "]}]
	}
	set data [$object subsql $sql {*}$extrafields]
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
	set definition {
		{coverage-testNA19240chr2122cg {}}
		{{comp=compare(testNA19239chr2122cg,testNA19240chr2122cg)} {}}
		{coverage-testNA19240chr2122cg count}
	}
	foreach {view(summary_rows) view(summary_cols) view(summary_cells)} $definition break
}
