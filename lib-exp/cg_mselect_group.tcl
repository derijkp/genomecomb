proc mselect_expandcode {header code trans} {
	# variable preprocessor first, to expand *
	# check and exchange variables needed
	set tokens [tsv_select_tokenize $header $code neededfields]
	# detokenize, making necessary changes
	mselect_detokenize $tokens $header neededfields $trans
}

proc monetdb_makesql_group_subsql {db supersql query} {
	set fsql "with \"temp\" as ($supersql) $query"
	exec mclient -d $db -f tab -s $fsql
}

proc monetdb_makesql_mfield {qfields field sample} {
	if {$sample eq ""} {
		set mfield [monetdb_convfield $field]
	} else {
		set mfield [monetdb_convfield $field-$sample]
		set pos [lsearch $qfields $mfield]
		if {$pos == -1} {
			set mfield [monetdb_convfield $field]
		}
	}
	return $mfield
}

proc monetdb_makesql_group {db table header theader query qfieldsVar {inverse 0} {trans {}} {group {}} {groupcols {}}} {
	upvar $qfieldsVar qfields
	set keepqfields $qfields
	# basic query
	set qfields [list_union $qfields [cg_monetdb_fields $db $table]]
	set supersql [monetdb_makesql $table $theader $query qfields {} $inverse $trans]
	# more than one groupcol not supported (yet)
	if {![llength $groupcols]} {
		set groupcol count
	} else {
		set groupcol [lindex $groupcols 0]
	}
	set grouptypelist [split [list_pop groupcol] ,]
	# check groupcols for presence of sample field
	set gsamples [select_parse_for_samples $groupcol $header]
	# parse grouptypes (aggregate results), and see which functions are needed
	set grouptypes [select_parse_grouptypes $grouptypelist]
	# make groupby
	set colnames {}
	set rowname {}
	set groupby {}
	set extrafields {}
	foreach {field} $group {
		set mfield [monetdb_convfield $field]
		lappend groupby \"$mfield\"
		lappend rowname \"$mfield\"
		# lappend colnames [list $colfield]
	}
	set num 0
	set resultfields [list "[join $rowname " + \'_\' + "] as \"col$num\""]
	foreach sample $gsamples {
		# expand columns (fieldnames will contain the name of the column in the result)
		if {![llength $groupcol]} {
			set fieldnames all
			set fieldqueries [list {}]
		} else {
			set fieldnames {{}}
			set fieldqueries {{}}
			foreach {field values} $groupcol {
				if {$field eq "sample"} {
					set tempqueries {}
					set tempnames {}
					foreach fieldquery $fieldqueries fieldname $fieldnames {
						lappend tempqueries $fieldquery
						lappend tempnames [list_union $fieldname [list $sample]]
					}
					set fieldqueries $tempqueries
					set fieldnames $tempnames
				} else {
					set tempqueries {}
					set tempnames {}
					foreach fieldquery $fieldqueries fieldname $fieldnames {
						set mfield [monetdb_makesql_mfield $qfields $field $sample]
						if {![llength $values]} {
							set values [monetdb_makesql_group_subsql $db $supersql [subst {select distinct("$mfield") from "temp"}]]
						}
						foreach value $values {
							lappend tempqueries [list_union $fieldquery [list "\"$mfield\" = \'$value\'"]]
							lappend tempnames [list_union $fieldname [list $value]]
						}
					}
					set fieldqueries $tempqueries
					set fieldnames $tempnames
				}
			}
		}
		foreach fieldquery $fieldqueries fieldname $fieldnames {
			foreach {func field} $grouptypes {
				if {$field eq ""} {
					set mfield 1
					lappend fieldname ${func}
				} else {
					set mfield \"[monetdb_makesql_mfield $qfields $field $sample]\"
					lappend fieldname ${func}_$field
				}
				if {$fieldquery eq ""} {
					set temp $mfield
				} else {
					set temp "case when [join $fieldquery " and "] then $mfield else null end"
				}
				if {$func eq "avg"} {
					lappend resultfields "round(avg($temp),3) as \"col$num\""
				} elseif {$func eq "percent"} {
					lappend resultfields "count($temp) as \"col$num\""
				} else {
					lappend resultfields "${func}($temp) as \"col$num\""
				}
				lappend colnames $fieldname
				incr num
			}
		}
	}
	set sql [subst {select [join $resultfields " , \n"] from "temp" }]
	if {[llength $groupby]} {
		append sql \n [subst {group by [join $groupby " , "]}]
		append sql \n [subst {order by [join $groupby " , "]}]
	}
	set data [monetdb_makesql_group_subsql $db $supersql $sql]
	return [list $colnames $data]
}

