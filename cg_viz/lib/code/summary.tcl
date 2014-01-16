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
		foreach {rowdef coldef celdef} $definition break
		if {$celdef ne ""} {
			lappend coldef $celdef
		}
		if {[file exists $tdata(indexdir)/info.tsv]} {
			set info [infofile_read $tdata(indexdir)/info.tsv]
			catch {set numlines [dict get $info size]}
		}
		if {![info exists numlines]} {
			progress start 1
			progress message "Creating Summary, please be patient (no progress shown)"
			update
			cg select -f $tdata(fields) -q $tdata(query) -g $rowdef -gc $coldef $tdata(file) $tdata(indexdir)/summary_results.tsv
			progress stop
		} else {
			set step [expr {$numlines/10}]
			if {$step > 50000} {set step 50000} elseif {$step < 1} {set step 1}
			progress start [expr {$numlines + 1}] "Making Summary" "Making Summary"
			bgcg [list $object queryprogress] [privatevar $object bgexechandle] \
				select -v $step -f $tdata(fields) -q $tdata(query) -g $rowdef -gc $coldef $tdata(file) $tdata(indexdir)/summary_results.tsv
			progress stop
		}
		if {$file ne ""} {
			file copy $tdata(indexdir)/summary_results.tsv $file
			return
		} else {
			set f [open $tdata(indexdir)/summary_results.tsv]
			set result [csv_file $f \t]
			return $result
		}
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
		set fieldsamples [list {}]
	} else {
		set fieldnames {{}}
		set fieldqueries {{}}
		set fieldsamples {{}}
	}
	set fields [$object fields]
	foreach {field values} $coldef {
		if {$field eq "sample"} {
			if {[llength [list_remove $fieldsamples {}]]} {
				error "error in definition, you can only use field sample once"
			}
			if {![llength $values]} {
				set values [samples $fields]
			}
			set tempqueries {}
			set tempnames {}
			set tempsamples {}
			foreach fieldquery $fieldqueries fieldname $fieldnames fieldsample $fieldsamples {
				foreach value $values {
					lappend tempqueries $fieldquery
					lappend tempnames [list_union $fieldname [list $value]]
					lappend tempsamples $value
				}
			}
			set fieldqueries $tempqueries
			set fieldnames $tempnames
			set fieldsamples $tempsamples
		} else {
			set tempqueries {}
			set tempnames {}
			set tempsamples {}
			foreach fieldquery $fieldqueries fieldname $fieldnames fieldsample $fieldsamples {
				if {$fieldsample eq ""} {
					set mfield [$object mfield $field extrafields colfield]
				} else {
					set mfield [$object mfield ${field}-$fieldsample extrafields colfield]
				}
				# lappend groupby \"$mfield\"
				if {![llength $values]} {
					set values [$object subsql [subst {select distinct("$mfield") from "temp"}] {*}$extrafields]
				}
				foreach value $values {
					lappend tempqueries [list_union $fieldquery [list "\"$mfield\" = \'$value\'"]]
					lappend tempnames [list_union $fieldname [list $value]]
					lappend tempsamples $fieldsample
				}
			}
			set fieldqueries $tempqueries
			set fieldnames $tempnames
			set fieldsamples $tempsamples
		}
	}
	set colnum 0
	set num 0
	foreach fieldquery $fieldqueries fieldname $fieldnames fieldsample $fieldsamples {
		if {$fieldsample eq ""} {
			set mfield [$object mfield $cellfield extrafields colfield]
		} else {
			set mfield [$object mfield ${cellfield}-$fieldsample extrafields colfield]
		}
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
		foreach line [ssort -natural [split $data \n]] {
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
