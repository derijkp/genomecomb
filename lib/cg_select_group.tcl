proc select_parse_grouptypes {grouptypelist} {
	set typetodoa {max max min min count {} percent total gpercent gtotal avg {avg} stddev {avg m2} stdev {avg m2} distinct distinct list list sum sum}
	set grouptypes {}
	foreach grouptype $grouptypelist {
		if {$grouptype eq "count"} {
			lappend grouptypes count {}
		} elseif {$grouptype eq "percent"} {
			lappend grouptypes percent {}
		} elseif {$grouptype eq "gpercent"} {
			lappend grouptypes gpercent {}
		} elseif {[regexp {^([^()]+)\(([^()]+)\)$} $grouptype temp func field]} {
			if {![dict exists $typetodoa $func]} {
				error "aggregate function $func unknown"
			}
			lappend grouptypes $func $field
		} else {
			error "aggregate function (last element in groupcol argument) must be of the form function(value) or count"
		}
	}
	return $grouptypes
}

proc select_parse_for_samples {group groupcol header} {
	set gsamples {}
	foreach {field values} [list_concat $group $groupcol] {
		if {$field eq "sample"} {
			if {[llength $values]} {
				if {[string first * $values] == -1} {
					set gsamples $values
				} else {
					set gsamples {} 
					foreach pattern $values {
						lappend gsamples {*}[samples $header $pattern]
					}
					set gsamples [list_remdup $gsamples]
				}
			} else {
				set gsamples [samples $header]
			}
			break
		}
	}
	# if no sample field found, do not actually use samples by making gsamples a list with only one empty sample
	if {![llength $gsamples]} {
		set gsamples {{}}
	}
	return $gsamples
}

proc tsv_select_makecol {name code {arg @neededfield@}} {
	subst -nocommands {
		proc $name {$arg} {
			if {[catch {expr {$code}} e]} {
				switch \$e {
					{domain error: argument not in valid range} {return NaN}
					{divide by zero} {return NaN}
				}
			}
			return \$e
		}
	}
}

proc tsv_select_sampleusefield {header field sample calccolsVar {neededfieldsVar {}}} {
	if {$neededfieldsVar ne ""} {upvar $neededfieldsVar neededfields}
	upvar $calccolsVar calccols
	if {[info exists calccols($field)]} {
		set fieldused $field
	} elseif {$sample ne "" && [info exists calccols($field-$sample)]} {
		set fieldused ${field}-$sample
	} elseif {[inlist $header $field]} {
		set fieldused $field
		lappend neededfields $fieldused
	} elseif {$sample ne "" && [inlist $header ${field}-$sample]} {
		set fieldused ${field}-$sample
		lappend neededfields $fieldused
	} else {
		return ""
	}
	return $fieldused
}

proc tsv_select_addaggregatecalc {todolist} {
	set colactions {}
	# add calculations for everything needed for aggregates to colactions
	append colactions \t\t\t\t\t {set resultdatacols($_colname) 1} \n
	append colactions \t\t\t\t\t {incr resultcount($_groupname,$_colname)} \n
	foreach {item todo} $todolist {
		foreach {field fieldused} $item break
		append colactions \t\t\t\t\t [string_change {set _val {@val@}} [list @val@ $field]] \n
		if {[inlist $todo total]} {
			append colactions \t\t\t\t\t {incr resultdata($_colname,t)} \n
		}
		if {[inlist $todo gtotal]} {
			append colactions \t\t\t\t\t {incr resultdata($_groupname,gt)} \n
		}
		if {[inlist $todo max]} {
			append colactions [string_change {
					if {![info exists resultdata($_groupname,$_colname,$_val,max)] || ${@val@} > $resultdata($_groupname,$_colname,$_val,max)} {
						set resultdata($_groupname,$_colname,$_val,max) ${@val@}
					}
			} [list @val@ $fieldused]]
		}
		if {[inlist $todo min]} {
			append colactions [string_change {
					if {![info exists resultdata($_groupname,$_colname,$_val,min)] || ${@val@} < $resultdata($_groupname,$_colname,$_val,min)} {
						set resultdata($_groupname,$_colname,$_val,min) ${@val@}
					}
			} [list @val@ $fieldused]]
		}
		if {[inlist $todo avg]} {
			append colactions [string_change {
					if {![info exists resultdata($_groupname,$_colname,$_val,avg)]} {
						set resultdata($_groupname,$_colname,$_val,avg) 0.0
					}
					set _delta [expr {${@val@}-$resultdata($_groupname,$_colname,$_val,avg)}]
					set resultdata($_groupname,$_colname,$_val,avg) [expr {$resultdata($_groupname,$_colname,$_val,avg) + $_delta/$resultcount($_groupname,$_colname)}]
			} [list @val@ $fieldused]]
		}
		if {[inlist $todo m2]} {
			append colactions [string_change {
					if {![info exists resultdata($_groupname,$_colname,$_val,m2)]} {
						set resultdata($_groupname,$_colname,$_val,m2) 0.0
					}
					set resultdata($_groupname,$_colname,$_val,m2) [expr {$resultdata($_groupname,$_colname,$_val,m2) + $delta*(${@val@}-$resultdata($_groupname,$_colname,$_val,avg))}]
			} [list @val@ $fieldused]]
		}
		if {[inlist $todo distinct]} {
			append colactions [string_change {
					dict set resultdata($_groupname,$_colname,$_val,d) ${@val@} 1
			} [list @val@ $fieldused]]
		}
		if {[inlist $todo list]} {
			append colactions [string_change {
					lappend resultdata($_groupname,$_colname,$_val,d) ${@val@}
			} [list @val@ $fieldused]]
		}
		if {[inlist $todo sum]} {
			append colactions [string_change {
					if {[isdouble ${@val@}]} {
						if {![info exists resultdata($_groupname,$_colname,$_val,s)]} {
							set resultdata($_groupname,$_colname,$_val,s) ${@val@}
						} else {
							set resultdata($_groupname,$_colname,$_val,s) [expr {$resultdata($_groupname,$_colname,$_val,s) + ${@val@}}]
						}
					}
			} [list @val@ $fieldused]]
		}
	}
	#	set stddev [expr {sqrt($m2/($n - 1))}]
	#	set stddev [expr {sqrt($m2/$n)}]
	return $colactions
}

proc tsv_select_addaggregateresult {grouptypes header sample calccolsVar} {
	upvar $calccolsVar calccols
	set calcresults ""
	foreach {func field} $grouptypes {
		if {$func eq "count"} {
			append calcresults {
				lappend result [get resultcount($_groupname,$col) 0]
			}
		} elseif {$func eq "percent"} {
			append calcresults {
				lappend result [format %.2f [expr {100.0*[get resultcount($_groupname,$col) 0]/[get resultdata($col,t) 1]}]]
			}
		} elseif {$func eq "gpercent"} {
			append calcresults {
				lappend result [format %.2f [expr {100.0*[get resultcount($_groupname,$col) 0]/[get resultdata($_groupname,gt) 1]}]]
			}
		} elseif {$func eq "max"} {
			append calcresults [string_change {
				lappend result [get resultdata($_groupname,$col,@field@,max) ""]
			} [list @field@ $field]]
		} elseif {$func eq "min"} {
			append calcresults [string_change {
				lappend result [get resultdata($_groupname,$col,@field@,min) ""]
			} [list @field@ $field]]
		} elseif {$func eq "avg"} {
			append calcresults [string_change {
				lappend result [get resultdata($_groupname,$col,@field@,avg) ""]
			} [list @field@ $field]]
		} elseif {$func eq "stddev" || $func eq "stdev"} {
			append calcresults [string_change {
				if {[info exists resultdata($_groupname,$col,@field@,m2)} {
					lappend result [expr {$resultdata($_groupname,$col,@field@,m2)/$resultcount($_groupname,$col)}]
				} else {
					lappend result ""
				}
			} [list @field@ $field]]
		} elseif {$func eq "distinct"} {
			append calcresults [string_change {
				lappend result [join [ssort -natural [dict keys [get resultdata($_groupname,$col,@field@,d) ""]]] ,]
			} [list @field@ $field]]
		} elseif {$func eq "list"} {
			append calcresults [string_change {
				lappend result [join [get resultdata($_groupname,$col,@field@,d) ""] ,]
			} [list @field@ $field]]
		} elseif {$func eq "sum"} {
			append calcresults [string_change {
				lappend result [get resultdata($_groupname,$col,@field@,s) ""]
			} [list @field@ $field]]
		}
	}
	return $calcresults
}

proc tsv_select_group {header pquery qposs qfields group groupcols neededfields {verbose 0}} {
# putsvars header pquery qposs qfields group groupcols neededfields
	# outcols not used in group
	# start making code
	# neededfields will contain all fields needed by the query proc and calculated fields
	# calculated fields will be calculated in the begining of the query proc
	# these are gathered in precalc
	# precalc is run for every match (sets some variables used in query, etc.)
	regsub -all \n [string trim $group] { } group
	if {[llength $group] == 1} {lappend group {}}
	set typetodoa {max max min min count {} percent total gpercent gtotal avg {avg} stddev {avg m2} stdev {avg m2} distinct distinct list list sum sum}
	unset -nocomplain calccols
	# more than one groupcol not supported (yet)
	set groupcol [lindex $groupcols 0]
	regsub -all \n [string trim $groupcol] { } groupcol
	if {![llength $groupcol]} {
		set groupcol count
	}
	set grouptypelist [split [list_pop groupcol] ,]
	# check groupcol for presence of sample field
	set gsamples [select_parse_for_samples $group $groupcol $header]
	# parse grouptypes (aggregate results), and see which functions are needed
	set grouptypes [select_parse_grouptypes $grouptypelist]
	# check for calculated fields in group, groupcol and grouptypes, add to qposs and qfields for making precalc
	set curpos 0
	set newgroup {}
	foreach {el values} $group {
		set pos [string first = $el]
		if {$pos != -1} {
			set field [string range $el 0 [expr {$pos-1}]]
			set code [string range $el [expr {$pos+1}] end]
			lset group $curpos $field
			lappend qposs [list code $code]
			lappend qfields $field
			lappend newgroup $field $values
		} else {
			lappend newgroup $el $values
		}
		incr curpos
	}
	set group $newgroup
	set curpos 0
	foreach {el values} $groupcol {
		set pos [string first = $el]
		if {$pos != -1} {
			set field [string range $el 0 [expr {$pos-1}]]
			set code [string range $el [expr {$pos+1}] end]
			lset groupcol $curpos $field
			lappend qposs [list code $code]
			lappend qfields $field
		}
		incr curpos 2
	}
	set curpos 1
	foreach {func el} $grouptypes {
		set pos [string first = $el]
		if {$pos != -1} {
			set field [string range $el 0 [expr {$pos-1}]]
			set code [string range $el [expr {$pos+1}] end]
			lset grouptypes $curpos $field
			lappend qposs [list code $code]
			lappend qfields $field
		}
		incr curpos 2
	}
	# make precalc, add functions for calculated columns
	set precalc {}
	set num 0
	set tclcode "package require genomecomb\n"
	foreach el $qposs field $qfields {
		if {![isint $el] && $el ne "" && ![info exists calccols($field)]} {
			set code [tsv_select_expandcode $header [lindex $el 1] tempneededfields]
			set tempneededfields [list_remdup $tempneededfields]
			lappend neededfields {*}$tempneededfields
			lappend outcols make_col$num
			append tclcode [tsv_select_makecol make_col$num $code $tempneededfields]
			lappend precalc "\t\t\t\tset \{$field\} \[make_col$num \$\{[join $tempneededfields \}\ \$\{]\}\]"
			set calccols($field) "make_col$num"
		}
		incr num
	}
	# make colactions, which will be executed only when the colquery and rowquery is true
	set addcols {}
	foreach sample $gsamples {
		# also make query for skipping data for which no cols will be made (colquery)
		set groupname {}
		set col {}
		set rowquery {}
		set colquery {}
		set colactions {}
		# calculate groupname
		foreach {field filter} $group {
			if {$field eq "sample"} {
				lappend groupname $sample
				continue
			} else {
				set fieldused [tsv_select_sampleusefield $header $field $sample calccols neededfields]
			}
			if {$fieldused ne ""} {
				lappend groupname \$\{$fieldused\}
				if {[llength $filter]} {
					if {[string first * $filter] == -1} {
						lappend rowquery "\[inlist \{$filter\} \$\{$fieldused\}\]"
					} else {
						foreach temp $filter {
							lappend rowquery "\[string match \{$temp\} \$\{$fieldused\}\]"
						}
					}
				}
			} else {
				lappend groupname $field
			}
		}
#		append colactions \t\t\t\t\t "set _groupname \[list [join $groupname { }]\]" \n
		append colactions \t\t\t\t\t "set _groupname \"[join $groupname \t]\"" \n
		# first see what I need to calculate requested aggregates
		unset -nocomplain todoa
		foreach {func field} $grouptypes {
			foreach item [dict get $typetodoa $func] {
				if {[inlist {percent gpercent} $func]} {
					set fieldused {}
					list_addnew todoa([list $field $fieldused]) $item
				} else {
					set fieldused [tsv_select_sampleusefield $header $field $sample calccols neededfields]
					if {$fieldused ne ""} {
						list_addnew todoa([list $field $fieldused]) $item
					}
				}
			}
		}
		set todolist [array get todoa]
		# make col and colquery variable
		foreach {field filter} $groupcol {
			if {$field eq "sample"} {
				lappend col $sample
				continue
			} else {
				set fieldused [tsv_select_sampleusefield $header $field $sample calccols neededfields]
			}
			if {$fieldused ne ""} {
				lappend col \$\{$fieldused\}
				if {[llength $filter]} {
					if {[string first * $filter] == -1} {
						lappend colquery "\[inlist \{$filter\} \$\{$fieldused\}\]"
					} else {
						foreach temp $filter {
							lappend colquery "\[string match \{$temp\} \$\{$fieldused\}\]"
						}
					}
				}
			}
		}
		append colactions \t\t\t\t\t {set resultgroups($_groupname) 1} \n
		append colactions \t\t\t\t\t "set _colname \"[join $col -]\"\n"
		# add calculations for everything needed for aggregates to colactions
		append colactions [tsv_select_addaggregatecalc $todolist]
		# create calcresults, that calculate the final agregate results for each group (runs after looping through the file)
		set calcresults [tsv_select_addaggregateresult $grouptypes $header $sample calccols]
		set q [list_concat $rowquery $colquery]
		if {[llength $q]} {
			append addcols \t\t\t\t "if \{[join $q { && }]\} \{\n[string trimright $colactions]\n\t\t\t\t\}\n"
		} else {
			append addcols $colactions
		}
	}
	set neededfields [list_remdup $neededfields]
	set neededcols [list_cor $header $neededfields]
	append tclcode [string_change {
		proc tsv_selectc_query {@neededfields@} {
			global resultdata resultgroups resultcount resultdatacols
			if {@pquery@} {
@precalc@
@addcols@
			}
			return 0
		}
		tsv_selectc tsv_selectc_query {@neededcols@} {} @verbose@
		set cols [ssort -natural [array names resultdatacols]]
		set header [list @grouph@]
		foreach col $cols {
			foreach {func val} @grouptypes@ {
				if {$val ne ""} {append func _$val} 
				lappend header [join [list {*}${col} $func] -]
			}
		}
		puts [join $header \t]
		foreach _groupname [ssort -natural [array names resultgroups]] {
			set result {}
			# puts "$name\t$resultdata($name)"
			foreach col $cols {
				@calcresults@
			}
			puts "$_groupname\t[join $result \t]"
		}
		exit
	} [list @neededfields@ $neededfields @pquery@ $pquery \
		@precalc@ [join [list_remdup $precalc] \n] @addcols@ $addcols \
		@neededcols@ $neededcols @calcresults@ $calcresults \
		@grouptypes@ [list $grouptypes] @grouph@ [list_unmerge $group] @verbose@ $verbose]
	]]
	return $tclcode
}
