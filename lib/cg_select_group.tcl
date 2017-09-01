set typetodoa {max max min min count {} percent total gpercent gtotal avg {avg} stddev {avg m2} stdev {avg m2} distinct distinct ucount distinct list list sum sum median list q1 list q3 list}

proc select_parse_grouptypes {grouptypelist} {
	global typetodoa
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

proc select_parse_for_samples {group groupcol header {calccolsVar {}}} {
	if {$calccolsVar ne ""} {upvar $calccolsVar calccols}
	set gsamples {}
	set usesamples 0
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
			set usesamples 1
			break
		}
		if {![inlist $header $field] && ![info exists calccols($field)]} {
			if {[lsearch -glob $header ${field}-*] != -1
				|| [llength [array names calccols ${field}-*]]
				|| (![catch {tsv_select_sampleinfo_wildcard ${field}-* $header} temp] && [llength $temp])
			} {
				# hidden sample
				set usesamples 1
			}
		}
	}
	if {$usesamples && ![llength $gsamples]} {set gsamples [samples $header]}
	# if no sample field found, do not actually use samples by making gsamples a list with only one empty sample
	if {![llength $gsamples]} {
		set gsamples {{}}
	}
	return $gsamples
}

proc tsv_select_makecol {name code {arg @neededfield@} {prequery {}}} {
	subst -nocommands {
		proc $name {$arg} {
			$prequery
			if {[catch {expr {$code}} e]} {
				switch -regexp \$e {
					{domain error: argument not in valid range} {return NaN}
					{divide by zero} {return NaN}
					{invalid command name "tcl::mathfunc::.*"} {
						regexp {tcl::mathfunc::([^"]*)} \$e temp temp
						error "unknown function \$temp"
					}
					default {error \$e}
				}
			}
			return \$e
		}
	}
}

proc tsv_select_sampleusefield {header field sample calccolsVar {neededfieldsVar {}}} {
	if {$neededfieldsVar ne ""} {upvar $neededfieldsVar neededfields}
	upvar $calccolsVar calccols
	if {$sample ne "" && [info exists calccols($field-$sample)]} {
		set fieldused ${field}-$sample
	} elseif {$sample ne "" && [inlist $header ${field}-$sample]} {
		set fieldused ${field}-$sample
		lappend neededfields $fieldused
	} elseif {[info exists calccols($field)]} {
		set fieldused $field
	} elseif {[inlist $header $field]} {
		set fieldused $field
		lappend neededfields $fieldused
	} elseif {[inlist $header sample] && ![catch {tsv_select_sampleinfo $field $header}]} {
		# long format sampleinfo
		set calccols($field) [list "\t\t\t\tset \{$field\} \[tsv_select_sampleinfo_long \"$field\" \$sample\]\n"]
		lappend neededfields sample
		set fieldused $field
	} else {
		return ""
	}
	return $fieldused
}

proc tsv_select_addforeach {varVar typesVar code} {
	upvar $varVar loopsa
	upvar $typesVar loopstypea
	set loops [array names loopsa]
	if {![llength $loops]} {
		return $code
	}
	set pre {} ; set post {}
	set pre {}
	foreach loop $loops {
		set values $loopsa($loop)
		if {$values eq ""} {
			if {$loopstypea($loop) == 1} {
				append pre "foreach _$loop \[split \$\{$loop\} {;,}\] \{\n"
			} else {
				append pre "foreach _$loop \[list_remdup \[split \$\{$loop\} {;,}\]\] \{\n"
			}
		} else {
			# this is for when values is known beforehand, e.g. coming from sampleinfo
			if {$loopstypea($loop) == 1} {
				append pre "foreach _$loop [list $values] \{\n"
			} else {
				append pre "foreach _$loop [list [list_remdup $values]] \{\n"
			}
		}
		append post \}\n
	}
	return $pre$code$post
}

proc tsv_select_matchfilter {filter value} {
	if {$filter in {{} *}} {return 1}
	if {[string first * $filter] == -1} {
		return [inlist $filter $value]
	} else {
		foreach temp $filter {
			if {[string match $temp $value]} {return 1}
		}
	}
	return 0
}

proc tsv_select_applyfilter {filter values} {
	if {$filter in {{} *}} {
		return $values
	} else {
		set poss [list_find -regexp $filter {\*}]
		set patterns [list_sub $filter $poss]
		set dofilter [list_sub $filter -exclude $poss]
		set list {}
		foreach el [split $values {;, }] {
			if {[inlist $dofilter $el]} {
				lappend list $el
				continue
			}
			foreach pattern $patterns {
				if {[string match $pattern $el]} {
					lappend list $el
					continue
				}
			}
		}
		return $list
	}
}

proc tsv_select_addaggregatecalc {todolist} {
	set colactions {}
	# add calculations for everything needed for aggregates to colactions
	append colactions \t\t\t\t\t {set resultdatacols($_colname) 1} \n
	append colactions \t\t\t\t\t {incr resultcount($_groupname,$_colname)} \n
	foreach {item todo} $todolist {
		foreach {field fieldused valuetype} $item break
		append colactions \t\t\t\t\t [string_change {set _val {@val@}} [list @val@ $field]] \n
		if {[inlist $todo total]} {
			append colactions \t\t\t\t\t {incr resultdata($_colname,t)} \n
		}
		if {[inlist $todo gtotal]} {
			append colactions \t\t\t\t\t {incr resultdata($_groupname,gt)} \n
		}
		if {[inlist $todo max]} {
			if {$valuetype} {
				if {[isdouble $fieldused]} {
					append colactions [string_change {
						if {![info exists resultdata($_groupname,$_colname,$_val,max)] || @val@ > $resultdata($_groupname,$_colname,$_val,max)} {
							set resultdata($_groupname,$_colname,$_val,max) {@val@}
						}
					} [list @val@ $fieldused]]
				}
			} else {
				append colactions [string_change {
					if {[isdouble ${@val@}] && (![info exists resultdata($_groupname,$_colname,$_val,max)] || ${@val@} > $resultdata($_groupname,$_colname,$_val,max))} {
						set resultdata($_groupname,$_colname,$_val,max) ${@val@}
					}
				} [list @val@ $fieldused]]
			}
		}
		if {[inlist $todo min]} {
			if {$valuetype} {
				if {[isdouble $fieldused]} {
					append colactions [string_change {
						if {![info exists resultdata($_groupname,$_colname,$_val,min)] || @val@ < $resultdata($_groupname,$_colname,$_val,min)} {
							set resultdata($_groupname,$_colname,$_val,min) {@val@}
						}
					} [list @val@ $fieldused]]
				}
			} else {
				append colactions [string_change {
					if {[isdouble ${@val@}] && (![info exists resultdata($_groupname,$_colname,$_val,min)] || ${@val@} < $resultdata($_groupname,$_colname,$_val,min))} {
						set resultdata($_groupname,$_colname,$_val,min) ${@val@}
					}
				} [list @val@ $fieldused]]
			}
		}
		if {[inlist $todo avg]} {
			append colactions {
					if {![info exists resultdata($_groupname,$_colname,$_val,avg)]} {
						set resultdata($_groupname,$_colname,$_val,avg) 0.0
					}
			}
			if {$valuetype} {
				if {[isdouble $fieldused]} {
					append colactions [string_change {
						set _delta [expr {@val@-$resultdata($_groupname,$_colname,$_val,avg)}]
						set resultdata($_groupname,$_colname,$_val,avg) [expr {$resultdata($_groupname,$_colname,$_val,avg) + $_delta/$resultcount($_groupname,$_colname)}]
					} [list @val@ $fieldused]]
				}
			} else {
				append colactions [string_change {
					if {[string is double ${@val@}]} {
						set _delta [expr {${@val@}-$resultdata($_groupname,$_colname,$_val,avg)}]
						set resultdata($_groupname,$_colname,$_val,avg) [expr {$resultdata($_groupname,$_colname,$_val,avg) + $_delta/$resultcount($_groupname,$_colname)}]
					}
				} [list @val@ $fieldused]]
			}
		}
		if {[inlist $todo m2]} {
			append colactions {
					if {![info exists resultdata($_groupname,$_colname,$_val,m2)]} {
						set resultdata($_groupname,$_colname,$_val,m2) 0.0
					}
			}
			if {$valuetype} {
				if {[isdouble $fieldused]} {
					append colactions [string_change {
						set resultdata($_groupname,$_colname,$_val,m2) [expr {$resultdata($_groupname,$_colname,$_val,m2) + $delta*(@val@-$resultdata($_groupname,$_colname,$_val,avg))}]
					} [list @val@ $fieldused]]
				}
			} else {
				append colactions [string_change {
					if {[string is double ${@val@}]} {
						set resultdata($_groupname,$_colname,$_val,m2) [expr {$resultdata($_groupname,$_colname,$_val,m2) + $delta*(${@val@}-$resultdata($_groupname,$_colname,$_val,avg))}]
					}
				} [list @val@ $fieldused]]
			}
		}
		if {[inlist $todo distinct]} {
			if {!$valuetype} {
				set fieldused \$\{$fieldused\}
			} elseif {$fieldused eq ""} {
				set fieldused {{}}
			}
			append colactions [string_change {
					dict set resultdata($_groupname,$_colname,$_val,d) @val@ 1
			} [list @val@ $fieldused]]
		}
		if {[inlist $todo list]} {
			if {!$valuetype} {
				set fieldused \$\{$fieldused\}
			} elseif {$fieldused eq ""} {
				set fieldused {{}}
			}
			append colactions [string_change {
					lappend resultdata($_groupname,$_colname,$_val,d) @val@
			} [list @val@ $fieldused]]
		}
		if {[inlist $todo sum]} {
			if {$valuetype} {
				if {[isdouble $fieldused]} {
					append colactions [string_change {
						if {![info exists resultdata($_groupname,$_colname,$_val,s)]} {
							set resultdata($_groupname,$_colname,$_val,s) {@val@}
						} else {
							set resultdata($_groupname,$_colname,$_val,s) [expr {$resultdata($_groupname,$_colname,$_val,s) + @val@}]
						}
					} [list @val@ $fieldused]]
				}
			} else {
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
		} elseif {$func eq "ucount"} {
			append calcresults [string_change {
				lappend result [join [llength [dict keys [get resultdata($_groupname,$col,@field@,d) ""]]] ,]
			} [list @field@ $field]]
		} elseif {$func eq "list"} {
			append calcresults [string_change {
				lappend result [join [get resultdata($_groupname,$col,@field@,d) ""] ,]
			} [list @field@ $field]]
		} elseif {$func eq "median"} {
			append calcresults [string_change {
				lappend result [tcl::mathfunc::median {*}[get resultdata($_groupname,$col,@field@,d) ""]]
			} [list @field@ $field]]
		} elseif {$func eq "q1"} {
			append calcresults [string_change {
				lappend result [tcl::mathfunc::q1 {*}[get resultdata($_groupname,$col,@field@,d) ""]]
			} [list @field@ $field]]
		} elseif {$func eq "q3"} {
			append calcresults [string_change {
				lappend result [tcl::mathfunc::q3 {*}[get resultdata($_groupname,$col,@field@,d) ""]]
			} [list @field@ $field]]
		} elseif {$func eq "sum"} {
			append calcresults [string_change {
				lappend result [get resultdata($_groupname,$col,@field@,s) ""]
			} [list @field@ $field]]
		}
	}
	return $calcresults
}

proc tsv_select_combinations {list} {
	set temp [lindex $list 0]
	foreach values [lrange $list 1 end] {
		set newtemp {}
		foreach el $temp {
			foreach value $values {
				lappend newtemp [list {*}$el $value]
			}
		}
		set temp $newtemp
	}
	return $temp
}

proc tsv_select_group {header pquery qposs qfields group groupcols neededfields sortfields} {
# putsvars header pquery qposs qfields group groupcols neededfields sortfields
	global typetodoa
	# outcols not used in group
	# start making code
	# neededfields will contain all fields needed by the query proc and calculated fields
	# calculated fields will be calculated in the begining of the query proc
	# these are gathered in precalc
	# precalc is run for every match (sets some variables used in query, etc.)
	regsub -all \n [string trim $group] { } group
	if {[llength $group] == 1} {lappend group {}}
	# start making code
	set tclcode "package require genomecomb\n"
	unset -nocomplain calccols
	# more than one groupcol not supported (yet)
	set groupcol [lindex $groupcols 0]
	regsub -all \n [string trim $groupcol] { } groupcol
	if {![llength $groupcol]} {
		set groupcol count
	}
	set grouptypelist [split [list_pop groupcol] ,]
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
	# add functions needed for calculated columns to calccols, will be put later into precalc (or prequery)
	set num 0
	foreach el $qposs field $qfields {
		if {![isint $el] && $el ne "" && ![info exists calccols($field)]} {
			if {[lindex $el 0] eq "code"} {
				set tempneededfields {}
				set code [tsv_select_expandcode $header [lindex $el 1] tempneededfields prequery calccols]
				set tempneededfields [list_remdup $tempneededfields]
				lappend neededfields {*}$tempneededfields
				# code for query
				append tclcode [tsv_select_makecol make_col$num $code $tempneededfields $prequery]
				if {[llength $tempneededfields]} {
					set calccols($field) [list "\t\t\t\tset \{$field\} \[make_col$num \$\{[join $tempneededfields \}\ \$\{]\}\]\n" $tempneededfields]
				} else {
					set calccols($field) [list "\t\t\t\tset \{$field\} \[make_col$num\]\n" $tempneededfields]
				}
			} elseif {[lindex $el 0] eq "directcode"} {
				set calccols($field) [list "\t\t\t\tset \{$field\} \[[lindex $el 1]\]\n"]
			}
		}
		incr num
	}
	if {[inlist $header sample]} {
		# if the file contains a sample field, treat it as any other field
		# gsamples is list of one (empty) element to go over the setup once, without incorporating sample specific code
		set gsamples {{}}
		set sampleinheader 1
	} else {
		# if no sample field, we have a wide format with sample info in field-sample format
		# check group and groupcol for presence of sample field to handle the query using it
		set gsamples [select_parse_for_samples $group $groupcol $header calccols]
		set sampleinheader 0
	}
	# make colactions, which will be executed only when the colquery and rowquery is true
	set addcols {}
	set oneok 0
	set fieldsnotfound {}
	unset -nocomplain loopsa
	unset -nocomplain loopstypea
	foreach sample $gsamples {
		# also make query for skipping data for which no cols will be made (colquery)
		set groupname {}
		set col {}
		set rowquery {}
		set colquery {}
		set colactions {}
		# calculate groupname
		unset -nocomplain sloopsa
		unset -nocomplain sloopstypea
		set found 1
		foreach {field filter} $group {
			if {!$sampleinheader && $field eq "sample"} {
				if {$sample ne ""} {
					lappend groupname $sample
				} else {
					lappend groupname \$sample
					lappend neededfields sample
				}
				continue
			}
			if {[string index $field 0] eq "+"} {
				set field [string range $field 1 end]
				set loop 1
			} elseif {[string index $field 0] eq "-" && [string length $field] > 1} {
				set field [string range $field 1 end]
				set loop 2
			} else {
				set loop 0
			}
			set fieldused [tsv_select_sampleusefield $header $field $sample calccols neededfields]
			if {$fieldused ne ""} {
				if {$loop} {
					set sloopsa($fieldused) {}
					set sloopstypea($fieldused) $loop
					set fieldused _$fieldused
				}
				lappend groupname \$\{$fieldused\}
				if {[llength $filter]} {
					if {[string first * $filter] == -1} {
						lappend rowquery "\[inlist \{$filter\} \$\{$fieldused\}\]"
					} else {
						set temprowquery {}
						foreach temp $filter {
							lappend temprowquery "\[string match \{$temp\} \$\{$fieldused\}\]"
						}
						lappend rowquery [join $temprowquery " || "]
					}
				}
			} elseif {![catch {tsv_select_sampleinfo $field-$sample $header} value] || ![catch {tsv_select_sampleinfo $field $header} value]} {
				if {$loop} {
					set list [tsv_select_applyfilter $filter [split $value {;,}]]
					set sloopsa($field) [list_concat [get sloopsa($field) ""] $list]
					set sloopstypea($field) $loop
					lappend groupname \$\{_$field\}
				} else {
					lappend groupname $value
					if {![tsv_select_matchfilter $filter $value]} {lappend rowquery 0}
				}
			} elseif {$field in {all - _}} {
				lappend groupname $field
			} else {
				lappend fieldsnotfound $field
				set found 0
				# error "field \"$field\" not present in file (or sampleinfo)"
			}
		}
		if {!$found} continue
		set oneok 1
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
						list_addnew todoa([list $field $fieldused 0]) $item
					} elseif {![catch {tsv_select_sampleinfo $field-$sample $header} value] || ![catch {tsv_select_sampleinfo $field $header} value]} {
						list_addnew todoa([list $field $value 1]) $item
					}
				}
			}
		}
		set todolist [array get todoa]
		# make col and colquery variable
		foreach {field filter} $groupcol {
			if {$field eq "sample"} {
				if {$sample ne ""} {
					lappend col $sample
				} else {
					lappend col \$sample
					lappend neededfields sample
				}
				continue
			}
			if {[string index $field 0] eq "+"} {
				set field [string range $field 1 end]
				set loop 1
			} elseif {[string index $field 0] eq "-"} {
				set field [string range $field 1 end]
				set loop 2
			} else {
				set loop 0
			}
			set fieldused [tsv_select_sampleusefield $header $field $sample calccols neededfields]
			if {$fieldused ne ""} {
				if {$loop} {
					set sloopsa($fieldused) {}
					set sloopstypea($fieldused) $loop
					set fieldused _$fieldused
				}
				lappend col \$\{$fieldused\}
				if {[llength $filter]} {
					if {[string first * $filter] == -1} {
						lappend colquery "\[inlist \{$filter\} \$\{$fieldused\}\]"
					} else {
						set tempcolquery {}
						foreach temp $filter {
							lappend tempcolquery "\[string match \{$temp\} \$\{$fieldused\}\]"
						}
						lappend colquery [join $colquery " || "]
					}
				}
			} elseif {![catch {tsv_select_sampleinfo $field-$sample $header} value] || ![catch {tsv_select_sampleinfo $field $header} value]} {
				if {$loop} {
					set list [tsv_select_applyfilter $filter [split $value {;,}]]
					set sloopsa($field) [list_concat [get sloopsa($field) ""] $list]
					set sloopstypea($field) $loop
					lappend col \$\{_$field\}
				} else {
					lappend col $value
					if {![tsv_select_matchfilter $filter $value]} {lappend colquery 0}
				}
			}
		}
		append colactions \t\t\t\t\t {set resultgroups($_groupname) 1} \n
		append colactions \t\t\t\t\t "set _colname \"[join $col -]\"\n"
		# add calculations for everything needed for aggregates to colactions
		append colactions [tsv_select_addaggregatecalc $todolist]
		# create calcresults, that calculate the final agregate results for each group (runs after looping through the file)
		set calcresults [tsv_select_addaggregateresult $grouptypes $header $sample calccols]
		set q [list_remove [list_concat $rowquery $colquery] {}]
		if {[llength $q]} {
			set colactions "\t\t\t\tif \{[join $q { && }]\} \{\n[string trimright $colactions]\n\t\t\t\t\}\n"
		}
		append addcols [tsv_select_addforeach sloopsa sloopstypea $colactions]
	}
	if {!$oneok} {
		error "some fields ([join [list_remdup $fieldsnotfound] ,]) needed were not found (in any of the samples)"
	}
	set neededfields [list_remdup $neededfields]
	if {[get ::tsv_select_sampleinfo_islong 0]} {
		append tclcode "[list set ::tsv_select_sampleinfofile $::tsv_select_sampleinfofile]\n"
	}
	# see what we need of calculated fields
	set calcfieldsquery [list_lremove $neededfields $header]
	set calcfieldsquery [list_remove $calcfieldsquery ROW]
	set neededfields [list_lremove $neededfields $calcfieldsquery]
	set prequery {}
	foreach field $calcfieldsquery {
		if {[info exists calccols($field)]} {
			append prequery [lindex $calccols($field) 0]
			lappend neededfields {*}[lindex $calccols($field) 1]
			incr num
		} else {
			# tsv_select_sampleinfo gives not present error if field also not found in sampleinfo
			set value [tsv_select_sampleinfo $field $header]
			if {[string first - $field] != -1} {
				append prequery "\t\t\tset \{$field\} \"$value\"\n"
			} elseif {[inlist $header sample]} {
				append prequery "\t\t\tset \{$field\} \[tsv_select_sampleinfo_long \"$field\" \$sample\]\n"
				lappend neededfields sample
			} else {
				error "field \"$field\" not present"
			}
		}
	}
	# actually make precalc: all calculated fields needed for group, but not for query (these are already in prequery)
	set precalc {}
	foreach field [array names calccols] {
		append precalc [lindex $calccols($field) 0]
	}
	set neededcols [list_cor $header $neededfields]
	set addcols [tsv_select_addforeach loopsa sloopstypea $addcols]
	# pregroup and precol
	set pregroup {}
	foreach {el values} $group {
		if {[info exists pregroup]} {
			if {$el eq "all"} {
				lappend pregroup all
			} elseif {$el eq "sample"} {
				if {$gsamples eq {{}}} {
					if {[llength [set temp [list_sub $values -exclude [list_find -regexp $values {\*}]]]]} {
						lappend pregroup $temp
					} else {
						unset pregroup; break
					}
				} else {
					lappend pregroup $gsamples
				}
			} elseif {[llength [set temp [list_sub $values -exclude [list_find -regexp $values {\*}]]]]} {
				lappend pregroup $temp
			} else {
				unset pregroup; break
			}
		}
	}
	if {[info exists pregroup]} {
		foreach _group [tsv_select_combinations $pregroup] {
			append tclcode "[list set resultgroups([join $_group \t]) 1]\n"
		}
	}
	set precol {}
	foreach {el values} $groupcol {
		if {[info exists precol]} {
			if {$el eq "all"} {
				lappend precol all
			} elseif {$el eq "sample"} {
				if {$gsamples eq {{}}} {
					if {[llength [set temp [list_sub $values -exclude [list_find -regexp $values {\*}]]]]} {
						lappend precol $temp
					} else {
						unset precol; break
					}
				} else {
					lappend precol $gsamples
				}
			} elseif {[llength [set temp [list_sub $values -exclude [list_find -regexp $values {\*}]]]]} {
				lappend precol $temp
			} else {
				unset precol; break
			}
		}
	}
	if {[info exists precol]} {
		foreach _cols [tsv_select_combinations $precol] {
			append tclcode "[list set resultdatacols([join $_cols -]) 1]\n"
		}
	}
	# add rest to tclcode (tsv_selectc_query, calculate aggregates)
	append tclcode {
		proc tsv_selectc_query {@neededfields@} {
			global resultdata resultgroups resultcount resultdatacols
			@prequery@
			if {@pquery@} {
@precalc@
@addcols@
			}
			return 0
		}
		tsv_selectc tsv_selectc_query {@neededcols@} {} @verbose@
		set cols [ssort -natural [array names resultdatacols]]
		if {![llength $cols]} {set cols {{}}}
		set header [list @grouph@]
		foreach col $cols {
			foreach {func val} @grouptypes@ {
				if {$val ne ""} {append func _$val} 
				lappend header [join [list $func {*}${col}] -]
			}
		}
		if {![llength @sortfields@]} {
			set o stdout
		} else {
			set tempfile [tempfile]
			set o [open $tempfile w]
		}
		puts $o [join $header \t]
		foreach _groupname [ssort -natural [array names resultgroups]] {
			set result {}
			# puts $o "$name\t$resultdata($name)"
			foreach col $cols {
				@calcresults@
			}
			puts $o "$_groupname\t[join $result \t]"
		}
		if {[llength @sortfields@]} {
			close $o
			cg select -s @sortfields@ $tempfile >@ stdout
		}
		exit
	}
	set groupheader {}
	foreach el [list_unmerge $group] {
		regsub {^[+-]} $el {} el
		lappend groupheader $el
	}
	set tclcode [string_change $tclcode [list @neededfields@ $neededfields @neededfieldsvals@ \$\{[join $neededfields \}\ \$\{]\} \
		@pquery@ $pquery @prequery@ $prequery\
		@precalc@ $precalc @addcols@ $addcols \
		@neededcols@ $neededcols @calcresults@ $calcresults \
		@grouptypes@ [list $grouptypes] @grouph@ $groupheader @verbose@ [get ::verbose 0] \
		@sortfields@ [list $sortfields]]
	]]
	# file_write /tmp/temp.txt $tclcode\n
	return $tclcode
}
