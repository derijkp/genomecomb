monetdbinit

proc mselect_idtopos {header id fields} {
	set id [string trim $id "\"\' "]
	set poss [list_cor $header $fields]
	if {[inlist $poss -1]} {error "sample $id not present"}
	return [lmath_calc $poss + 1]
}

proc mselect_samegeno {a11 a12 a21 a22} {
	return "(($a11 = $a21 and $a12 = $a22) or ($a11 = $a22 and $a12 = $a21))"
}

proc mselect_compare {header ids trans} {
	if {[llength $ids]  != 2} {error "compare only with 2 ids"}
	foreach {id1 id2} $ids break
	set id1 [string trim $id1 "\"\' "]
	set id2 [string trim $id2 "\"\' "]
	set a11 \"[trans $trans alleleSeq1-$id1]\"
	set a21 \"[trans $trans alleleSeq2-$id1]\"
	set a12 \"[trans $trans alleleSeq1-$id2]\"
	set a22 \"[trans $trans alleleSeq2-$id2]\"
	set s1 [trans $trans sequenced-$id1]
	set s2 [trans $trans sequenced-$id2]
	set sql [subst {(case
		when "$s1" = 'r' and "$s2" = 'r' then 'r'
		when "$s1" = 'u' or "$s2" = 'u' then 'un'
		when "$s1" = 'r' and "$s2" = 'v' then 'df'
		when "$s1" = 'v' and "$s2" = 'r' then 'df'
		when "$s1" = 'v' and "$s2" = 'v' and [mselect_samegeno $a11 $a21 $a12 $a22] then 'sm'
		when "$s1" = 'v' and "$s2" = 'v' then 'mm'
		else '?' end)
	}]
}

proc mselect_zyg {header ids trans} {
	if {[llength $ids]  != 1} {error "wrong # of args for zyg: only one allowed"}
	foreach {id1} $ids break
	set id1 [string trim $id1 "\"\' "]
	set a11 \"[trans $trans alleleSeq1-$id1]\"
	set a21 \"[trans $trans alleleSeq2-$id1]\"
	set s1 [trans $trans sequenced-$id1]
	set sql [subst {(case
		when "$s1" = 'r' then 'r'
		when "$s1" = 'u' then 'u'
		when "$s1" = 'v' and "$a11" = "alt" and "$a21" = "alt" then 'm'
		when "$s1" = 'v' and (("$a11" = "ref" and "$a21" = "alt") or ("$a11" = "alt" and "$a21" = "ref") then 't'
		when "$s1" = 'v' and "$a11" = "alt" or "$a21" = "alt" then 'c'
		else 'o' end)
	}]
}

proc mselect_sm {header ids trans} {
	set id1 [list_pop ids]
	set id1 [string trim $id1 "\"\' "]
	set a11 \"[trans $trans alleleSeq1-$id1]\"
	set a21 \"[trans $trans alleleSeq2-$id1]\"
	set temp [list "(\"[trans $trans sequenced-$id1]\" = \'v\')"]
	foreach id $ids {
		set id [string trim $id "\"\' "]
		set a12 \"[trans $trans alleleSeq1-$id]\"
		set a22 \"[trans $trans alleleSeq2-$id]\"
		lappend temp "(\"[trans $trans sequenced-$id]\" = \'v\')"  "[mselect_samegeno $a11 $a21 $a12 $a22]"
	}
	set temp "([join $temp " and "])"
}

proc mselect_same {header ids trans} {
	set id1 [list_pop ids]
	set id1 [string trim $id1 "\"\' "]
	set a11 \"[trans $trans alleleSeq1-$id1]\"
	set a21 \"[trans $trans alleleSeq2-$id1]\"
	set sequenced \"[trans $trans sequenced-$id1]\"
	set seqlist [list "$sequenced <> \'u\'"]
	set temp {}
	foreach id $ids {
		set id [string trim $id "\"\' "]
		lappend seqlist "$sequenced <> \'u\'"
		set a12 \"[trans $trans alleleSeq1-$id]\"
		set a22 \"[trans $trans alleleSeq2-$id]\"
		set sequenced \"[trans $trans sequenced-$id]\"
		lappend temp "[mselect_samegeno $a11 $a21 $a12 $a22]"
	}
	set temp "([join $seqlist " && "] && [join $temp " && "])"
}

proc mselect_df {header ids trans} {
	set temp1 {}
	set temp2 {}
	set seqlist {}
	foreach id $ids {
		set id [string trim $id "\"\' "]
		set a1 \"[trans $trans alleleSeq1-$id]\"
		set a2 \"[trans $trans alleleSeq2-$id]\"
		set sequenced \"[trans $trans sequenced-$id]\"
		lappend seqlist "($sequenced <> \'u\')"
		lappend temp1 "($sequenced = \'v\')"
		lappend temp2 "($sequenced = \'r\')"
	}
	set temp "(([join $seqlist " and " ]) and ([join $temp1 " or " ]) and ([join $temp2 " or "]))"
}

proc mselect_mm {header ids trans} {
	set temp1 {}
	set temp2 {}
	set list {}
	set seqlist {}
	foreach id $ids {
		set id [string trim $id "\"\' "]
		set a1 \"[trans $trans alleleSeq1-$id]\"
		set a2 \"[trans $trans alleleSeq2-$id]\"
		set sequenced \"[trans $trans sequenced-$id]\"
		lappend seqlist "($sequenced = \'v\')"
		lappend list [list $a1 $a2]
	}
	set ref [lsearch $header reference]
	if {$ref != -1} {	
		set ref "reference"
	} else {
		set ref "ref"
	}
	while {[llength $list]} {
		foreach {a1 a2} [list_pop list] break
		lappend temp1 "(($a1 <> $ref) or ($a2 <> $ref))"
		list_foreach {a12 a22} $list {
			lappend temp2 "not [mselect_samegeno $a1 $a2 $a12 $a22]"
		}
	}
	set temp "(([join $seqlist " and " ]) and ([join $temp1 " and " ]) and ([join $temp2 " or "]))"
}

proc mselect_un {header ids trans} {
	set temp1 {}
	set temp2 {}
	foreach id $ids {
		set id [string trim $id "\"\' "]
		foreach {a1 a2 sequenced} [mselect_idtopos $header $id [list [trans $trans alleleSeq1-$id] [trans $trans alleleSeq2-$id] [trans $trans sequenced-$id]]] break
		lappend temp1 "(\$$sequenced == \"v\")"
		lappend temp2 "(\$$sequenced == \"u\")"
	}
	set temp "(([join $temp2 " || "]) && ([join $temp1 " || " ]))"
}

proc mselect_count {arguments header trans} {
	# We do *1 to force to a correct type (otherwise we get an error in monetdb)
	set test [list_pop arguments]
	set temp {}
	foreach q $arguments {
		lappend q {*}$test
		set q [tsv_select_precedence $q]
		lappend temp "([mselect_detokenize $q $header neededfields $trans])*1"
	}
	return "([join $temp " + "])"
}
# cg select -q 'count($alleleSeq1,$alleleSeq2, == "G") == 1' annotvar.tsv

proc mselect_if {arguments header trans} {
	if {[llength $arguments] != 3} {
		error "wrong # args for function if, must be: if(condition,true,false)"
	}
	foreach field {condition true false} q $arguments {
		set q [tsv_select_precedence $q]
		set $field "([mselect_detokenize $q $header neededfields $trans])"
	}
	return "(case when ($condition) then $true else $false end)"
}

#proc mselect_lmin {} {
#	upvar awkfunctions awkfunctions
#	lappend awkfunctions {
#		function lmin(list,def) {
#			if (def == nill) {def = 999999999}
#		        split(list,a,/[,;]/);
#		        minv = a[1];
#		        for (i in a) {
#		                if (a[i] != a[i]+0) {a[i] = def}
#		                if (a[i] < minv) {minv = a[i]}
#		        }
#		        return minv
#		}
#	}
#}
#
#proc mselect_lmax {} {
#	upvar awkfunctions awkfunctions
#	lappend awkfunctions {
#		function lmax(list,def) {
#			if (def == nill) {def = -999999999}
#		        split(list,a,/[,;]/);
#		        maxv = a[1];
#		        for (i in a) {
#		                if (a[i] != a[i]+0) {a[i] = def}
#		                if (a[i] > maxv) {maxv = a[i]}
#		        }
#		        return maxv
#		}
#	}
#}
#
#proc mselect_min {num} {
#	incr num
#	set temp [list_fill $num 1 1]
#	set result [subst {
#		function min(a[join $temp ,a]) \{
#			if (a1 ~ /def=(.*)/) {
#				def = substr(a1,5)
#				if (a2 != a2 + 0) {a2 = def}
#				minv = a2
#			} else {
#				if (a1 != a1 + 0) {a1 = def}
#				def = 999999999
#				minv = a1
#			}
#	}]
#	foreach field [lrange $temp 1 end] {
#		append result [subst {
#			if (a$field == nill) {return minv}
#			if (a$field != a$field + 0) {a$field = def}
#			if (a$field < minv) {minv = a$field}
#		}]
#	}
#	append result [subst {
#			return minv
#		\}
#	}]
#	return $result
#}
#
#proc mselect_max {num} {
#	incr num
#	set temp [list_fill $num 1 1]
#	set result [subst {
#		function max(a[join $temp ,a]) \{
#			if (a1 ~ /def=(.*)/) {
#				def = substr(a1,5)
#				if (a2 != a2 + 0) {a2 = def}
#				maxv = a2
#			} else {
#				if (a1 != a1 + 0) {a1 = def}
#				def = -999999999
#				maxv = a1
#			}
#	}]
#	foreach field [lrange $temp 1 end] {
#		append result [subst {
#			if (a$field == nill) {return maxv}
#			if (a$field != a$field + 0) {a$field = def}
#			if (a$field > maxv) {maxv = a$field}
#		}]
#	}
#	append result [subst {
#			return maxv
#		\}
#	}]
#	return $result
#}
#
#proc mselect_counthasone {ids} {
#	upvar awkfunctions awkfunctions
#	upvar tsv_funcnum tsv_funcnum
#	set test [list_pop ids]
#	lappend awkfunctions [subst -nocommands {
#		function tsvfunc${tsv_funcnum}(list) {
#		        split(list,a,/[,;]/);
#		        for (i in a) {
#		                if (a[i] $test) {return 1}
#		        }
#		        return 0
#		}
#	}]
#	set temp {}
#	foreach id $ids {
#		lappend temp "tsvfunc${tsv_funcnum}($id)"
#	}
#	return "([join $temp " + "])"
#}
#
#proc mselect_counthasall {ids} {
#	upvar awkfunctions awkfunctions
#	upvar tsv_funcnum tsv_funcnum
#	set test [list_pop ids]
#	lappend awkfunctions [subst -nocommands {
#		function tsvfunc${tsv_funcnum}(list) {
#		        split(list,a,/[,;]/);
#		        for (i in a) {
#		                if (!(a[i] $test)) {return 0}
#		        }
#		        return 1
#		}
#	}]
#	set temp {}
#	foreach id $ids {
#		lappend temp "tsvfunc${tsv_funcnum}($id)"
#	}
#	return "([join $temp " + "])"
#}

proc mselect_min {ids header} {
putsvars ids header
	set result "(case \n"
	set else [lindex $ids end]
	set test [lrange $ids 0 end-1]
	set pos 1
	foreach id $test {
		set rest [list_sub $ids -exclude $pos]
		incr pos
		set temp {}
		foreach rid $rest {
			lappend temp "($id) <= ($rid)"
		}
		append result "when [join $temp " and "] then ($id)\n"
	}
	append result "else $else \n end)"
	return $result
}

proc mselect_max {ids header} {
	set result "(case \n"
	set else [lindex $ids end]
	set test [lrange $ids 0 end-1]
	set pos 1
	foreach id $test {
		set rest [list_sub $ids -exclude $pos]
		incr pos
		set temp {}
		foreach rid $rest {
			lappend temp "($id) >= ($rid)"
		}
		append result "when [join $temp " and "] then ($id)\n"
	}
	append result "else $else \n end)"
	return $result
}

proc mselect_asint {ids header} {
	return "cast ([lindex $ids 0] as integer)"
}

proc mselect_asdouble {ids header} {
	return "cast ([lindex $ids 0] as double)"
}

proc mselect_region {ids header} {
	set ids [split $ids ":-_ ,\{\}\(\)\"'"]
	set ids [list_remove $ids {}]
	set poss [tsv_basicfields $header 3]
	set fields [list_sub $header $poss]
	foreach {rchr rbegin rend} $fields break
	lappend neededfields {*}$fields
	set result {}
	foreach {chr begin end} $ids {
		set q [list "(\"$rchr\" = '$chr' \
			or (left(\"$rchr\",3) = 'chr' and right(\"$rchr\",length(\"$rchr\")-3) = '$chr') \
			or (left('$chr',3) = 'chr' and right('$chr',length('$chr')-3) = \"$rchr\") \
		)"]
		set temp "(case when left(\"$rchr\",3) = 'chr' then right(\"$rchr\",length(\"$rchr\")-3) else \"$rchr\" end)"
		append temp " = (case when left('$chr',3) = 'chr' then right('$chr',length('$chr')-3) else '$chr' end)"
		set q [list $temp]
		if {$begin ne ""} {
			lappend q "\"$rend\" > ($begin) "
		}
		if {$end ne ""} {
			lappend q "\"$rbegin\" < ($end)"
		}
		lappend result "([join $q " and "])"
	}
	return "(([join $result "\) or \("]))"
}

proc mselect_expandfield {header field qpossVar trans} {
	upvar $qpossVar qposs
	set qposs [list_find -glob $header $field]
	if {![llength $qposs]} {
		error "no fields matched \"$field\""
	}
	set result {}
	foreach field [list_sub $header $qposs] {
		if {[dict exists $trans $field]} {set field [dict get $trans $field]}
		lappend result $field
	}
	set qposs [lmath_calc $qposs + 1]
	return $result
}

proc mselect_expandfields {header qfields qcodeVar trans} {
	upvar $qcodeVar qcode
	upvar tsv_funcnum tsv_funcnum
	set qcode {}
	set rfields {}
	foreach field $qfields {
		if {$field eq "ROW"} {
			lappend rfields ROW
			lappend qcode {("rowid"-1)}
			continue
		}
		set pos [string first = $field]
		if {$pos != -1} {
			lappend rfields [string range $field 0 [expr {$pos-1}]]
			set code [string range $field [expr {$pos+1}] end]
			lappend qcode [mselect_expandcode $header $code $trans]
		} elseif {[string first * $field] != -1} {
			lappend rfields {*}[mselect_expandfield $header $field poss $trans]
			foreach pos $poss {
				lappend qcode {}
			}
		} else {
			set pos [lsearch $header $field]
			if {$pos == -1} {
				error "field \"$field\" not present"
			}
			if {[dict exists $trans $field]} {set field [dict get $trans $field]}
			lappend rfields $field
			lappend qcode {}
		}
	}
	return $rfields
}

array set mselect_transop {
	== = != <>
	eq = ne <>
	in in ni {not in}
	&& and || or
	! not
}

proc mselect_detokenize {tokens header neededfieldsVar trans} {
	global mselect_transop
	upvar $neededfieldsVar neededfields
	set result {}
	foreach line $tokens {
		foreach {type val} $line break
		switch -exact $type {
			@num - @val {
				set pre [string index $val 0]
				set post [string index $val end]
				if {($pre eq "\{" && $post eq "\}") || ($pre eq "\"" && $post eq "\"")} {
					set val \'[string range $val 1 end-1]\'
				}
				lappend result $val
			}
			@var {
				if {[dict exists $trans $val]} {set val [dict get $trans $val]}
				lappend result \"$val\"
			}
			@op {
				if {[info exists mselect_transop($val)]} {
					lappend result $mselect_transop($val)
				} else {
					lappend result $val
				}
			}
			@newop {
				lappend result $val
			}
			@function {
				set arguments [lrange $line 2 end]
				set ids {}
				foreach el $arguments {
					lappend ids [mselect_detokenize $el $header neededfields $trans]
				}
				switch $val {
					sm {
						set temp [mselect_sm $header $ids $trans]
					}
					same {
						set temp [mselect_same $header $ids $trans]
					}
					df {
						set temp [mselect_df $header $ids $trans]
					}
					mm {
						set temp [mselect_mm $header $ids $trans]
					}
					un {
						set temp [mselect_un $header $ids $trans]
					}
					compare {
						set temp [mselect_compare $header $ids $trans]
					}
					zyg {
						set temp [mselect_zyg $header $ids $trans]
					}
					hovar {
						set temp [mselect_hovar $header $ids $trans]
					}
					count {
						set temp [mselect_count $arguments $header $trans]
					}
					if {
						set temp [mselect_if $arguments $header $trans]
					}
					region {
						set temp [mselect_region $ids $header]
					}
					min {
						set temp [mselect_min $ids $header]
					}
					max {
						set temp [mselect_max $ids $header]
					}
					asint {
						set temp [mselect_asint $ids $header]
					}
					asdouble {
						set temp [mselect_asdouble $ids $header]
					}
					counthasone {
						foreach {operand value} [lindex $line end] break
						set operand [mselect_detokenize [list $operand] $header needefields $trans]
						set value [mselect_detokenize [list $value] $header needefields $trans]
						set temp [mselect_counthasone [lrange $ids 0 end-1] $operand $value]
					}
					counthasall {
						foreach {operand value} [lindex $line end] break
						set operand [mselect_detokenize [list $operand] $header needefields $trans]
						set value [mselect_detokenize [list $value] $header needefields $trans]
						set temp [mselect_counthasall [lrange $ids 0 end-1] $operand $value]
					}
					hasone {
						if {[llength $ids] == 2} {
							foreach {operand value} [lindex $line end] break
							set operand [mselect_detokenize [list $operand] $header needefields $trans]
							set value [mselect_detokenize [list $value] $header needefields $trans]
							set temp "${val}\([lindex $ids 0], \"$operand\", $value\)"
						} else {
							set temp "${val}\([join $ids ", "]\)"
						}
					}
					oneof {
						set temp [mselect_oneof $ids]
					}
					regexp {
						set temp "pcre_match\([join $ids ", "]\)"
					}
					default {
						set temp "${val}\([join $ids ", "]\)"
					}
				}
				lappend result $temp
			}
			@braces {
				lappend result \([mselect_detokenize [lrange $line 1 end] $header needefields $trans]\)
			}
			@command {
				append result \[$val\]
			}
			default {
				error "unkown token $type"
			}
		}
	}
	return [join $result " "]
}

proc monetdb_makesql {table theader query qfieldsVar {sortfields {}} {inverse 0} {trans {}} {offset {}} {limit {}}} {
	upvar $qfieldsVar qfields
# putsvars table theader query qfieldsVar sortfields inverse trans offset limit qfields
	set sqlfields {}
	set sqlwhere {}
	set sqlsort {}
	if {$qfields eq ""} {
		set sqlfields $theader
	}
	set qfields [mselect_expandfields $theader $qfields qcode $trans]
	if {$inverse} {
		set qfields [list_lremove $theader $qfields]
	}
	foreach field $qfields code $qcode {
		if {$code ne ""} {
			lappend sqlfields "$code as \"$field\""
		} else {
			lappend sqlfields \"$field\"
		}
	}
	set sqlfields [join $sqlfields ,]
	if {[llength $sortfields]} {
		set sqlsort "order by \"[join $sortfields "\",\""]\""
	}
	if {$query ne ""} {
		set tquery [mselect_expandcode $theader $query $trans]
		set sqlwhere "where $tquery"
	}
	set sql "select $sqlfields from \"$table\" $sqlwhere $sqlsort"
	if {[isint $limit]} {
		append sql " limit $limit"
	}
	if {[isint $offset]} {
		append sql " offset $offset"
	}
	return $sql
}

proc monetdb_select {db table query {qfields {}} {sortfields {}} {newheader {}} {sepheader {}} {out stdout} {inverse 0} {group {}} {groupcols {}}} {
#putsvars db table header query qfields sortfields newheader sepheader out hc inverse group groupcols
	set theader [list_remove [cg_monetdb_fields $db $table] rowid]
	if {![catch {
		set trans [cg_monetdb_sql $db [subst {select "value" from "_genomecomb_info" where "table" = '$table' and "key" = 'fieldtrans'}]]
	}]} {
		set header [dict keys $trans]
	} else {
		set trans {}
		set header $theader
	}
	if {$out ne ""} {
		fconfigure $out -buffering none
	}
	set keepheader {}
	if {![inlist $sortfields rowid]} {
		lappend sortfields rowid
	}
	if {$group ne ""} {
		foreach {resultheader result} [monetdb_makesql_group $db $table $header $theader $query qfields $inverse $trans $group $groupcols] break
		set nh ""
		if {$sepheader ne ""} {
			file_write $sepheader [join $resultheader \t]\n
		} else {
			puts $out [join $resultheader \t]
		}
		puts $out $result
	} else {
		set sql [monetdb_makesql $table $theader $query qfields $sortfields $inverse $trans]
		if {$qfields ne ""} {
			set nh $qfields
		} else {
			set nh $header
		}
		if {[llength $newheader]} {
			if {[llength $newheader] != [llength $nh]} {error "new header (-nh) of wrong length for query results"}
			set nh $newheader
		}
		#putslog -------------sql--------------------
		#putslog sql:\n$sql
		#putslog ------------------------------------
		if {$out ne ""} {
			if {$sepheader ne ""} {
				file_write $sepheader ${keepheader}[join $nh \t]\n
			} else {
				puts $out ${keepheader}[join $nh \t]
			}
			exec mclient -d $db -f tab -s $sql >@ $out
		} else {
			exec mclient -d $db -f tab -s $sql
		}
	}
}

proc cg_mselect {args} {
	if {[llength $args] == 0} {
		errorformat mselect
		exit 1
	}
	set query {}; set fields {}; set sortfields {}; set newheader {}; set sepheader ""; set hc 0; set inverse 0; set printheader 0; set group {}; set groupcols {}
	set pos 0
	foreach {key value} $args {
		switch -- $key {
			-q {
				set query $value
			}
			-qf {
				error "not yet"
				set f [gzopen $value]
				set header [tsv_open $f]
				set data [csv_file $f \t]
				gzclose $f
				set query {}
				foreach line $data {
					set el ""
					foreach field $header v $line {
						lappend el "\$$field == \"$v\""
					}
					lappend query "\( [join $el " && "] \)"
				}
				set query [join $query " || "]
			}
			-f {set fields $value}
			-rf {
				set fields $value
				set inverse 1
			}
			-g {set group $value}
			-gc {lappend groupcols $value}
			-nh {set newheader $value}
			-sh {set sepheader $value}
			-hc {set hc 1}
			-hf {set hc $value}
			-s {set sortfields $value}
			-n {
				set header [cg_monetdb_fields $db $table]
				set names {}
				foreach col $header {
					set split [split $col _]
					if {[llength $split] > 1} {
						lappend names [lindex $split end]
					}
				}
				puts stdout [join [list_remdup $names] \n]
				exit 0
			}
			-h {
				incr pos
				set args [lrange $args $pos end]
				set len [llength $args]
				if {$len < 2 || $len > 3} {
					errorformat mselect
					exit 1
				}
				foreach {db table outfile} $args break
				set header [cg_monetdb_fields $db $table]
				puts stdout [join $header \n]
				exit 0
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	if {[llength $groupcols] && ![llength $group]} {
		error "cannot use -gc option without -g option"
	}
	set args [lrange $args $pos end]
	set len [llength $args]
	if {$len < 2 || $len > 3} {
		errorformat mselect
		exit 1
	}
	foreach {db table outfile} $args break
	regsub -all {\n#[^\n]*} $fields {} fields
	regsub -all {\n#[^\n]*} $query {} query
	regsub -all {\n|\t} $query { } query
	set query [string trim $query]
	#puts stderr [list fields=$fields query=$query]
	if {$len == 3} {
		set o [open $outfile w]
	} else {
		set o stdout
	}
	monetdb_select $db $table $query $fields $sortfields $newheader $sepheader $o $hc $inverse $group $groupcols
#	set error [catch {monetdb_select $db $table $query $fields $sortfields $newheader $sepheader $o $hc $inverse $group $groupcols} result]
#	if {$o ne "stdout"} {catch {close $o}}
#	if {$error} {
#		puts stderr $result
#		exit 1
#	}
}

if 0 {
	set dbfarm ~/my-dbfarm
	set database test
	set port 50000

set db test
set table compar
set tsvfile /media/wd2t/complgen/projects/dlb1/dlb_compar.tsv
set tsvfile /complgen/projects/amg_md200/compar/annotamg_md200_compar.tsv


cg monetdb sql test 'select count(*) from compar'
cg monetdb sql test 'select * from compar'
cg monetdb sql test $'select ("begin"/1000000) as rpos,count(*) as count from compar where chromosome = \'chr1\' group by rpos order by rpos'
cg monetdb sql test $'select (min("begin","end")) as "min" from compar where chromosome = \'chr1\' and "begin" < 100000'

select "chromosome","begin","end" from "compar" where "knownGenextype" like 'ex%';

select ("begin"/1000000) as rpos,count(*) as count from compar where chromosome = 'chr1' group by rpos;
select ("begin"/10000) AS rpos,"begin" from compar where chromosome = 'chr1' and "begin" < 100000;
select ("begin"/10000) AS rpos,"begin" from compar where chromosome = 'chr1' and "begin" < 100000 and rpos <= 4;

cg monetdb sql test $'select (case when "begin"<"end" then "begin" else "end" end) as "minimum","begin","end" from compar where chromosome = \'chr1\' and (case when "begin"<"end" then "begin" else "end" end) = 31843'

}
