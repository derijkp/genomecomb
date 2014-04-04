#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc tsv_select_compare {ids header neededfieldsVar} {
	upvar $neededfieldsVar neededfields
	set id1 [string trim [list_pop ids] "\"\' "]
	set fields [list sequenced-$id1 alleleSeq1-$id1 alleleSeq2-$id1]
	foreach id $ids {
		set id [string trim $id "\"\' "]
		lappend fields sequenced-$id alleleSeq1-$id alleleSeq2-$id
	}
	lappend needed {*}$fields
	set needed [list_remdup $needed]
	set poss [list_find [list_cor $header $needed] -1]
	if {[llength $poss]} {
		error "Could not find fields needed for compare(\"$id1\",[join $ids ,]): [list_sub $needed $poss]"
	}
	lappend neededfields {*}$needed
	set temp "\[compare \$\{[join $fields "\} \$\{"]\}\]"
}

proc tsv_select_zyg {ids header neededfieldsVar} {
	upvar $neededfieldsVar neededfields
	set id1 [string trim [list_pop ids] "\"\' "]
	if {$id1 eq ""} {
		set needed [list sequenced alleleSeq1 alleleSeq2 ref alt]
	} else {
		set needed [list sequenced-$id1 alleleSeq1-$id1 alleleSeq2-$id1 ref alt]
	}
	set poss [list_cor $header $needed]
	if {[lindex $poss 0] == -1} {
		set needed [lrange $needed 1 end]
		set poss [lrange $poss 1 end]
	}
	set poss [list_find $poss -1]
	if {[llength $poss]} {
		error "Could not find fields needed for zyg(\"$id1\",[join $ids ,]): [list_sub $needed $poss]"
	}
	lappend neededfields {*}$needed
	set temp "\[zyg \$\{[join $needed "\} \$\{"]\}\]"
}

proc tsv_select_sm {ids header neededfieldsVar} {
	upvar $neededfieldsVar neededfields
	set id1 [string trim [list_pop ids] "\"\' "]
	set temp [list "(\$\{sequenced-$id1\} == \"v\")"]
	lappend needed sequenced-$id1 alleleSeq1-$id1 alleleSeq2-$id1
	foreach id $ids {
		set id [string trim $id "\"\' "]
		lappend needed sequenced-$id alleleSeq1-$id alleleSeq2-$id
		lappend temp "(\$\{sequenced-$id\} == \"v\")"  "samegeno(\$\{alleleSeq1-$id1\},\$\{alleleSeq2-$id1\},\$\{alleleSeq1-$id\},\$\{alleleSeq2-$id\})"
	}
	set needed [list_remdup $needed]
	set poss [list_find [list_cor $header $needed] -1]
	if {[llength $poss]} {
		error "Could not find fields needed for sm(\"$id1\",[join $ids ,]): [list_sub $needed $poss]"
	}
	lappend neededfields {*}$needed
	set temp "([join $temp " && "])"
}

proc tsv_select_same {ids header neededfieldsVar} {
	upvar $neededfieldsVar neededfields
	set id1 [string trim [list_pop ids] "\"\' "]
	set seqlist [list "\$\{sequenced-$id1\} != \"u\""]
	lappend needed sequenced-$id1 alleleSeq1-$id1 alleleSeq2-$id1
	set temp {}
	foreach id $ids {
		set id [string trim $id "\"\' "]
		lappend needed sequenced-$id alleleSeq1-$id alleleSeq2-$id
		lappend seqlist "\$\{sequenced-$id\} != \"u\""
		lappend temp "samegeno(\$\{alleleSeq1-$id1\},\$\{alleleSeq2-$id1\},\$\{alleleSeq1-$id\},\$\{alleleSeq2-$id\})"
	}
	set needed [list_remdup $needed]
	set poss [list_find [list_cor $header $needed] -1]
	if {[llength $poss]} {
		error "Could not find fields needed for same(\"$id1\",[join $ids ,]): [list_sub $needed $poss]"
	}
	lappend neededfields {*}$needed
	set temp "([join $seqlist " && "] && [join $temp " && "])"
}

proc tsv_select_df {ids header neededfieldsVar} {
	upvar $neededfieldsVar neededfields
	set temp1 {}
	set temp2 {}
	set seqlist {}
	foreach id $ids {
		set id [string trim $id "\"\' "]
		lappend needed sequenced-$id
		lappend seqlist "(\$\{sequenced-$id\} != \"u\")"
		lappend temp1 "(\$\{sequenced-$id\} == \"v\")"
		lappend temp2 "(\$\{sequenced-$id\} == \"r\")"
	}
	set needed [list_remdup $needed]
	set poss [list_find [list_cor $header $needed] -1]
	if {[llength $poss]} {
		error "Could not find fields needed for df([join $ids ,]): [list_sub $needed $poss]"
	}
	lappend neededfields {*}$needed
	set temp "(([join $seqlist " && " ]) && ([join $temp1 " || " ]) && ([join $temp2 " || "]))"
}

proc tsv_select_mm {ids header neededfieldsVar} {
	upvar $neededfieldsVar neededfields
	set temp1 {}
	set temp2 {}
	set list {}
	set seqlist {}
	foreach id $ids {
		set id [string trim $id "\"\' "]
		lappend seqlist "(\$\{sequenced-$id\} == \"v\")"
		lappend list [list \{alleleSeq1-$id\} \{alleleSeq2-$id\}]
		lappend needed sequenced-$id alleleSeq1-$id alleleSeq2-$id
	}
	if {[lsearch $header reference] != -1} {
		set ref reference
	} else {
		set ref ref
	}
	lappend needed $ref
	while {[llength $list]} {
		foreach {a1 a2} [list_pop list] break
		lappend temp1 "((\$$a1 != \$$ref) || (\$$a2 != \$$ref))"
		list_foreach {a12 a22} $list {
			lappend temp2 "!samegeno(\$$a1,\$$a2,\$$a12,\$$a22)"
		}
	}
	set needed [list_remdup $needed]
	set poss [list_find [list_cor $header $needed] -1]
	if {[llength $poss]} {
		error "Could not find fields needed for mm([join $ids ,]): [list_sub $needed $poss]"
	}
	lappend neededfields {*}$needed
	set temp "(([join $seqlist " && " ]) && ([join $temp1 " && " ]) && ([join $temp2 " || "]))"
}

proc tsv_select_un {ids header neededfieldsVar} {
	upvar $neededfieldsVar neededfields
	set temp1 {}
	set temp2 {}
	foreach id $ids {
		set id [string trim $id "\"\' "]
		# foreach {a1 a2 sequenced} [tsv_select_idtopos $header $id [list alleleSeq1-$id alleleSeq2-$id sequenced-$id]] break
		lappend temp1 "(\$\{sequenced-$id\} == \"v\")"
		lappend temp2 "(\$\{sequenced-$id\} == \"u\")"
		lappend needed sequenced-$id
	}
	set needed [list_remdup $needed]
	set poss [list_find [list_cor $header $needed] -1]
	if {[llength $poss]} {
		error "Could not find fields needed for un([join $ids ,]): [list_sub $needed $poss]"
	}
	lappend neededfields {*}$needed
	set temp "(([join $temp2 " || "]) && ([join $temp1 " || " ]))"
}

proc tsv_select_hovar {ids header neededfieldsVar} {
	upvar $neededfieldsVar neededfields
	set temp {}
	foreach id $ids {
		set id [string trim $id "\"\' "]
		lappend needed sequenced-$id alleleSeq1-$id alleleSeq2-$id
		lappend temp "\$\{sequenced-$id\} == \"v\" && \$\{alleleSeq1-$id\} == \$\{alleleSeq2-$id\}"
	}
	set needed [list_remdup $needed]
	set poss [list_find [list_cor $header $needed] -1]
	if {[llength $poss]} {
		error "Could not find fields needed for hovar([join $ids ,]): [list_sub $needed $poss]"
	}
	lappend neededfields {*}$needed
	set temp "([join $temp {) && (}])"
	return $temp
}

proc tsv_select_count {arguments header neededfieldsVar} {
	upvar $neededfieldsVar neededfields
	set test [list_pop arguments]
	set temp {}
	foreach q $arguments {
		lappend q {*}$test
		set q [tsv_select_precedence $q]
		lappend temp "([tsv_select_detokenize $q $header neededfields])"
	}
	return "([join $temp " + "])"
}

proc tsv_select_if {arguments header neededfieldsVar} {
	upvar $neededfieldsVar neededfields
	if {[llength $arguments] < 3} {
		error "wrong # args for function if, must be: if(condition1,true1,?condition2?,?true2?,...,false)"
	}
	set q [list_pop arguments]
	set q [tsv_select_precedence $q]
	set false "([tsv_select_detokenize $q $header neededfields])"
	set result {}
	foreach {true condition} [list_reverse $arguments] {
		set true [tsv_select_precedence $true]
		set true "([tsv_select_detokenize $true $header neededfields])"
		set condition [tsv_select_precedence $condition]
		set condition "([tsv_select_detokenize $condition $header neededfields])"
		set result "($condition ? $true : $false)"
		set false $result
	}
	return $result
}

proc tsv_select_catch {arguments header neededfieldsVar} {
	upvar $neededfieldsVar neededfields
	if {[llength $arguments] != 1} {
		error "wrong # args for function if, must be: catch(expression)"
	}
	return "\[catch \{expr \{[tsv_select_detokenize [lindex $arguments 0] $header neededfields]\}\}\]"
}

proc tsv_select_counthasone {ids operator value} {
	set temp {}
	foreach id $ids {
		lappend temp "hasone($id, \"$operator\", $value)"
	}
	return "([join $temp " + "])"
}

proc tsv_select_counthasall {ids operator value} {
	set temp {}
	foreach id $ids {
		lappend temp "hasall($id, \"$operator\", $value)"
	}
	return "([join $temp " + "])"
}

proc tsv_select_region {ids header neededfieldsVar} {
	upvar $neededfieldsVar neededfields
	set ids [split $ids ":-_ ,\{\}\""]
	set ids [list_remove $ids {}]
	set poss [tsv_basicfields $header 3]
	set fields [list_sub $header $poss]
	foreach {rchr rbegin rend} $fields break
	lappend neededfields {*}$fields
	set result {}
	foreach {chr begin end} $ids {
		set q "samechr(\"$chr\",\$$rchr)"
		if {$begin ne ""} {
			lappend q "($begin) < \$$rend"
		}
		if {$end ne ""} {
			lappend q "($end) > \$$rbegin"
		}
		lappend result "([join $q " && "])"
	}
	return "(([join $result "\) || \("]))"
}

proc tsv_select_expandfield {header field} {
	if {$field eq "ROW"} {return "ROW"}
	set qposs [list_find -glob $header $field]
	if {![llength $qposs]} {
#		error "no fields matched \"$field\""
return [list $field]
	}
	set result [list_sub $header $qposs]
	return [list_remdup $result]
}

proc tsv_select_expandcalcfield {header fielddef} {
	set pos [string first = $fielddef]
	set calcfield [string trim [string range $fielddef 0 [expr {$pos-1}]]]
	set code [string range $fielddef [expr {$pos+1}] end]
	set fields [tsv_select_extractvars $code]
	unset -nocomplain ha
	foreach field $fields {
		set wildvars [regexp -inline -all {\*+} $field]
		if {![llength $wildvars]} continue
		set poss [list_find -glob $header $field]
		set hfields [list_sub $header $poss]
		regsub -all {\*+} $field {(.+)} pattern
		set pattern ^$pattern\$
		unset -nocomplain a
		foreach hfield $hfields {
			set vals [lrange [regexp -inline $pattern $hfield] 1 end]
			foreach var $wildvars val $vals {
				lappend a($var) $val
			}
		}
		foreach var [array names a] {
			if {![info exists ha($var)]} {
				set ha($var) [list_remdup $a($var)]
			} else {
				set ha($var) [list_common $ha($var) [list_remdup $a($var)]]
			}
		}
	}
	set wildvars [regexp -inline -all {\*+} $calcfield]
	set vars [array names ha]
	set notpresent [list_lremove $wildvars $vars]
	if {[llength $notpresent]} {
		error "some patterns ([join $notpresent ,]) in calcfield were not matched in code"
	}
	set result [list $fielddef]
	foreach var [lsort -decreasing $wildvars] {
		set newresult {}
		foreach val $ha($var) {
			foreach line $result {
				lappend newresult [string_change $line [list $var $val]]
			}
		}
		set result $newresult
	}
	return $result
}

proc tsv_select_expandfields {header qfields qpossVar} {
	upvar $qpossVar qposs
	upvar tsv_funcnum tsv_funcnum
	set qposs {}
	set rfields {}
	foreach field $qfields {
		if {$field eq "ROW"} {
			lappend rfields ROW
#			lappend qposs {}
#			# empty for code will output ROW directly later
			lappend qposs ROW
			continue
		}
		set pos [string first = $field]
		if {$pos != -1} {
			set fieldname [string trim [string range $field 0 [expr {$pos-1}]]]
			if {[string first * $fieldname] != -1} {
				set list [tsv_select_expandcalcfield $header $field]
			} else {
				set list [list $field]
			}
			foreach field $list {
				set pos [string first = $field]
				set fieldname [string trim [string range $field 0 [expr {$pos-1}]]]
				set fpos [lsearch $rfields $fieldname]
				set code [string range $field [expr {$pos+1}] end]
				if {$fpos == -1} {
					lappend rfields $fieldname
					lappend qposs [list code $code]
				} else {
					lset qposs $fpos [list code $code]
				}
			}
		} elseif {[string first * $field] != -1} {
			set efields [list_lremove [tsv_select_expandfield $header $field] $rfields]
			lappend rfields {*}$efields
			foreach pos [list_cor $header $efields] {
				lappend qposs $pos
			}
		} else {
			if {[inlist $rfields $field]} continue
			set pos [lsearch $header $field]
			if {$pos == -1} {
				error "field \"$field\" not present"
			}
			lappend rfields $field
			lappend qposs $pos
		}
	}
	return $rfields
}

array set tsv_select_tokenize_opsa {
	~ {u 1} ! {u 1} 
	** {d 2}
	@** {d 2}
	* {d 3} / {d 3} % {d 3}
	@* {d 3} @/ {d 3} @% {d 3}
	- {d 4} + {d 4} 
	@- {d 4} @+ {d 4} 
	<< {d 5} >> {d 5}
	< {d 6} > {d 6} <= {d 6} >= {d 6}
	@< {d 6} @> {d 6} @<= {d 6} @>= {d 6}
	== {d 7} = {d 7} != {d 7}
	@== {d 7} @!= {d 7}
	eq {d 8} ne {d 8}
	in {s 9} ni {s 9}
	vin {s 9} vni {s 9}
	& {d 10} ^ {d 11} | {d 12}
	&& {d 13} || {d 14}
	@&& {d 13} @|| {d 14}
	vand {d 13} vor {d 14}
	? {t 15} : {t 16}
}

array set tsv_select_tokenize_newopsa {
	> gtchecked >= gtechecked < ltchecked <= ltechecked
	@** vpower
	@* vtimes @/ vdivide @% vmod
	@- vminus @+ vplus 
	@> vgt @< vlt @>= vgte @<= vlte
	@== veq @!= vne
	vin vin vni vni
	@&& vand @|| vor
	vand vand vor vor
}

proc tsv_select_tokenize_opsdata op {
	global tsv_select_tokenize_opsa
	if {$op eq "and"} {set op &&} elseif {$op eq "or"} {set op ||}
	if {[info exists tsv_select_tokenize_opsa($op)]} {
		foreach {numargs prec} $tsv_select_tokenize_opsa($op) break
		if {[info exists ::tsv_select_tokenize_newopsa($op)]} {
			set type @newop
			set op $::tsv_select_tokenize_newopsa($op)
		} else {
			set type @op
		}
	} else {
		set type @newop
		set prec 0
		set numargs 2
	}
	list $type $op $prec
}

proc tsv_select_extractvars {code} {
	set len [string length $code]
	set prevpos 0
	set pos 0
	set result {}
	while {$pos < $len} {
		set pos [lindex [regexp -start $pos -inline -indices {\$} $code] 0 0]
		if {$pos eq ""} break
		incr pos
		set prevpos $pos
		set char [string index $code $pos]
		if {$char eq "\{"} {
			set pos [lindex [regexp -start $pos -inline -indices \} $code] 0 0]
			lappend result [string range $code [expr {$prevpos+1}] [expr {$pos-1}]]
			set prevpos [expr {$pos+1}]
		} else {
			while 1 {
				set pos [lindex [regexp -start $pos -inline -indices {[^A-Za-z0-9._*]|$} $code] 0 0]
				set char [string index $code $pos]
				if {$char ne "-"} break
				incr pos
				set char [string index $code $pos]
				if {![regexp {[A-Za-z.*]} $char]} break
			}
			lappend result [string range $code $prevpos [expr {$pos-1}]]
			set prevpos $pos
		}
	}
	return $result
}

proc tsv_select_tokenize {header code neededfieldsVar} {
	global tsv_select_tokenize_opsa
	upvar $neededfieldsVar neededfields
	upvar tsv_funcnum tsv_funcnum
	# variable preprocessor first, to expand *
	# check and exchange variables needed
	# variables are all changed to quoted format: ${variable name}
	# subsequent code expects this quoting!
	set code [string trim $code]
	set newcode {}
	set len [string length $code]
	set prevpos 0
	set pos 0
	while {$pos < $len} {
		set pos [lindex [regexp -start $pos -inline -indices {\$} $code] 0 0]
		if {$pos eq ""} {
			append newcode [string range $code $prevpos end]
			break
		}
		append newcode [string range $code $prevpos $pos]
		incr pos
		set prevpos $pos
		set char [string index $code $pos]
		if {$char eq "\{"} {
			set pos [lindex [regexp -start $pos -inline -indices \} $code] 0 0]
			set field [string range $code [expr {$prevpos+1}] [expr {$pos-1}]]
			set fields [tsv_select_expandfield $header $field]
			if {![llength $fields]} {
				error "field \"$field\" not present"
			}
			append newcode \{[join $fields \},\$\{]\}
			lappend neededfields {*}[list_common $header $fields]
			set prevpos [expr {$pos+1}]
		} else {
			while 1 {
				set pos [lindex [regexp -start $pos -inline -indices {[^A-Za-z0-9._*]|$} $code] 0 0]
				set char [string index $code $pos]
				if {$char ne "-"} break
				incr pos
				set char [string index $code $pos]
				if {![regexp {[A-Za-z.*]} $char]} break
			}
			set field [string range $code $prevpos [expr {$pos-1}]]
			set fields [tsv_select_expandfield $header $field]
			if {![llength $fields]} {
				error "field \"$field\" not present"
			}
			append newcode \{[join $fields \},\$\{]\}
			lappend neededfields {*}[list_common $header $fields]
			set prevpos $pos
		}
	}
	set code $newcode
	#
	# tokenize
	# returns list of values and ops
	# elements of this list can contain sublists (braces, functions)
	# pos indicates current position, and goes over the string, detecting values, ops, etc 
	# and lappending the codes to curstack
	# braces, functions will put curstack to stack and start new curstack
	# on a closing brace sublists are managed, and prev curstack is returned
	set len [string length $code]
	set prevpos 0
	set pos 0
	set stack {}
	set curstack {}
	set curtype {}
	set work {}
	while {$pos < $len} {
		set char [string index $code $pos]
		while {[regexp {[ \t\n]} $char]} {
			incr pos
			set char [string index $code $pos]
		}
		set prevpos $pos
		if {[regexp {[0-9.]} $char]} {
			# number
			incr pos
			set pos [lindex [regexp -start $pos -inline -indices {[^0-9.]|$} $code] 0 0]
			lappend curstack [list @num [string range $code $prevpos [expr {$pos-1}]]]
		} elseif {[regexp {[+!*/%<>&^|@?:=-]} $char]} {
			# operator
			incr pos
			set pos [lindex [regexp -start $pos -inline -indices {[^+!*/%<>&^|?@:=-]|$} $code] 0 0]
			set op [string range $code $prevpos [expr {$pos-1}]]
			lappend curstack [tsv_select_tokenize_opsdata $op]
		} elseif {[regexp {[A-Za-z]} $char]} {
			# function or text operator
			set pos [lindex [regexp -start $pos -inline -indices {[^A-Za-z0-9_]|$} $code] 0 0]
			set char [string index $code $pos]
			if {$char eq "\("} {
				# function
				lappend curstack [list @function [string range $code $prevpos [expr {$pos-1}]]]
				lappend stack $curstack
				set curstack {}
				incr pos
			} else {
				# text operator
				set op [string range $code $prevpos [expr {$pos-1}]]
				lappend curstack [tsv_select_tokenize_opsdata $op]
			}
		} elseif {$char eq "\("} {
			lappend curstack [list @braces]
			lappend stack $curstack
			set curstack {}
			incr pos
			set prevpos $pos
		} elseif {$char eq "\$"} {
			# variable
			set prevpos $pos
			incr pos 2
			set pos [string first \} $code $pos]
			lappend curstack [list @var [string range $code [expr {$prevpos+2}] [expr {$pos-1}]]]
			incr pos
		} elseif {$char eq ","} {
			set prevstack [list_pop stack]
			set temp [list_pop prevstack]
			if {[lindex $temp 0] ne "@function"} {
				error "Error: unexpected \",\" outside function argument list on position $pos in $code"
			}
			set curstack [tsv_select_precedence $curstack]
			lappend temp $curstack
			lappend prevstack $temp
			lappend stack $prevstack
			set curstack {}
			incr pos
		} elseif {$char eq "\)"} {
			set prevstack [list_pop stack]
			set prevtype [lindex $prevstack end 0]
			if {$prevtype eq "@function"} {
				set temp [list_pop prevstack]
				set curstack [tsv_select_precedence $curstack]
				lappend temp $curstack
				lappend prevstack $temp
				set curstack $prevstack
				incr pos
			} elseif {$prevtype eq "@braces"} {
				set temp [list_pop prevstack]
				set curstack [tsv_select_precedence $curstack]
				lappend temp {*}$curstack
				lappend prevstack $temp
				set curstack $prevstack
				incr pos
			} else {
				error "Error: unexpected \")\""
			}
		} elseif {$char eq "\{"} {
			incr pos
			set pos [string first \} $code $pos]
			lappend curstack [list @val [string range $code $prevpos $pos]]
			incr pos
		} elseif {$char eq "\""} {
			set prevpos $pos
			incr pos
			while {1} {
				set pos [string first \" $code $pos]
				if {$pos == -1} {
					set quoted [string range $code $prevpos end]
					break
				}
				set quoted [string range $code $prevpos $pos]
				if {[info complete $quoted]} break
				incr pos
			}
			if {![info complete $quoted]} {
				error "error: incomplete quoted expression: $quoted"
			}
			lappend curstack [list @val $quoted]
			incr pos
		} elseif {$char eq "\["} {
			incr pos
			set prevpos $pos
			while {1} {
				set pos [string first \] $code $pos]
				if {$pos == -1} break
				set command [string range $code $prevpos [expr {$pos-1}]]
				if {[info complete $command]} break
				incr pos
			}
			if {![info complete $command]} {
				error "error: incomplete command in expression"
			}
			lappend curstack [list @command $command]
			incr pos
		} elseif {$char eq "~"} {
			# awk style regexp
			incr pos
			set char [string index $code $pos]
			while {[regexp {[ \t\n]} $char]} {
				incr pos
				set char [string index $code $pos]
			}
			if {$char eq "/"} {
				incr pos
				set prevpos $pos
				lappend curstack [list @newop regexp]
				set pos [string first / $code $pos]
				lappend curstack [list @val \{[string range $code $prevpos [expr {$pos-1}]]\}]
				incr pos
			} else {
				lappend curstack [list @op ~]
			}
		} else {
			error "????"
		}
	}
	if {$stack ne ""} {
		error "unbalanced expression"
	}
	set curstack [tsv_select_precedence $curstack]
	return $curstack
}

proc tsv_select_precedence {curstack} {
	if {[llength $curstack] < 3} {return $curstack}
	set poss [list_find -glob $curstack @newop*]
	# only need to do precedence myself using braces if there are @newops
	if {![llength $poss]} {return $curstack}
	lappend poss {*}[list_find -glob $curstack @op*]
	set poss [lsort -integer -decreasing $poss]
	set prev -1
	unset -nocomplain todoops
	foreach pos $poss {
		if {$prev - $pos == 1} {
			set temp [expr {$prev+1}]
			set curstack [lreplace $curstack $prev $temp [list @braces {*}[lrange $curstack $prev $temp]]]
		} elseif {$prev != -1} {
			set temp [lindex $curstack $prev]
			set level [lindex $temp 2]
			if {![isint $level]} {set level 0}
			lappend todoops($level) $temp
		}
		set prev $pos
	}
	if {$prev == 0} {
		set curstack [lreplace $curstack 0 1 [list @braces {*}[lrange $curstack 0 1]]]
	} else {
		set temp [lindex $curstack $prev]
		set level [lindex $temp 2]
		if {![isint $level]} {set level 0}
		lappend todoops($level) $temp
	}
	# join $curstack \n
	set levels [array names todoops]
	if {[llength $levels] > 0} {
		set levels [lsort -integer $levels]
		foreach level $levels {
			set poss {}
			foreach todoop [list_remdup $todoops($level)] {
				lappend poss {*}[list_find $curstack $todoop]
			}
			set poss [lsort -integer $poss]
			set shift 0
			foreach pos $poss {
				set from [expr {$pos+$shift-1}]
				set to [expr {$pos+$shift+1}]
				foreach {pre op post} [lrange $curstack $from $to] break
				if {[lindex $op 0] eq "@op"} {
					set curstack [lreplace $curstack $from $to [list @braces $pre $op $post]]
				} else {
					set op [lindex $op 1]
					set curstack [lreplace $curstack $from $to [list @function $op [list $pre] [list $post]]]
				}
				incr shift -2
			}
		}
	}
	return $curstack
}

proc tsv_select_detokenize {tokens header neededfieldsVar} {
	global newoptransa
	upvar $neededfieldsVar neededfields
	set result {}
	foreach line $tokens {
		foreach {type val} $line break
		switch -exact $type {
			@num - @val {
				lappend result $val
			}
			@var {
				lappend result \$\{$val\}
				lappend neededfields $val
			}
			@op {
				if {$val eq "="} {
					set val ==
				}
				lappend result $val
			}
			@newop {
				lappend result $val
			}
			@function {
				set saggra {
					scount {}
					slist {vector slist_cond_}
					sdistinct {distinct sdistinct_cond_}
					sucount {ucount sucount_cond_}
					smin {lmin smin_cond_}
					smax {lmax smax_cond_}
					ssum {lsum ssum_cond_}
					savg {lavg savg_cond_}
					sstdev {lstdev sstdev_cond_}
					smedian {lmedian smedian_cond_}
					smode {lmode smode_cond_}
					spercent {}
				}
				set arguments [lrange $line 2 end]
				if {[dict exists $saggra $val]} {
					switch $val {
						scount {
							set temp [tsv_select_scount $arguments $header neededfields]
							lappend result $temp
							continue
						}
						spercent {
							set temp [tsv_select_spercent $arguments $header neededfields]
							lappend result $temp
							continue
						}
						default {
							set temp [tsv_select_saggr {*}[dict get $saggra $val] $arguments $header neededfields]
							lappend result $temp
							continue
						}
					}
				}
				set ids {}
				foreach el $arguments {
					lappend ids [tsv_select_detokenize $el $header neededfields]
				}
				switch $val {
					sm {
						set temp [tsv_select_sm $ids $header neededfields]
					}
					same {
						set temp [tsv_select_same $ids $header neededfields]
					}
					df {
						set temp [tsv_select_df $ids $header neededfields]
					}
					mm {
						set temp [tsv_select_mm $ids $header neededfields]
					}
					un {
						set temp [tsv_select_un $ids $header neededfields]
					}
					hovar {
						set temp [tsv_select_hovar $ids $header neededfields]
					}
					compare {
						set temp [tsv_select_compare $ids $header neededfields]
					}
					zyg {
						set temp [tsv_select_zyg $ids $header neededfields]
					}
					count {
						set temp [tsv_select_count $arguments $header neededfields]
					}
					if {
						set temp [tsv_select_if $arguments $header neededfields]
					}
					catch {
						set temp [tsv_select_catch $arguments $header neededfields]
					}
					region {
						set temp [tsv_select_region $ids $header neededfields]
					}
					counthasone {
						foreach {operator value} [lindex $line end] break
						set operator [tsv_select_detokenize [list $operator] $header neededfields]
						set value [tsv_select_detokenize [list $value] $header neededfields]
						set temp [tsv_select_counthasone [lrange $ids 0 end-1] $operator $value]
					}
					counthasall {
						foreach {operator value} [lindex $line end] break
						set operator [tsv_select_detokenize [list $operator] $header neededfields]
						set value [tsv_select_detokenize [list $value] $header neededfields]
						set temp [tsv_select_counthasall [lrange $ids 0 end-1] $operator $value]
					}
					hasone {
						if {[llength $ids] == 2} {
							foreach {operator value} [lindex $line end] break
							set operator [tsv_select_detokenize [list $operator] $header neededfields]
							set value [tsv_select_detokenize [list $value] $header neededfields]
							set temp "${val}\([lindex $ids 0], \"$operator\", $value\)"
						} else {
							set temp "${val}\([join $ids ", "]\)"
						}
					}
					default {
						set temp "${val}\([join $ids ", "]\)"
					}
				}
				lappend result $temp
			}
			@braces {
				lappend result \([tsv_select_detokenize [lrange $line 1 end] $header neededfields]\)
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

proc tsv_select_expandcode {header code neededfieldsVar} {
	upvar $neededfieldsVar neededfields
	# variable preprocessor first, to expand *
	# check and exchange variables needed
	set tokens [tsv_select_tokenize $header $code neededfields]
	# detokenize, making necessary changes
	tsv_select_detokenize $tokens $header neededfields
}

proc tsv_select {query {qfields {}} {sortfields {}} {newheader {}} {sepheader {}} {f stdin} {out stdout} {hc 0} {inverse 0} {group {}} {groupcols {}} {index {}} {verbose -1} {samplingskip 0}} {
#putsvars query qfields sortfields newheader sepheader f out hc inverse group groupcols index verbose samplingskip
	fconfigure $f -buffering none
	fconfigure $out -buffering none
	if {$hc ni {0 1 2}} {
		set hf [gzopen $hc]
		set header [tsv_open $hf keepheader]
		close $hf
	} else {
		set header [tsv_open $f keepheader]
		if {$hc > 0} {
			tsv_hcheader $f keepheader header
			if {$hc == 2} {
				append keepheader \#
			}
		}
	}
	set neededfields {}
	set sort ""
	set cut ""
	set tsv_funcnum 1
	set qfields [tsv_select_expandfields $header $qfields qposs]
	if {$inverse} {
		set qfields [list_lremove $header $qfields]
		set qfields [tsv_select_expandfields $header $qfields qposs]
	}
	if {![file exists $index/cols]} {set index {}}
	if {[llength $qfields] == [llength $header]} {set index {}}
	if {[llength $sortfields]} {
		set index {}
		if {$sortfields eq "-"} {
			set poss [list_sub [tsv_basicfields $header 6 0] -exclude 4]
			set poss [list_remove $poss -1]
		} else {
			set poss [list_cor $header $sortfields]
		}
		if {[lsearch $poss -1] != -1} {error "fields [join [list_sub $sortfields [list_find $poss -1]] ,] not found"}
		if {[llength $poss]} {
			set poss [lmath_calc $poss + 1]
			set keys {}
			foreach pos $poss {
				lappend keys $pos,$pos
			}
			set sort "gnusort8 -T \"[scratchdir]\" -t \\t -N -s -k[join $keys " -k"]"
		}
	}
	set query [string trim $query]
	if {$query ne ""} {
		set pquery [tsv_select_expandcode $header $query neededfields]
	} else {
		set pquery 1
	}
	set pipe {}
	if {$sort ne ""} {
		lappend pipe $sort
	}
# putslog stderr ----------\n$query\n----------
	if {$group ne ""} {
		append tclcode \n [tsv_select_group $header $pquery $qposs $qfields $group $groupcols $neededfields $verbose]
		#file_write /tmp/temp.tcl $tclcode\n
		#putsvars tclcode
		if {[string length $tclcode] > 2000} {
			set tempfile [tempfile]
			file_write $tempfile $tclcode\n
			#file copy $tempfile /tmp/temp.txt
			lappend pipe [list cg source $tempfile]
		} else {
			#putsvars tclcode
			lappend pipe [list cg exec $tclcode]
		}
	} elseif {$query ne "" || [llength $qfields]} {
		set outcols {}
		set num 0
		# start making code
		# neededfields will contain all fields needed by the query proc and by the output procs
		# as the same args will be used for all, start with a placeholder (@neededfields@) here 
		# that will be filled in at the end
		# outcols contains the column number for plain fields and a proc name for calculated fields
		set tclcode "package require genomecomb\n"
		foreach el $qposs field $qfields {
			if {[string index $field 0] eq "-"} {set rfield [string range $field 1 end]} else {set rfield $field}
			if {[isint $el]} {
				lappend outcols $el
			} elseif {$el eq "ROW"} {
				# empty instead of code or int will directly output ROW
				lappend outcols {}
			} elseif {$el eq ""} {
				error "field $rfield not present"
			} elseif {[info exists calccols($rfield)]} {
				if {$rfield eq $field} {lappend outcols $calccols($rfield)}
			} else {
				set code [tsv_select_expandcode $header [lindex $el 1] neededfields]
				if {$rfield eq $field} {lappend outcols make_col$num}
				append tclcode [tsv_select_makecol make_col$num $code]
				set calccols($rfield) make_col$num
			}
			incr num
		}
		# neededcols contains the column numbers used to supply neededfields by tsv_selectc to the run proc
		set neededfields [list_remdup $neededfields]
		# see what we need of calculated fields
		set calcfieldsquery [list_lremove $neededfields $header]
		set calcfieldsquery [list_remove $calcfieldsquery ROW]
		set neededfields [list_lremove $neededfields $calcfieldsquery]
		set prequery {}
		foreach field $calcfieldsquery {
			if {[info exists calccols($field)]} {
				append prequery "\t\t\tset \{$field\} \[$calccols($field) \$\{[join $neededfields \}\ \$\{]\}\]\n"
			} else {
				error "field \"$field\" not present in file"
			}
			
		}
		set neededcols [list_cor $header $neededfields]
		if {$index eq ""} {
			append tclcode [subst {
				proc tsv_selectc_query {$neededfields} {
					$prequery
					expr {$pquery}
				}
				tsv_selectc tsv_selectc_query [list $neededcols] [list $outcols] $verbose $samplingskip
				exit
			}]
		} else {
			set index [file_absolute $index]
			set indexcols [list_remove $neededcols -1]
			foreach el $outcols {
				if {[isint $el]} {lappend indexcols $el}
			}
			set indexcols [list_remdup $indexcols]
			set temp {}
			foreach el $outcols {
				if {[isint $el]} {
					lappend temp [lsearch $indexcols $el]
				} else {
					lappend temp $el
				}
			}
			set outcols $temp
			set neededcols [list_cor $indexcols $neededcols]
			set indexfiles {}
			foreach col $indexcols {
				lappend indexfiles $index/cols/[lindex $header $col].col
			}
			append tclcode [subst {
				proc tsv_selectc_query {$neededfields} {
					$prequery
					expr {$pquery}
				}
				tsv_selectc_indexed tsv_selectc_query [list $neededcols] [list $outcols] [list $indexfiles]
				exit
			}]
		}
		set tclcode [string_change $tclcode [list @neededfield@ $neededfields]]
		if {[string length $tclcode] > 2000} {
			set tempfile [tempfile]
			file_write $tempfile $tclcode\n
			#file copy $tempfile /tmp/temp.txt
			lappend pipe [list cg source $tempfile]
		} else {
			#putsvars tclcode
			lappend pipe [list cg exec $tclcode]
		}
	}
#putslog -------------pipe-------------------
#putslog pipe:[join $pipe " | "]
#putslog ------------------------------------
	if {$qfields ne ""} {
		set nh [list_sub $qfields -exclude [list_find -glob $qfields -*]]
	} else {
		set nh $header
	}
	if {$group ne ""} {
	} elseif {$sepheader ne ""} {
		file_write $sepheader ${keepheader}[join $header \t]\n
	} elseif {[llength $newheader]} {
		if {[llength $newheader] != [llength $nh]} {error "new header (-nh) of wrong length for query results"}
		puts $out ${keepheader}[join $newheader \t]
	} else	{
		puts $out ${keepheader}[join $nh \t]
	}
	if {![llength $pipe]} {
		if {[info exists ::filebuffer($f)]} {
			foreach line $::filebuffer($in) {
				puts $o $line
			}
			unset ::filebuffer($f)
		}
		fcopy $f $out
	} elseif {$index ne ""} {
		exec {*}[join $pipe " | "] >@ $out 2>@ stderr
	} else {
		chanexec $f $out [join $pipe " | "]
	}
}

proc tsv_hcheader {f keepheaderVar headerVar} {
	upvar $keepheaderVar keepheader
	upvar $headerVar header
	if {[catch {
		# this does not work on a stream
		seek $f [expr {-[string length $header] - 1}] current
	}]} {
		# this has to be explicitely supported downstream, only tsv_select does this!
		set ::filebuffer($f) [list [join $header \t]]
	}
	set temp [split [string trimright $keepheader] \n]
	set header [split [string range [list_pop temp] 1 end] \t]
	set keepheader [join $temp \n]\n
}

proc cg_select {args} {
	if {[llength $args] == 0} {
		errorformat select
		exit 1
	}
	set query {}; set fields {}; set sortfields {}; set newheader {}; set sepheader ""; set hc 0; set inverse 0; set group {}; set groupcols {} ; set verbose -1; set samplingskip 0
	set pos 0
	foreach {key value} $args {
		switch -- $key {
			-q {
				set query $value
				if {[regexp {[^=!><\\]=[^=]} $query]} {puts stderr "you may have used = instead of == in query"}
				regsub -all {\\=} $query = query
			}
			-qf {
				set f [gzopen $value]
				set header [tsv_open $f]
				set data [csv_file $f \t]
				close $f
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
			-hc {
				if {$value ni {0 1 2}} {error "-hc must be 0, 1 or 2"}
				set hc $value
			}
			-hf {set hc $value}
			-s {set sortfields $value}
			-n {
				if {$value eq ""} {
					set header [tsv_open stdin]
				} else {
					set f [gzopen $value]
					set header [tsv_open $f]
					catch {close $f}
				}
				puts stdout [join [samples $header] \n]
				exit 0
			}
			-h {
				if {$value eq ""} {
					set header [tsv_open stdin]
				} else {
					set f [gzopen $value]
					set header [tsv_open $f]
					catch {close $f}
				}
				puts stdout [join $header \n]
				exit 0
			}
			-v {
				set verbose $value
			}
			-samplingskip {
				set samplingskip $value
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
	if {[llength $group] && $samplingskip} {
		error "cannot use -samplingskip option with -g option"
	}
	set args [lrange $args $pos end]
	# clean fields and query: remove comments \n anf \t to space
	regsub -all {\n#[^\n]*} $fields {} fields
	regsub -all {\n#[^\n]*} $query {} query
	regsub -all {\n|\t} $query { } query
	set query [string trim $query]
#puts stderr [list fields=$fields query=$query]
	if {[llength $args] > 0} {
		set filename [lindex $args 0]
		set index $filename.index
		set f [gzopen $filename]
	} else {
		set index {}
		set f stdin
	}
	if {[llength $args] > 1} {
		set outfile [lindex $args 1]
		set o [open $outfile w]
	} else {
		set o stdout
	}
	set error [catch {tsv_select $query $fields $sortfields $newheader $sepheader $f $o $hc $inverse $group $groupcols $index $verbose $samplingskip} result]
	if {$f ne "stdin"} {catch {close $f}}
	if {$o ne "stdout"} {catch {close $o}}
	if {$error} {
		puts stderr $result
		exit 1
	}
}

if {[info exists argv0] && [file tail [info script]] eq [file tail $argv0]} {
	package require pkgtools
	set appdir [file dir [pkgtools::startdir]]
	lappend auto_path $appdir/lib
	append env(PATH) :[file dir [file dir $appdir]]/bin:$appdir/bin
	package require Extral
	set ::base [file tail [info script]]
	cg_select {*}$argv
}
