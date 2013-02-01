#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc tsv_select_sm {ids neededfieldsVar} {
	upvar $neededfieldsVar neededfields
	set id1 [string trim [list_pop ids]]
	set temp [list "(\$\{sequenced-$id1\} == \"v\")"]
	lappend neededfields sequenced-$id1 alleleSeq1-$id1 alleleSeq2-$id1
	foreach id $ids {
		set id [string trim $id]
		lappend neededfields sequenced-$id alleleSeq1-$id alleleSeq2-$id
		lappend temp "(\$\{sequenced-$id\} == \"v\")"  "samegeno(\$\{alleleSeq1-$id1\},\$\{alleleSeq2-$id1\},\$\{alleleSeq1-$id\},\$\{alleleSeq2-$id\})"
	}
	set temp "([join $temp " && "])"
}

proc tsv_select_same {ids neededfieldsVar} {
	upvar $neededfieldsVar neededfields
	set id1 [string trim [list_pop ids]]
	set seqlist [list "\$\{sequenced-$id1\} != \"u\""]
	lappend neededfields sequenced-$id1 alleleSeq1-$id1 alleleSeq2-$id1
	set temp {}
	foreach id $ids {
		set id [string trim $id]
		lappend neededfields sequenced-$id alleleSeq1-$id alleleSeq2-$id
		lappend seqlist "\$\{sequenced-$id\} != \"u\""
		lappend temp "samegeno(\$\{alleleSeq1-$id1\},\$\{alleleSeq2-$id1\},\$\{alleleSeq1-$id\},\$\{alleleSeq2-$id\})"
	}
	set temp "([join $seqlist " && "] && [join $temp " && "])"
}

proc tsv_select_df {ids neededfieldsVar} {
	upvar $neededfieldsVar neededfields
	set temp1 {}
	set temp2 {}
	set seqlist {}
	foreach id $ids {
		set id [string trim $id]
		lappend neededfields sequenced-$id
		lappend seqlist "(\$\{sequenced-$id\} != \"u\")"
		lappend temp1 "(\$\{sequenced-$id\} == \"v\")"
		lappend temp2 "(\$\{sequenced-$id\} == \"r\")"
	}
	set temp "(([join $seqlist " && " ]) && ([join $temp1 " || " ]) && ([join $temp2 " || "]))"
}

proc tsv_select_mm {header ids neededfieldsVar} {
	upvar $neededfieldsVar neededfields
	set temp1 {}
	set temp2 {}
	set list {}
	set seqlist {}
	foreach id $ids {
		set id [string trim $id]
		lappend seqlist "(\$\{sequenced-$id\} == \"v\")"
		lappend list [list \{alleleSeq1-$id\} \{alleleSeq2-$id\}]
		lappend neededfields sequenced-$id alleleSeq1-$id alleleSeq2-$id
	}
	if {[lsearch $header reference] != -1} {
		set ref reference
	} else {
		set ref ref
	}
	lappend neededfields $ref
	while {[llength $list]} {
		foreach {a1 a2} [list_pop list] break
		lappend temp1 "((\$$a1 != \$$ref) || (\$$a2 != \$$ref))"
		list_foreach {a12 a22} $list {
			lappend temp2 "!samegeno(\$$a1,\$$a2,\$$a12,\$$a22)"
		}
	}
	set temp "(([join $seqlist " && " ]) && ([join $temp1 " && " ]) && ([join $temp2 " || "]))"
}

proc tsv_select_un {ids neededfieldsVar} {
	upvar $neededfieldsVar neededfields
	set temp1 {}
	set temp2 {}
	foreach id $ids {
		set id [string trim $id]
		foreach {a1 a2 sequenced} [tsv_select_idtopos $header $id [list alleleSeq1-$id alleleSeq2-$id sequenced-$id]] break
		lappend temp1 "(\$\{sequenced-$id\} == \"v\")"
		lappend temp2 "(\$\{sequenced-$id\} == \"u\")"
		lappend neededfields sequenced-$id
	}
	set temp "(([join $temp2 " || "]) && ([join $temp1 " || " ]))"
}

proc tsv_select_hovar {ids neededfieldsVar} {
	upvar $neededfieldsVar neededfields
	set temp {}
	foreach id $ids {
		lappend neededfields sequenced-$id alleleSeq1-$id alleleSeq2-$id
		lappend temp "\$\{sequenced-$id\} == \"v\" && \$\{alleleSeq1-$id\} == \$\{alleleSeq2-$id\}"
	}
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

proc tsv_select_counthasone {ids operand value} {
	set temp {}
	foreach id $ids {
		lappend temp "hasone($id, \"$operand\", $value)"
	}
	return "([join $temp " + "])"
}

proc tsv_select_counthasall {ids operand value} {
	set temp {}
	foreach id $ids {
		lappend temp "hasall($id, \"$operand\", $value)"
	}
	return "([join $temp " + "])"
}

proc tsv_select_oneof {ids} {
	set value [list_shift ids]
	set temp {}
	foreach id $ids {
		lappend temp "$value == $id"
	}
	return "([join $temp " || "])"
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
		error "no fields matched \"$field\""
	}
	set result [list_sub $header $qposs]
	return [list_remdup $result]
}

proc tsv_select_expandfields {header qfields qpossVar} {
	upvar $qpossVar qposs
	upvar tsv_funcnum tsv_funcnum
	set qposs {}
	set rfields {}
	foreach field $qfields {
		if {$field eq "ROW"} {
			lappend rfields ROW
			lappend qposs {}
			# empty for code will output ROW directly later
			continue
		}
		set pos [string first = $field]
		if {$pos != -1} {
			set fieldname [string range $field 0 [expr {$pos-1}]]
			lappend rfields $fieldname
			set code [string range $field [expr {$pos+1}] end]
			lappend qposs [list code $code]
		} elseif {[string first * $field] != -1} {
			set efields [tsv_select_expandfield $header $field]
			lappend rfields {*}$efields
			foreach pos [list_cor $header $efields] {
				lappend qposs $pos
			}
		} else {
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
	== {d 7} != {d 7}
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
	set escape 0
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
			lappend neededfields {*}$fields
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
			lappend neededfields {*}$fields
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
			# operand
			incr pos
			set pos [lindex [regexp -start $pos -inline -indices {[^+!*/%<>&^|?@:=-]|$} $code] 0 0]
			set op [string range $code $prevpos [expr {$pos-1}]]
			lappend curstack [tsv_select_tokenize_opsdata $op]
		} elseif {[regexp {[A-Za-z]} $char]} {
			# function or text operand
			set pos [lindex [regexp -start $pos -inline -indices {[^A-Za-z0-9_]|$} $code] 0 0]
			set char [string index $code $pos]
			if {$char eq "\("} {
				# function
				lappend curstack [list @function [string range $code $prevpos [expr {$pos-1}]]]
				lappend stack $curstack
				set curstack {}
				incr pos
			} else {
				# text operand
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
			incr pos
			set pos [string first \" $code $pos]
			lappend curstack [list @val [string range $code $prevpos $pos]]
			incr pos
		} elseif {$char eq "\["} {
			incr pos
			set prevpos $pos
			while {$pos < $len} {
				set pos [string first \] $code $pos]
				set command [string range $code $prevpos [expr {$pos-1}]]
				if {[info complete $command]} break
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

#	upvar $neededfieldsVar neededfields
#	if {[llength $tokens] >= 3} {
#		set ops [list_subindex $tokens 0]
#		set poss [lsort -integer -decreasing [list_find $ops @newop]]
#		foreach pos $poss {
#			foreach {pre op post} [lrange $tokens [expr {$pos-1}] [expr {$pos+1}]] break
#			set tokens [lreplace $tokens [expr {$pos-1}] [expr {$pos+1}] [list @function [lindex $op end] [list $pre] [list $post]]]
#		}
#	}

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
			}
			@op {
				lappend result $val
			}
			@newop {
				lappend result $val
			}
			@function {
				set arguments [lrange $line 2 end]
				set ids {}
				foreach el $arguments {
					lappend ids [tsv_select_detokenize $el $header neededfields]
				}
				switch $val {
					sm {
						set temp [tsv_select_sm $ids neededfields]
					}
					same {
						set temp [tsv_select_same $ids neededfields]
					}
					df {
						set temp [tsv_select_df $ids neededfields]
					}
					mm {
						set temp [tsv_select_mm $header $ids neededfields]
					}
					un {
						set temp [tsv_select_un $ids neededfields]
					}
					hovar {
						set temp [tsv_select_hovar $ids neededfields]
					}
					count {
						set temp [tsv_select_count $arguments $header neededfields]
					}
					region {
						set temp [tsv_select_region $ids $header neededfields]
					}
					counthasone {
						foreach {operand value} [lindex $line end] break
						set operand [tsv_select_detokenize [list $operand] $header neededfields]
						set value [tsv_select_detokenize [list $value] $header neededfields]
						set temp [tsv_select_counthasone [lrange $ids 0 end-1] $operand $value]
					}
					counthasall {
						foreach {operand value} [lindex $line end] break
						set operand [tsv_select_detokenize [list $operand] $header neededfields]
						set value [tsv_select_detokenize [list $value] $header neededfields]
						set temp [tsv_select_counthasall [lrange $ids 0 end-1] $operand $value]
					}
					hasone {
						if {[llength $ids] == 2} {
							foreach {operand value} [lindex $line end] break
							set operand [tsv_select_detokenize [list $operand] $header neededfields]
							set value [tsv_select_detokenize [list $value] $header neededfields]
							set temp "${val}\([lindex $ids 0], \"$operand\", $value\)"
						} else {
							set temp "${val}\([join $ids ", "]\)"
						}
					}
					oneof {
						set temp [tsv_select_oneof $ids]
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
	upvar tsv_funcnum tsv_funcnum
	# variable preprocessor first, to expand *
	# check and exchange variables needed
	set tokens [tsv_select_tokenize $header $code neededfields]
	# detokenize, making necessary changes
	tsv_select_detokenize $tokens $header neededfields
}

proc tsv_select {query {qfields {}} {sortfields {}} {newheader {}} {sepheader {}} {f stdin} {out stdout} {hc 0} {inverse 0}} {
# putsvars query qfields sortfields newheader sepheader f stdin out stdout h inverse
	fconfigure $f -buffering none
	fconfigure $out -buffering none
	if {$hc ne "0" && $hc ne "1"} {
		set hf [gzopen $hc]
		set header [tsv_open $hf keepheader]
		close $hf
	} else {
		set header [tsv_open $f keepheader]
		if {$hc eq "1"} {
			tsv_hcheader $f keepheader header
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

	if {[llength $sortfields]} {
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
	if {$query ne "" || [llength $qfields]} {
		set outcols {}
		set todo {}
		set num 0
		# neededfields will be the same for all procs, so first gather all we need using todo list
		# then make the procs from todo list later
		foreach el $qposs {
			if {[isint $el]} {
				lappend outcols $el
			} elseif {$el eq ""} {
				# empty instead of code or int will directly output ROW
				lappend outcols {}
			} else {
				set code [tsv_select_expandcode $header [lindex $el 1] neededfields]
				lappend outcols make_col$num
				lappend todo make_col$num $code
			}
			incr num
		}
		set neededfields [list_remdup $neededfields]
		set neededcols [list_cor $header $neededfields]
		set tclcode "package require genomecomb\n"
		foreach {name code} $todo {
			append tclcode [subst -nocommands {
				proc $name {$neededfields} {
					if {[catch {expr {$code}} e]} {
						switch \$e {
							{domain error: argument not in valid range} {return NaN}
							{divide by zero} {return NaN}
						}
					}
					return \$e
				}
			}]
		}
		append tclcode [subst {
			proc tsv_selectc_query {$neededfields} {
				expr {$pquery}
			}
			tsv_selectc tsv_selectc_query [list $neededcols] [list $outcols]
			exit
		}]
		lappend pipe [list cg exec $tclcode]
	}
#putslog -------------pipe-------------------
#putslog pipe:[join $pipe " | "]
#putslog ------------------------------------
	if {$qfields ne ""} {
		set nh $qfields
	} else {
		set nh $header
	}
	if {$sepheader ne ""} {
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
	} else {
		chanexec $f $out [join $pipe " | "]
	}
}

proc tsv_hcheader {f keepheaderVar headerVar} {
	upvar $keepheaderVar keepheader
	upvar $headerVar header
	set ::filebuffer($f) [list [join $header \t]]
	set temp [split [string trimright $keepheader] \n]
	set header [split [string range [list_pop temp] 1 end] \t]
	set keepheader [join $temp \n]\n
}

proc cg_select {args} {
	if {[llength $args] == 0} {
		errorformat select
		exit 1
	}
	set query {}; set fields {}; set sortfields {}; set newheader {}; set sepheader ""; set hc 0; set inverse 0
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
			-nh {set newheader $value}
			-sh {set sepheader $value}
			-hc {set hc 1}
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
				set names {}
				foreach col $header {
					set split [split $col -]
					if {[llength $split] > 1} {
						lappend names [lindex $split end]
					}
				}
				puts stdout [join [list_remdup $names] \n]
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
			-- break
			default {
				break
			}
		}
		incr pos 2
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
		set f [gzopen $filename]
	} else {
		set f stdin
	}
	if {[llength $args] > 1} {
		set outfile [lindex $args 1]
		set o [open $outfile w]
	} else {
		set o stdout
	}
	set error [catch {tsv_select $query $fields $sortfields $newheader $sepheader $f $o $hc $inverse} result]
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


