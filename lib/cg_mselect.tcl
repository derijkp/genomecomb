monetdbinit

proc mselect_idtopos {header id fields} {
	set id [string trim $id \"]
	set poss [list_cor $header $fields]
	if {[inlist $poss -1]} {error "sample $id not present"}
	return [lmath_calc $poss + 1]
}

proc mselect_samegeno {a11 a12 a21 a22} {
	return "(($a11 = $a21 and $a12 = $a22) or ($a11 = $a22 and $a12 = $a21))"
}

proc mselect_sm {header ids} {
	set id1 [list_pop ids]
	set id1 [string trim $id1 \"]
	set a11 \"alleleSeq1_$id1\"
	set a21 \"alleleSeq2_$id1\"
	set temp [list "(\"sequenced_$id1\" = \'v\')"]
	foreach id $ids {
		set id [string trim $id \"]
		set a12 \"alleleSeq1_$id\"
		set a22 \"alleleSeq2_$id\"
		lappend temp "(\"sequenced_$id\" = \'v\')"  "[mselect_samegeno $a11 $a21 $a12 $a22]"
	}
	set temp "([join $temp " and "])"
}

proc mselect_same {header ids} {
	set id1 [list_pop ids]
	set id1 [string trim $id1 \"]
	set a11 \"alleleSeq1_$id1\"
	set a21 \"alleleSeq2_$id1\"
	set sequenced \"sequenced_$id1\"
	set seqlist [list "$sequenced <> \'u\'"]
	set temp {}
	foreach id $ids {
		lappend seqlist "$sequenced <> \'u\'"
		set a12 \"alleleSeq1_$id\"
		set a22 \"alleleSeq2_$id\"
		set sequenced \"sequenced_$id\"
		lappend temp "[mselect_samegeno $a11 $a21 $a12 $a22]"
	}
	set temp "([join $seqlist " && "] && [join $temp " && "])"
}

proc mselect_df {header ids} {
	set temp1 {}
	set temp2 {}
	set seqlist {}
	foreach id $ids {
		set a1 \"alleleSeq1_$id\"
		set a2 \"alleleSeq2_$id\"
		set sequenced \"sequenced_$id\"
		lappend seqlist "($sequenced <> \'u\')"
		lappend temp1 "($sequenced = \'v\')"
		lappend temp2 "($sequenced = \'r\')"
	}
	set temp "(([join $seqlist " and " ]) and ([join $temp1 " or " ]) and ([join $temp2 " or "]))"
}

proc mselect_mm {header ids} {
	set temp1 {}
	set temp2 {}
	set list {}
	set seqlist {}
	foreach id $ids {
		set a1 \"alleleSeq1_$id\"
		set a2 \"alleleSeq2_$id\"
		set sequenced \"sequenced_$id\"
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

proc mselect_un {header ids} {
	set temp1 {}
	set temp2 {}
	foreach id $ids {
		foreach {a1 a2 sequenced} [mselect_idtopos $header $id [list alleleSeq1-$id alleleSeq2-$id sequenced-$id]] break
		lappend temp1 "(\$$sequenced == \"v\")"
		lappend temp2 "(\$$sequenced == \"u\")"
	}
	set temp "(([join $temp2 " || "]) && ([join $temp1 " || " ]))"
}

proc mselect_count {ids} {
	set test [list_pop ids]
	set temp {}
	foreach id $ids {
		lappend temp "($id $test)"
	}
	return "([join $temp " + "])"
}
# cg select -q 'count($alleleSeq1,$alleleSeq2, == "G") == 1' annotvar.tsv

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

proc mselect_expandfield {header field qpossVar} {
	upvar $qpossVar qposs
	set qposs [list_find -glob $header $field]
	if {![llength $qposs]} {
		error "no fields matched \"$field\""
	}
	set result [list_sub $header $qposs]
	set qposs [lmath_calc $qposs + 1]
	return $result
}

proc mselect_expandfields {header qfields qcodeVar} {
	upvar $qcodeVar qcode
	upvar tsv_funcnum tsv_funcnum
	set qcode {}
	set rfields {}
	foreach field $qfields {
		set pos [string first = $field]
		if {$pos != -1} {
			lappend rfields [string range $field 0 [expr {$pos-1}]]
			set code [string range $field [expr {$pos+1}] end]
			lappend qcode [mselect_expandcode $header $code]
		} elseif {[string first * $field] != -1} {
			lappend rfields {*}[mselect_expandfield $header $field poss]
			foreach pos $poss {
				lappend qcode {}
			}
		} else {
			set pos [lsearch $header $field]
			if {$pos == -1} {
				error "field \"$field\" not present"
			}
			lappend rfields $field
			lappend qcode {}
		}
	}
	return $rfields
}

proc mselect_expandcode {header code} {
	upvar tsv_funcnum tsv_funcnum
	set indices [list_unmerge [regexp -all -indices -inline {["]?([*a-zA-z0-9_.-]*[*][*a-zA-z0-9_.-]*)["]?} $code]]
	set indices [list_reverse $indices]
	list_foreach {start end} $indices {
		set field [string trim [string range $code $start $end] \"]
		set temp [mselect_expandfield $header $field tposs]
		if {![llength $temp]} {error "field \"$field\" not present"}
		set new {}
		foreach pos $tposs {
			lappend new {}
		}
		set code [string_replace $code $start $end \"[join $temp \",\"]\"]
	}
	set indices [list_unmerge [regexp -all -indices -inline {([a-zA-z0-9_]+)\([^)]+\)} $code]]
	set indices [list_reverse $indices]
	list_foreach {start end} $indices {
		set full [string range $code [expr {$start}] $end]
		if {[regexp {^(.*)\((.*)\)$} $full temp func args]} {
			switch $func {
				sm {
					set ids [split $args ,]
					set temp [mselect_sm $header $ids]
				}
				same {
					set ids [split $args ,]
					set temp [mselect_same $header $ids]
				}
				df {
					set ids [split $args ,]
					set temp [mselect_df $header $ids]
				}
				mm {
					set ids [split $args ,]
					set temp [mselect_mm $header $ids]
				}
				un {
					set ids [split $args ,]
					set temp [mselect_un $header $ids]
				}
				count {
					set ids [split $args ,]
					set temp [mselect_count $ids]
				}
				counthasone {
					set ids [split $args ,]
					set temp [mselect_counthasone $ids]
				}
				counthasall {
					set ids [split $args ,]
					set temp [mselect_counthasall $ids]
				}
				lmin {
					mselect_lmin
				}
				lmax {
					mselect_lmax
				}
				cmin {
					set num [regexp -all , $args]
				}
				cmax {
					set num [regexp -all , $args]
				}
			}
			set code [string_replace $code $start $end $temp]
		} else {
			set pos [lsearch $header $field]
			if {$pos == -1} {error "field \"$field\" not present"}
			incr pos
			set code [string_replace $code $start $end \$$pos]
		}
	}
	return $code
}

proc monetdb_makesql {table header query qfieldsVar {sortfields {}} {inverse 0} {offset {}} {limit {}}} {
	upvar $qfieldsVar qfields
	set sqlfields {}
	set sqlwhere {}
	set sqlsort {}
	if {$qfields eq ""} {
		set sqlfields \"[join $header "\",\""]\"
	} else {
		set qfields [mselect_expandfields $header $qfields qcode]
		if {$inverse} {
			set qfields [list_lremove $header $qfields]
		}
		foreach field $qfields code $qcode {
			if {$code ne ""} {
				lappend sqlfields "$code as \"$field\""
			} else {
				lappend sqlfields \"$field\"
			}
		}
		set sqlfields [join $sqlfields ,]
	}
	if {[llength $sortfields]} {
		set sqlsort "order by \"[join $sortfields "\",\""]\""
	}
	if {$query ne ""} {
		set query [mselect_expandcode $header $query]
		set sqlwhere "where $query"
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

proc monetdb_select {db table header query {qfields {}} {sortfields {}} {newheader {}} {sepheader {}} {out stdout} {hc 0} {inverse 0}} {
#putsvars db table header query qfields sortfields newheader sepheader out hc inverse
	if {$out ne ""} {
		fconfigure $out -buffering none
	}
	set keepheader {}
#	if {$hc ne "0" && $hc ne "1"} {
#		set hf [gzopen $hc]
#		set header [tsv_open $hf keepheader]
#		close $hf
#	} else {
#		set header [tsv_open $f keepheader]
#		if {$hc eq "1"} {
#			tsv_hcheader $f keepheader header
#		}
#	}
	set sql [monetdb_makesql $table $header $query qfields $sortfields $inverse]
#putslog -------------sql--------------------
#putslog sql:\n$sql
#putslog ------------------------------------
	if {$qfields ne ""} {
		set nh $qfields
	} else {
		set nh $header
	}
	if {[llength $newheader]} {
		if {[llength $newheader] != [llength $nh]} {error "new header (-nh) of wrong length for query results"}
		set nh $newheader
	}
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

proc cg_mselect {args} {
	if {[llength $args] == 0} {
		errorformat mselect
		exit 1
	}
	set query {}; set fields {}; set sortfields {}; set newheader {}; set sepheader ""; set hc 0; set inverse 0; set printheader 0
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
				set header [cg_monetdb fields $db $table]
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
				set printheader 1
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	set len [llength $args]
	if {$len < 2 || $len > 3} {
		errorformat select
		exit 1
	}
	foreach {db table outfile} $args break
	set header [cg_monetdb fields $db $table]
	if {$printheader} {
		puts stdout [join $header \n]
		exit 0
	}
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
	set error [catch {monetdb_select $db $table $header $query $fields $sortfields $newheader $sepheader $o $hc $inverse} result]
	if {$o ne "stdout"} {catch {close $o}}
	if {$error} {
		puts stderr $result
return
		exit 1
	}
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
