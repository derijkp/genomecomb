#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc tsv_select_idtopos {header id fields} {
	set poss [list_cor $header $fields]
	if {[inlist $poss -1]} {error "sample $id not present"}
	return [lmath_calc $poss + 1]
}

proc tsv_select_f_samegeno {} {
	upvar awkfunctions awkfunctions
	lappend awkfunctions {
		function samegeno(a11,a12,a21,a22) {
			if ((a11 == a21) && (a12 == a22)) {return 1;}
			if ((a11 == a22) && (a12 == a21)) {return 1;}
			return 0;
		}
	}
}

proc tsv_select_sm {header ids} {
	upvar awkfunctions awkfunctions
	tsv_select_f_samegeno
	set id1 [list_pop ids]
	foreach {a11 a21 sequenced} [tsv_select_idtopos $header $id1 [list alleleSeq1-$id1 alleleSeq2-$id1 sequenced-$id1]] break
	set temp [list "(\$$sequenced == \"v\")"]
	foreach id $ids {
		foreach {a12 a22 sequenced} [tsv_select_idtopos $header $id1 [list alleleSeq1-$id alleleSeq2-$id sequenced-$id]] break
		lappend temp "(\$$sequenced == \"v\")"  "samegeno(\$$a11,\$$a21,\$$a12,\$$a22)"
	}
	set temp "([join $temp " && "])"
}

proc tsv_select_same {header ids} {
	upvar awkfunctions awkfunctions
	tsv_select_f_samegeno
	set id1 [list_pop ids]
	foreach {a11 a21 sequenced} [tsv_select_idtopos $header $id1 [list alleleSeq1-$id1 alleleSeq2-$id1 sequenced-$id1]] break
	set seqlist [list "\$$sequenced != \"u\""]
	set temp {}
	foreach id $ids {
		lappend seqlist "\$$sequenced != \"u\""
		foreach {a12 a22 sequenced} [tsv_select_idtopos $header $id1 [list alleleSeq1-$id alleleSeq2-$id sequenced-$id]] break
		lappend temp "samegeno(\$$a11,\$$a21,\$$a12,\$$a22)"
	}
	set temp "([join $seqlist " && "] && [join $temp " && "])"
}

proc tsv_select_df {header ids} {
	set temp1 {}
	set temp2 {}
	set seqlist {}
	foreach id $ids {
		foreach {a1 a2 sequenced} [tsv_select_idtopos $header $id [list alleleSeq1-$id alleleSeq2-$id sequenced-$id]] break
		lappend seqlist "(\$$sequenced != \"u\")"
		lappend temp1 "(\$$sequenced == \"v\")"
		lappend temp2 "(\$$sequenced == \"r\")"
	}
	set temp "(([join $seqlist " && " ]) && ([join $temp1 " || " ]) && ([join $temp2 " || "]))"
}

proc tsv_select_mm {header ids} {
	upvar awkfunctions awkfunctions
	tsv_select_f_samegeno
	set temp1 {}
	set temp2 {}
	set list {}
	set seqlist {}
	foreach id $ids {
		foreach {a1 a2 sequenced} [tsv_select_idtopos $header $id [list alleleSeq1-$id alleleSeq2-$id sequenced-$id]] break
		lappend seqlist "(\$$sequenced == \"v\")"
		lappend list [list $a1 $a2]
	}
	set ref [lsearch $header reference]
	if {$ref == -1} {
		set ref [lsearch $header ref]
	}
	incr ref
	while {[llength $list]} {
		foreach {a1 a2} [list_pop list] break
		lappend temp1 "((\$$a1 != \$$ref) || (\$$a2 != \$$ref))"
		list_foreach {a12 a22} $list {
			lappend temp2 "!samegeno(\$$a1,\$$a2,\$$a12,\$$a22)"
		}
	}
	set temp "(([join $seqlist " && " ]) && ([join $temp1 " && " ]) && ([join $temp2 " || "]))"
}

proc tsv_select_un {header ids} {
	set temp1 {}
	set temp2 {}
	foreach id $ids {
		foreach {a1 a2 sequenced} [tsv_select_idtopos $header $id [list alleleSeq1-$id alleleSeq2-$id sequenced-$id]] break
		lappend temp1 "(\$$sequenced == \"v\")"
		lappend temp2 "(\$$sequenced == \"u\")"
	}
	set temp "(([join $temp2 " || "]) && ([join $temp1 " || " ]))"
}

proc tsv_select_count {ids} {
	set test [list_pop ids]
	set temp {}
	foreach id $ids {
		lappend temp "($id $test)"
	}
	return "([join $temp " + "])"
}
# cg select -q 'count($alleleSeq1,$alleleSeq2, == "G") == 1' annotvar.tsv

proc tsv_select_lmin {} {
	upvar awkfunctions awkfunctions
	lappend awkfunctions {
		function lmin(list,def) {
			if (def == nill) {def = 999999999}
		        split(list,a,",");
		        minv = a[1];
		        for (i in a) {
		                if (a[i] != a[i]+0) {a[i] = def}
		                if (a[i] < minv) {minv = a[i]}
		        }
		        return minv
		}
	}
}

proc tsv_select_lmax {} {
	upvar awkfunctions awkfunctions
	lappend awkfunctions {
		function lmax(list,def) {
			if (def == nill) {def = -999999999}
		        split(list,a,",");
		        maxv = a[1];
		        for (i in a) {
		                if (a[i] != a[i]+0) {a[i] = def}
		                if (a[i] > maxv) {maxv = a[i]}
		        }
		        return maxv
		}
	}
}

proc tsv_select_min {num} {
	incr num
	set temp [list_fill $num 1 1]
	set result [subst {
		function min(a[join $temp ,a]) \{
			if (a1 ~ /def=(.*)/) {
				def = substr(a1,5)
				if (a2 != a2 + 0) {a2 = def}
				minv = a2
			} else {
				if (a1 != a1 + 0) {a1 = def}
				def = 999999999
				minv = a1
			}
	}]
	foreach field [lrange $temp 1 end] {
		append result [subst {
			if (a$field == nill) {return minv}
			if (a$field != a$field + 0) {a$field = def}
			if (a$field < minv) {minv = a$field}
		}]
	}
	append result [subst {
			return minv
		\}
	}]
	return $result
}

proc tsv_select_max {num} {
	incr num
	set temp [list_fill $num 1 1]
	set result [subst {
		function max(a[join $temp ,a]) \{
			if (a1 ~ /def=(.*)/) {
				def = substr(a1,5)
				if (a2 != a2 + 0) {a2 = def}
				maxv = a2
			} else {
				if (a1 != a1 + 0) {a1 = def}
				def = -999999999
				maxv = a1
			}
	}]
	foreach field [lrange $temp 1 end] {
		append result [subst {
			if (a$field == nill) {return maxv}
			if (a$field != a$field + 0) {a$field = def}
			if (a$field > maxv) {maxv = a$field}
		}]
	}
	append result [subst {
			return maxv
		\}
	}]
	return $result
}

proc tsv_select_counthasone {ids} {
	upvar awkfunctions awkfunctions
	upvar tsv_funcnum tsv_funcnum
	set test [list_pop ids]
	lappend awkfunctions [subst -nocommands {
		function tsvfunc${tsv_funcnum}(list) {
		        split(list,a,",");
		        for (i in a) {
		                if (a[i] $test) {return 1}
		        }
		        return 0
		}
	}]
	set temp {}
	foreach id $ids {
		lappend temp "tsvfunc${tsv_funcnum}($id)"
	}
	return "([join $temp " + "])"
}

proc tsv_select_counthasall {ids} {
	upvar awkfunctions awkfunctions
	upvar tsv_funcnum tsv_funcnum
	set test [list_pop ids]
	lappend awkfunctions [subst -nocommands {
		function tsvfunc${tsv_funcnum}(list) {
		        split(list,a,",");
		        for (i in a) {
		                if (!(a[i] $test)) {return 0}
		        }
		        return 1
		}
	}]
	set temp {}
	foreach id $ids {
		lappend temp "tsvfunc${tsv_funcnum}($id)"
	}
	return "([join $temp " + "])"
}

proc tsv_select_oneof {header ids} {
	set value [list_shift ids]
	set temp {}
	foreach id $ids {
		lappend temp "$value == $id"
	}
	return "([join $temp " || "])"
}
# cg select -q 'oneof($alleleSeq1-dlb_a_d390,"G","C")' annottest_compar.tsv

proc tsv_select_expandfield {header field qpossVar} {
	upvar $qpossVar qposs
	set qposs [list_find -glob $header $field]
	if {![llength $qposs]} {
		error "no fields matched \"$field\""
	}
	set result [list_sub $header $qposs]
	set qposs [lmath_calc $qposs + 1]
	return $result
}

proc tsv_select_expandfields {header qfields qpossVar awkfunctionsVar} {
	upvar $qpossVar qposs
	upvar $awkfunctionsVar awkfunctions
	set qposs {}
	set rfields {}
	foreach field $qfields {
		set pos [string first = $field]
		if {$pos != -1} {
			lappend rfields [string range $field 0 [expr {$pos-1}]]
			set code [string range $field [expr {$pos+1}] end]
			lappend qposs [tsv_select_expandcode $header $code awkfunctions]
		} elseif {[string first * $field] != -1} {
			lappend rfields {*}[tsv_select_expandfield $header $field poss]
			foreach pos $poss {
				lappend qposs \$$pos
			}
		} else {
			set pos [lsearch $header $field]
			if {$pos == -1} {
				error "field \"$field\" not present"
			}
			lappend rfields $field
			lappend qposs \$[expr {$pos+1}]
		}
	}
	return $rfields
}

proc tsv_select_expandcode {header code awkfunctionsVar} {
	upvar $awkfunctionsVar awkfunctions
	set indices [list_unmerge [regexp -all -indices -inline {[$]([*a-zA-z0-9_.-]+)} $code]]
	set indices [list_reverse $indices]
	list_foreach {start end} $indices {
		set field [string range $code [expr {$start+1}] $end]
		if {[string first * $field] == -1} {
			set pos [lsearch $header $field]
			if {$pos == -1} {error "field \"$field\" not present"}
			incr pos
			set code [string_replace $code $start $end \$$pos]
		} else {
			set temp [tsv_select_expandfield $header $field tposs]
			if {![llength $temp]} {error "field \"$field\" not present"}
			set new {}
			foreach pos $tposs {
				lappend new \$$pos
			}
			set code [string_replace $code $start $end [join $new ,]]
		}
	}
	set indices [list_unmerge [regexp -all -indices -inline {([a-zA-z0-9_]+)\([^)]+\)} $code]]
	set indices [list_reverse $indices]
	list_foreach {start end} $indices {
		set full [string range $code [expr {$start}] $end]
		if {[regexp {^(.*)\((.*)\)$} $full temp func args]} {
			switch $func {
				sm {
					set ids [split $args ,]
					set temp [tsv_select_sm $header $ids]
				}
				same {
					set ids [split $args ,]
					set temp [tsv_select_same $header $ids]
				}
				df {
					set ids [split $args ,]
					set temp [tsv_select_df $header $ids]
				}
				mm {
					set ids [split $args ,]
					set temp [tsv_select_mm $header $ids]
				}
				un {
					set ids [split $args ,]
					set temp [tsv_select_un $header $ids]
				}
				count {
					set ids [split $args ,]
					set temp [tsv_select_count $ids]
				}
				counthasone {
					set ids [split $args ,]
					set temp [tsv_select_counthasone $ids]
				}
				counthasall {
					set ids [split $args ,]
					set temp [tsv_select_counthasall $ids]
				}
				oneof {
					set ids [split $args ,]
					set temp [tsv_select_oneof $header $ids]
				}
				lmin {
					tsv_select_lmin
				}
				lmax {
					tsv_select_lmax
				}
				min {
					set num [regexp -all , $args]
					if {$num > [lindex $awkfunctions 0]} {lset awkfunctions 0 $num}
				}
				max {
					set num [regexp -all , $args]
					if {$num > [lindex $awkfunctions 1]} {lset awkfunctions 1 $num}
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

proc tsv_select {query {qfields {}} {sortfields {}} {newheader {}} {f stdin} {out stdout} {hc 0} {inverse 0}} {
	fconfigure $f -buffering none
	fconfigure $out -buffering none
	set header [tsv_open $f keepheader]
	if {$hc} {
		tsv_hcheader $f keepheader header
	}
	set awk ""
	set awkfunctions {0 0}
	set sort ""
	set cut ""
	set qfields [tsv_select_expandfields $header $qfields qposs awkfunctions]
	if {$inverse} {
		set qfields [list_lremove $header $qfields]
		set qfields [tsv_select_expandfields $header $qfields qposs awkfunctions]
	}

	set tsv_funcnum 1
	if {[llength $sortfields]} {
		set poss [list_cor $header $sortfields]
		if {[lsearch $poss -1] != -1} {error "fields [join [list_sub $sortfields [list_find $poss -1]] ,] not found"}
		set poss [lmath_calc $poss + 1]
		set keys {}
		foreach pos $poss {
			lappend keys $pos,$pos
		}
		set sort "gnusort8 -T \"[scratchdir]\" -t \\t -N -s -k[join $keys " -k"]"
	}
	if {($query ne "") || ($qfields ne "")} {
		append awk {BEGIN {FS="\t" ; OFS="\t"}}
		if {$query ne ""} {
			set query [tsv_select_expandcode $header $query awkfunctions]
			append awk " $query "
		}
		append awk " \{print [join $qposs ,]\}"
	}
	set pipe {}
	if {$sort ne ""} {
		lappend pipe $sort
	}
	if {[lindex $awkfunctions 0]} {lappend awkfunctions [tsv_select_min [lindex $awkfunctions 0]]}
	if {[lindex $awkfunctions 1]} {lappend awkfunctions [tsv_select_max [lindex $awkfunctions 1]]}
	set awk [join [list_remdup [lrange $awkfunctions 2 end]] \n]\n$awk
# putslog stderr ----------\n$awk\n----------
	if {[string trim $awk] ne ""} {
		lappend pipe [list awk $awk]
	}
#putslog -------------pipe-------------------
#putslog pipe:[join $pipe " | "]
#putslog ------------------------------------
	if {$qfields ne ""} {
		set nh $qfields
	} else {
		set nh $header
	}
	if {[llength $newheader]} {
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

if 0 {
	lappend auto_path ~/dev/completegenomics/lib
	package require Tclx
	signal -restart error SIGINT
	package require Extral
	cg select -h /complgen/multicompar/compar.tsv
	cg select -q 'same(GS102,GS103)' -f 'chromosome begin end reference type alleleSeq1-GS102 alleleSeq2-GS102 alleleSeq1-GS103 alleleSeq2-GS103' /complgen/multicompar/compar.tsv | less
	cg select -q 'count($coverage-GS102,$coverage-GS103,>20) == 2' -f 'chromosome begin end reference type alleleSeq1-GS102 alleleSeq2-GS102 alleleSeq1-GS103 alleleSeq2-GS103' /complgen/multicompar/compar.tsv | less

	set f [open GS102/ASM/var-GS000000078-ASM.tsv]
	set query "\$begin < 2000"
	set qfields "chromosome begin end"
	set sortfields "haplotype"
	set sortfields "chromosome begin"

cg select -q '' < GS102/ASM/var-GS000000078-ASM.tsv | less
cg select -f 'haplotype chromosome begin' < GS102/ASM/var-GS000000078-ASM.tsv | less
cg select -q '$begin < 2000' < GS102/ASM/var-GS000000078-ASM.tsv | less

cg select -q '$begin < 2000' < GS102/ASM/var-GS000000078-ASM.tsv > /tmp/test1
cg select -q '$begin < 2000' -f 'chromosome begin end' -s 'chromosome begin' < GS102/ASM/var-GS000000078-ASM.tsv > /tmp/test2
cg select -q '$begin < 2000' -f 'chromosome begin end' -s 'haplotype' < GS102/ASM/var-GS000000078-ASM.tsv > /tmp/test3

}

proc tsv_sort {filename fields} {
	set f [open $filename]
	set line [gets $f]
	close $f
	if {[string index $line 0] eq "#"} {set line [string range $line 1 end]}
	set header [split $line \t]
	set poss [list_cor $header $fields]
	if {[lsearch $poss -1] != -1} {error "fields [join [list_sub $fields [list_find $poss -1]] ,] not found"}
	set poss [lmath_calc $poss + 1]
	puts [join $header \t]
	set command "tail +2 [list $filename] | gnusort8 -T \"[scratchdir]\" -t \\t -N -s -k[join $poss " -k"] >@ stdout"
	eval exec $command
}

proc tsv_hcheader {f keepheaderVar headerVar} {
	upvar $keepheaderVar keepheader
	upvar $headerVar header
	set ::filebuffer($f) [list [join $header \t]]
	set temp [split [string trimright $keepheader] \n]
	set header [split [string range [list_pop temp] 1 end] \t]
	set keepheader [join $temp \n]\n
}

if 0 {
	lappend auto_path ~/dev/completegenomics/lib
	package require Tclx
	signal -restart error SIGINT
	package require Extral
	cd /complgen/compar

	set filename /data/db/ucsc_ori/_data_db_ucsc-exapted_repeats.tsv
	set f [open $filename]
	set fields {chrom chromStart chromEnd}

	# set f [open /complgen/compar/78vs79_compar-filter-sc.tsv]
	set f [open /complgen/compar/78vs79_compar_pvt.tsv]
	set query {compar df sample "|79 78,79" refcons "" ns "" lowscore "" trf "" str "" rp "" sd "" sc "" dbsnp "" loc EXON}

}

proc cg_select {args} {
	if {[llength $args] == 0} {
		errorformat select
		exit 1
	}
	set query {}; set fields {}; set sortfields {}; set newheader {}; set hc 0; set inverse 0
	set pos 0
	foreach {key value} $args {
		switch -- $key {
			-q {
				set query $value
				if {[regexp {[^=!><\\]=[^=]} $query]} {puts stderr "you may have used = instead of == in query"}
				regsub -all {\\=} $query = query
			}
			-f {set fields $value}
			-rf {
				set fields $value
				set inverse 1
			}
			-nh {set newheader $value}
			-hc {set hc 1}
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
	set error [catch {tsv_select $query $fields $sortfields $newheader $f $o $hc $inverse} result]
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

if 0 {

	lappend auto_path ~/dev/completegenomics/lib
	lappend auto_path /complgen/bin/complgen/apps/cg/lib
	lappend auto_path ~/bin/complgen/apps/cg/lib
	package require Extral
	package require Tclx
	signal -restart error SIGINT

	cd /complgen/1.8/cnvcg
	set file1 cnvcg-GS102.tsv
	set file2 cnvcg-GS103.tsv
	set joinfields1 {chr begin end}
	set joinfields2 {chr begin end}
	set postfix1 -GS102
	set postfix2 -GS103
	set method overlap

}

