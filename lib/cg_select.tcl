# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

proc tsv_select_idtopos {header id fields} {
	set poss [list_cor $header $fields]
	if {[inlist $poss -1]} {error "sample $id not present"}
	return [lmath_calc $poss + 1]
}

proc tsv_select_sm {header ids} {
	set id1 [list_pop ids]
	foreach {a11 a21 sequenced} [tsv_select_idtopos $header $id1 [list alleleSeq1-$id1 alleleSeq2-$id1 sequenced-$id1]] break
	set temp [list "(\$$sequenced == \"v\")"]
	foreach id $ids {
		foreach {a12 a22 sequenced} [tsv_select_idtopos $header $id1 [list alleleSeq1-$id alleleSeq2-$id sequenced-$id]] break
		lappend temp "(\$$sequenced == \"v\")"  "((\$$a11 == \$$a12 && \$$a21 == \$$a22) || (\$$a11 == \$$a22 && \$$a21 == \$$a12))"
	}
	set temp "([join $temp " && "])"
}

proc tsv_select_same {header ids} {
	set id1 [list_pop ids]
	foreach {a11 a21 sequenced} [tsv_select_idtopos $header $id1 [list alleleSeq1-$id1 alleleSeq2-$id1 sequenced-$id1]] break
	set seqlist [list "\$$sequenced != \"u\""]
	set temp {}
	foreach id $ids {
		lappend seqlist "\$$sequenced != \"u\""
		foreach {a12 a22 sequenced} [tsv_select_idtopos $header $id1 [list alleleSeq1-$id alleleSeq2-$id sequenced-$id]] break
		lappend temp "((\$$a11 == \$$a12 && \$$a21 == \$$a22) || (\$$a11 == \$$a22 && \$$a21 == \$$a12))"
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
			lappend temp2 "\$$a1 != \$$a12 || \$$a2 != \$$a22"
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

proc tsv_select {query {qfields {}} {sortfields {}} {newheader {}} {f stdin} {out stdout}} {
	fconfigure $f -buffering none
	fconfigure $out -buffering none
	set header [tsv_open $f]
	set awk ""
	set awkfunctions {0 0}
	set sort ""
	set cut ""
	set qfields [tsv_select_expandfields $header $qfields qposs awkfunctions]
	set tsv_funcnum 1
	if {[llength $sortfields]} {
		set poss [list_cor $header $sortfields]
		if {[lsearch $poss -1] != -1} {error "fields [join [list_sub $sortfields [list_find $poss -1]] ,] not found"}
		set poss [lmath_calc $poss + 1]
		set keys {}
		foreach pos $poss {
			lappend keys $pos,$pos
		}
		set sort "gnusort8 -t \\t -V -s -k[join $keys " -k"]"
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
# putslog stderr ----------\n$awk\n----------
	if {[lindex $awkfunctions 0]} {lappend awkfunctions [tsv_select_min [lindex $awkfunctions 0]]}
	if {[lindex $awkfunctions 1]} {lappend awkfunctions [tsv_select_max [lindex $awkfunctions 1]]}
	set awk [join [list_remdup [lrange $awkfunctions 2 end]] \n]\n$awk
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
		puts $out [join $newheader \t]
	} else	{
		puts $out [join $nh \t]
	}
	if {![llength $pipe]} {
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
	set command "tail +2 [list $filename] | gnusort8 -t \\t -V -s -k[join $poss " -k"] >@ stdout"
	eval exec $command
}

proc file_rootgz {filename} {
	set ext [file extension $filename]
	if {$ext eq ".gz"} {
		return [file root [file root $filename]]
	} elseif {$ext eq ".rz"} {
		return [file root [file root $filename]]
	} elseif {$ext eq ".bgz"} {
		return [file root [file root $filename]]
	} else {
		file root [file root $filename]
	}
}

proc tsv_open {f} {
	set keep 0
	set buffering [fconfigure $f -buffering]
	fconfigure $f -buffering line
	set vcf 0
	set first 1
	while {![eof $f]} {
		set line [gets $f]
		if {![string length $line]} continue
		set fchar [string index $line 0]
		if {$fchar ne "#"} break
		if {$first} {
			set fchar2 [string index $line 1]
			if {$fchar2 eq "#"} {set vcf 1}
		} elseif {$vcf} {
			if {$fchar2 ne "#"} break
		}
		set first 0
		
		set keep [tell $f]
		set header $line
	}
	fconfigure $f -buffering $buffering
	if {[string index $line 0] eq ">"} {
		return [split [string range $line 1 end] \t]
	}
	if {![info exists header]} {
		return [split $line \t]
	} elseif {[string index $header 0] eq "#"} {
		seek $f $keep
		return [split [string range $header 1 end] \t]
	} else {
		seek $f $keep
		return [split $header \t]
	}
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

proc tsv_next {f xpos next {shift 100000}} {
	# do a ~ binary search to get at next faster
	set start [tell $f]
	while 1 {
		seek $f $shift current
		gets $f
		set line [getnotempty $f]
		set x [lindex $line $xpos]
		if {![isdouble $x] || ($x >= $next)} {
			seek $f $start
			set shift [expr {$shift / 2}]
			if {$shift < 1000} break
		}
		set start [tell $f]
	}
	while {![eof $f]} {
		set fpos [tell $f]
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set x [lindex $line $xpos]
		if {$x >= $next} break
	}
	if {![eof $f]} {
		return $fpos
	} else {
		return $x
	}
}

proc tsv_nextline {f xpos next {shift 100000}} {
	# do a ~ binary search to get at next faster
	set start [tell $f]
	while 1 {
		seek $f $shift current
		gets $f
		set line [getnotempty $f]
		set x [lindex $line $xpos]
		if {![isdouble $x] || ($x >= $next)} {
			seek $f $start
			set shift [expr {$shift / 2}]
			if {$shift < 1000} break
		}
		set start [tell $f]
	}
	while {![eof $f]} {
		set fpos [tell $f]
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set x [lindex $line $xpos]
		if {$x >= $next} break
	}
	return $line
}

proc tsv_index {xfield file} {
	set indexname [gzroot $file].${xfield}_index
	if {[inlist {.rz} [file extension $file]]} {
		set tempfile [scratchfile]
		exec razip -d -c $file > $tempfile
		set f [open $tempfile]
	} elseif {[inlist {.bgz} [file extension $file]]} {
		set tempfile [scratchfile]
		exec bgzip -d -c $file > $tempfile
		set f [open $tempfile]
	} elseif {[inlist {.gz} [file extension $file]]} {
		set tempfile [scratchfile]
		exec gunzip -c $file > $tempfile
		set f [open $tempfile]
	} else {
		set f [open $file]
	}
	set header [tsv_open $f]
	set xpos [lsearch $header $xfield]
	if {$xpos == -1} {error "field $xfield not present in file $file"}
	set fstart [tell $f]
	set fpos $fstart
	set line [split [gets $f] \t]
	set xmin [lindex $line $xpos]
	set findex [expr {$xmin-$xmin%10000}]
	set prev $findex
	set next [expr {$prev + 10000}]
	set index [list $fpos]
	catch {Classy::Progress start [file size $file] "Making index"}
	while {![eof $f]} {
		set fpos [tsv_next $f $xpos $next]
		if {[eof $f]} {
			set xmax $fpos
			break
		}
		lappend index $fpos
		incr prev 10000
		incr next 10000
		catch {Classy::Progress set [tell $f]}
		if {![expr $next%1000000]} {putslog $next}
	}
	catch {Classy::Progress stop}
	close $f
	set size [file size $file]
	set f [open $file]
	seek $f [expr {$size-5000}] start
	gets $f
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set xmax [lindex $line $xpos]
	}
	close $f
	set o [open $indexname.temp w]
	puts $o 10000
	puts $o $findex
	puts $o $xmin
	puts $o $xmax
	puts $o [join $index \n]
	close $o
	file rename $indexname.temp $indexname
	if {[info exists tempfile]} {
		file delete $tempfile
	}
}

proc tsv_index_header {file} {
	global cache
	set file [file normalize $file]
	return $cache(tsv_index,$file,header)
}

proc tsv_index_open {file field {uncompress 0}} {
	global cache
	set file [file normalize $file]
	if {[info exists cache(tsv_index,$file,$field,step)]} return
	set ext [file extension $file]
	if {$ext eq ".gz" || $ext eq ".rz" || $ext eq ".bgz"} {set uncompress 1}
	set root [gzroot $file]
	set workfile $file
	set uncompressed 0
	set remove 0
	if {[inlist {.rz .gz .bgz} $ext]} {
		if {$uncompress} {
			set workfile [scratchfile]
			puts stderr "temporarily uncompressing $file to $workfile"
			file delete -force $workfile.temp
			file delete -force $workfile
			gunzip $file $workfile.temp
			file rename $workfile.temp $workfile
			set uncompressed 1
			set remove 1
		}
	} else {
		set uncompressed 1
	}
	set indexname $root.${field}_index
	if {![file exists $indexname]} {
		tsv_index $field $workfile
		if {$workfile ne $file} {
			file rename [gzroot $workfile].${field}_index $indexname
		}
	}
	set o [open $indexname]
	set cache(tsv_index,$file,$field,step) [gets $o]
	set cache(tsv_index,$file,$field,findex) [gets $o]
	set xmin [gets $o]
	set cache(tsv_index,$file,$field,xmin) $xmin
	set cache(tsv_index,$file,$field,xmax) [gets $o]
	set cache(tsv_index,$file,$field,fx) [expr {$xmin-$xmin%10000}]
	set cache(tsv_index,$file,$field,index) [split [string trim [read $o]] \n]
	set cache(tsv_index,$file,$field,workfile) $workfile
	set cache(tsv_index,$file,$field,uncompressed) $uncompressed
	set cache(tsv_index,$file,$field,remove) $remove
	close $o
	set f [gzopen $workfile]
	set cache(tsv_index,$file,header) [tsv_open $f]
	if {$uncompressed} {
		set cache(tsv_index,$file,$field,channel) $f
	} else {
		catch {close $f}
	}
}

proc tsv_index_close {file field} {
	global cache
	set file [file normalize $file]
	if {[info exists cache(tsv_index,$file,$field,workfile)] && [get cache(tsv_index,$file,$field,remove) 0] && ($cache(tsv_index,$file,$field,workfile) ne $file)} {
		close $cache(tsv_index,$file,$field,channel)
		# puts "remove $cache(tsv_index,$file,$field,workfile)"
		file delete $cache(tsv_index,$file,$field,workfile)
	}
	set indexname [gzroot $file].${field}_index
	unset -nocomplain cache(tsv_index,$file,$field,step)
	unset -nocomplain cache(tsv_index,$file,$field,findex)
	unset -nocomplain cache(tsv_index,$file,$field,xmin)
	unset -nocomplain cache(tsv_index,$file,$field,xmax)
	unset -nocomplain cache(tsv_index,$file,$field,fx)
	unset -nocomplain cache(tsv_index,$file,$field,index)
	unset -nocomplain cache(tsv_index,$file,$field,workfile)
	unset -nocomplain cache(tsv_index,$file,$field,uncompressed)
	unset -nocomplain cache(tsv_index,$file,$field,remove)
	unset -nocomplain cache(tsv_index,$file,$field,channel)
	unset -nocomplain cache(tsv_index,$file,header)
}

proc tsv_index_apprstop {file field} {
	global cache
	set uncompressed $cache(tsv_index,$file,$field,uncompressed)
	if {!$uncompressed} {
		catch {close $f}
	}
}

proc tsv_index_apprgoto {file field pos} {
	global cache
	set file [file normalize $file]
	set index $cache(tsv_index,$file,$field,index)
	set step $cache(tsv_index,$file,$field,step)
	set start [expr {round($pos)-round($pos)%$step}]
	set indexpos [expr {($start-$cache(tsv_index,$file,$field,findex))/$step}]
	if {$indexpos < 0} {set indexpos 0}
	if {$indexpos >= [llength $index]} {set indexpos end}
	set fpos [expr {round([lindex $index $indexpos])}]
	set uncompressed $cache(tsv_index,$file,$field,uncompressed)
	if {$uncompressed} {
		set f $cache(tsv_index,$file,$field,channel)
		seek $f $fpos
	} else {
		set f [gzopen $file $fpos]
	}
	return $f
}

proc tsv_index_get {file field pos} {
	global cache
	set file [file normalize $file]
	set header $cache(tsv_index,$file,header)
	set fieldpos [lsearch $header $field]
	set uncompressed [get cache(tsv_index,$file,$field,uncompressed) 0]
	set f [tsv_index_apprgoto $file $field $pos]
	if {$uncompressed} {
		set line [tsv_nextline $f $fieldpos $pos]
	} else {
		set line {}
		while 1 {
			# read in mem in chunks
			# do a ~ binary search to get at target faster
			set chunk [read $f 20480]
			append chunk [gets $f]
			if {![string length $chunk]} break
			set table [split $chunk \n]
			set temp [lindex $table end $fieldpos]
			if {$temp < $pos} continue
			if {$temp == $pos} break
			set ipos [binsearch $table $fieldpos $pos]
			set line [lindex $table $ipos]
			break
		}
	}
	tsv_index_apprstop $file $field
	set temp [lindex $line $fieldpos]
	if {$temp != $pos} {error "$pos not found in $file,$field"}
	return $line
}

# -1 if comp1 < comp2
# -2 if comp1 > comp2
# 0 if comp1 = comp2
# > 0 if comp1 ~ comp2
proc tsv_align_match {comp1 comp2} {
	foreach e1 $comp1 e2 $comp2 {
		if {$e1 ne $e2} {
			if {$e2 eq ""} {
				return -1
			} elseif {$e1 eq ""} {
				return -2
			}
			if {[lsort -dict [list $e1 $e2]] eq [list $e1 $e2]} {return -1} else {return -2}
		}
	}
	return 0
}

proc tsv_align_compareoverlap {comp1 comp2} {
	foreach {chr1 p11 p12} $comp1 break
	foreach {chr2 p21 p22} $comp2 break
	if {$chr1 ne $chr2} {
		if {$chr2 eq ""} {
			return -1
		} elseif {$chr1 eq ""} {
			return -2
		}
		if {[lsort -dict [list $chr1 $chr2]] eq [list $chr1 $chr2]} {return -1} else {return -2}
	}
	if {![isint $p11]} {
		return -1
	} elseif {![isint $p21]} {
		return -2
	} elseif {$p12 <= $p21} {
		return -1
	} elseif {$p22 <= $p11} {
		return -2
	} else {
		return 0
		return [overlap $p11 $p12 $p21 $p22]
	}
}

proc tsv_align {file1 file2 joinfields1 joinfields2 postfix1 postfix2} {
	catch {close $f1}; catch {close $f2}
	set o stdout
	# open files
	set f1 [open $file1]
	set header1 [split [gets $f1] \t]
	set comparposs1 [list_cor $header1 $joinfields1]
	if {[inlist $comparposs1 -1]} {
		error "$file1 does not contain all joinfields"
	}
	set f2 [open $file2]
	set header2 [split [gets $f2] \t]
	set comparposs2 [list_cor $header2 $joinfields2]
	if {[inlist $comparposs2 -1]} {
		error "$file2 does not contain all joinfields"
	}
	#
	set dummy1 [list_fill [llength $header1] {}]
	set dummy2 [list_fill [llength $header2] {}]
	# start making output files
	# make output header
	set oheader {}
	foreach field $header1 {
		lappend oheader ${field}$postfix1
	}
	# make output header
	foreach field $header2 {
		lappend oheader ${field}$postfix2
	}
	puts $o [join $oheader \t]
	set cur1 [split [gets $f1] \t]
	set comp1 [list_sub $cur1 $comparposs1]
	set cur2 [split [gets $f2] \t]
	set comp2 [list_sub $cur2 $comparposs2]
	set prevcomp1 {}
	set prevcomp2 {}
	# go over the files
	set num 1
	while {![eof $f1] || ![eof $f2]} {
		incr num
		if {![expr {$num % 100000}]} {putslog $num}
		set d [tsv_align_match $comp1 $comp2]
		if {$d == 0} {
			puts $o [join $cur1 \t]\t[join $cur2 \t]
			set cur1 [compare_annot_getline $f1]
			set comp1 [list_sub $cur1 $comparposs1]
			set cur2 [compare_annot_getline $f2]
			set comp2 [list_sub $cur2 $comparposs2]
		} elseif {$d == -1} {
			while {[tsv_align_match $comp1 $comp2] == -1} {
				puts $o [join $cur1 \t]\t[join $dummy2 \t]
				if {[eof $f1]} break
				set cur1 [compare_annot_getline $f1]
				set comp1 [list_sub $cur1 $comparposs1]
				if {![llength $cur1]} break
				incr num
				if {![expr {$num % 100000}]} {putslog $num}
			}
		} elseif {$d == -2} {
			while {[tsv_align_match $comp1 $comp2] == -2} {
				puts $o [join $dummy1 \t]\t[join $cur2 \t]
				if {[eof $f2]} break
				set cur2 [compare_annot_getline $f2]
				set comp2 [list_sub $cur2 $comparposs2]
				if {![llength $cur2]} break
				incr num
				if {![expr {$num % 100000}]} {putslog $num}
			}
		}
	}
	close $f1; close $f2
}

proc cg_select {args} {
	if {[llength $args] == 0} {
		errorformat select
		exit 1
	}
	set query {}; set fields {}; set sortfields {}; set newheader {}
	set pos 0
	foreach {key value} $args {
		switch -- $key {
			-q {
				set query $value
				if {[regexp {[^=!><]=[^=]} $query]} {puts stderr "you used = instead of == in query"; exit 1}
			}
			-f {set fields $value}
			-nh {set newheader $value}
			-s {set sortfields $value}
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
	set error [catch {tsv_select $query $fields $sortfields $newheader $f $o} result]
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

set params {compar-cnvcg.tsv ogtregions.tsv "chr-GS102 begin-GS102 end-GS102" "chromosome begin end" "" "-ogt" overlap}
set params {cnvcg-GS102.tsv cnvcg-GS103.tsv {chr begin end} {chr begin end} -GS102 -GS103 overlap}
set params {allvalsnps3.tsv ../multicompar/compar.tsv {chromosome begin end type} {chromosome begin end type} -val {}}
foreach {file1 file2 joinfields1 joinfields2 postfix1 postfix2 method} $params break

	tsv_align $file1 $file2 $joinfields1 $joinfields2 $postfix1 $postfix2 $method
	
paste acnvcg-GS102.tsv acnvcg-GS103.tsv acnvcg-GS102.tsv > cnvcg-compar.tsv

}

