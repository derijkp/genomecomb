proc xlong {file resultfile} {
	set o [open $resultfile.temp w]
	set f [gzopen $file]
	set header [tsv_open $f comment1]
	set poss [tsv_basicfields $header 6 0]
	set poss [list_remove $poss -1]
	set prefields [list_sub $header $poss]
	set rest [list_lremove [list_fill [llength $header] 0 1] $poss]
	set restfields [list_sub $header $rest]
	while {[gets $f line] != -1} {
		set line [split $line \t]
		set pre [list_sub $line $poss]
		set pre [join $pre -]
		foreach field $restfields value [list_sub $line $rest] {
			puts $o "$pre\t$field\t$value"
		}
	}
	gzclose $f
	close $o
	file rename $resultfile.temp $resultfile
	return 0
}

proc tsvdiff_file {file1 file2 rcomments type fields diffopts splitlines diffprog} {
	if {![catch {exec diff -q $file1 $file2}]} return
	set f1 [gzopen $file1]
	set header1 [tsv_open $f1 comment1]
	catch {close $f1}
	set f2 [gzopen $file2]
	set header2 [tsv_open $f2 comment2]
	catch {close $f2}
	set common [list_common $header1 $header2]
	if {[llength $fields]} {
		set h1 [list_common $header1 [expandfields $header1 $fields]]
		set h2 [list_common $header2 [expandfields $header2 $fields]]
	} else {
		set h1 $header1
		set h2 $header2
	}
	set common [list_common $h1 $h2]
	set error {}
	if {[llength $common] != [llength $h1] || [llength $common] != [llength $h2]} {
		append error "header diff\n"
		append error "<extrafields: [list_lremove $h1 $common]\n"
		append error "---\n"
		append error ">extrafields: [list_lremove $h2 $common]\n"
	}
	set temp1 [tempfile]
	set temp2 [tempfile]
	if {$type eq "xl"} {
		set error1 [xlong $file1 $temp1]
		set error2 [xlong $file2 $temp2]
	} elseif {$splitlines} {
		set error1 [catch {
			exec cg select -rc $rcomments -f $common $file1 | awk {{print $0"\n----"}} > $temp1
		}]
		set error2 [catch {
			exec cg select -rc $rcomments -f $common $file2 | awk {{print $0"\n----"}} > $temp2
		}]
	} else {
		set error1 [catch {
			cg select -rc $rcomments -f $common $file1 $temp1
		}]
		set error2 [catch {
			cg select -rc $rcomments -f $common $file2 $temp2
		}]
	}
	if {$error1 || $error2} {
		if {[catch {exec diff $file1 $file2} msg]} {
			puts stderr $msg
			return
		}
	}
	if {$diffprog ne ""} {
		exec {*}$diffprog $temp1 $temp2
	} else {
		if {[catch {exec diff {*}$diffopts $temp1 $temp2} result]} {
			append error "header\n  [join $common \t]\n"
			append error $result
		}
		if {$error ne ""} {
			puts stderr "diff $file1 $file2"
			puts stderr $error
		}
	}
}

proc tsvdiff_file_brief {file1 file2 rcomments type fields diffopts splitlines diffprog} {
	if {![catch {exec diff -q $file1 $file2}]} return
	set f1 [gzopen $file1]
	set header1 [tsv_open $f1 comment1]
	catch {close $f1}
	set f2 [gzopen $file2]
	set header2 [tsv_open $f2 comment2]
	catch {close $f2}
	set common [list_common $header1 $header2]
	set error {}
	if {[llength $common] != [llength $header1] || [llength $common] != [llength $header2]} {
		puts stderr "Files differ: $file1 $file2"; return
	}
	set temp1 [tempfile]
	set temp2 [tempfile]
	if {[catch {
		cg select -rc $rcomments -f $common $file1 $temp1
	}]} {
		if {[catch {exec diff -q $file1 $file2} msg]} {
			puts stderr "Files differ: $file1 $file2"; return
			return
		}
	}
	if {[catch {
		cg select -rc $rcomments -f $common $file2 $temp2
	}]} {
		if {[catch {exec diff -q $file1 $file2} msg]} {
			puts stderr "Files differ: $file1 $file2"; return
			return
		}
	}
	if {[catch {exec diff -q $temp1 $temp2} result]} {
		puts stderr "Files differ: $file1 $file2"; return
	}
}

proc tsvdiff {file1 file2 rcomments exclude brief type fields diffopts splitlines diffprog} {
	if {![file isdir $file1] && ![file isdir $file1]} {
		if {[file extension [gzroot $file1]] in ".tsv .sft .hsmetrics"} {
			if {$brief} {
				tsvdiff_file_brief $file1 $file2 $rcomments $type $fields $diffopts $splitlines $diffprog
			} else {
				tsvdiff_file $file1 $file2 $rcomments $type $fields $diffopts $splitlines $diffprog
			}
		} else {
			if {$brief} {
				if {[catch {exec diff -q $file1 $file2} msg]} {
					puts stderr "Files differ: $file1 $file2"
				}
			} else {
				if {[catch {exec diff $file1 $file2} msg]} {
					puts stderr "diff $file1 $file2"
					puts stderr $msg
				}
			}
		}
		return
	}
	if {![file isdir $file1] || ![file isdir $file1]} {
		puts stderr "Dir vs file: $file1 $file2"
	}
	set files1 [lsort -dict [dirglob $file1 *]]
	foreach pattern $exclude {
		set files1 [list_sub $files1 -exclude [list_find -regexp $files1 $pattern]]
	}
	set files2 [lsort -dict [dirglob $file2 *]]
	foreach pattern $exclude {
		set files2 [list_sub $files2 -exclude [list_find -regexp $files2 $pattern]]
	}
	set gzfiles1 {}
	foreach file $files1 {
		lappend gzfiles1 [gzroot $file]
	}
	set gzfiles2 {}
	foreach file $files2 {
		lappend gzfiles2 [gzroot $file]
	}
	set common [list_common $gzfiles1 $gzfiles2]
	foreach file [list_lremove $gzfiles1 $common] {
		puts stderr "Only in $file1: [file root [gzfile $file1/$file]]"
	}
	foreach file [list_lremove $gzfiles2 $common] {
		puts stderr "Only in $file2: [file root [gzfile $file2/$file]]"
	}
	foreach file $common {
		tsvdiff [gzfile $file1/$file] [gzfile $file2/$file] $rcomments $exclude $brief $type $fields $diffopts $splitlines $diffprog
	}
}

proc cg_tsvdiff args {
	set rcomments 1
	set fields {}
	set exclude {}
	set brief 0
	set type diff
	set splitlines 0
	set diffprog {}
	set diffopts {}
	cg_options tsvdiff args {
		-c {
			if {[true $value]} {set rcomments 0} else {set rcomments 1}
		}
		-f - --fields {
			set fields $value
		}
		-x - --exclude {
			lappend exclude $value
		}
		-t - --type {
			set type $value
		}
		-y - --side-by-side {
			set side-by-side $value
		}
		--suppress-common-lines {
			set suppress-common-lines $value
		}
		-s - --splitlines {
			set splitlines $value
		}
		-d - --diffprog {
			set diffprog $value
		}
		-w - --width {
			set width $value
		}
		-q - --brief {
			if {$value in "1 0"} {
				set brief $value
			} else {
				set brief 1
				incr pos -1
			}
		}
	} {file1 file2}
	if {$type eq "xl"} {
		if {![info exists side-by-side]} {set side-by-side 1}
		if {![info exists suppress-common-lines]} {set suppress-common-lines 1}
		if {![info exists width]} {set width 200}
	} elseif {$type eq "sd"} {
		set side-by-side 1
		set suppress-common-lines 1
		if {![info exists width]} {set width 200}
	}
	if {[get side-by-side 0]} {lappend diffopts --side-by-side}
	if {[get suppress-common-lines 0]} {lappend diffopts --suppress-common-lines}
	if {[info exists width]} {lappend diffopts --width=$width}
	tsvdiff $file1 $file2 $rcomments $exclude $brief $type $fields $diffopts $splitlines $diffprog
}

if 0 {
set file1 tmp/mastr_120477_MT_private/ceph1347_02_36640/sreg-gatk-crsbwa-ceph1347_02_36640.tsv.lz4
set file2 expected/mastr_120477_MT_private/ceph1347_02_36640/sreg-gatk-crsbwa-ceph1347_02_36640.tsv
cg tsvdiff tmp/mastr_120477_MT_private/ceph1347_02_36640/sreg-gatk-crsbwa-ceph1347_02_36640.tsv.lz4 expected/mastr_120477_MT_private/ceph1347_02_36640/sreg-gatk-crsbwa-ceph1347_02_36640.tsv

set dir1 tmp/mastr_120477_MT_private/ceph1347_02_36640
set dir2 expected/mastr_120477_MT_private/ceph1347_02_36640


}