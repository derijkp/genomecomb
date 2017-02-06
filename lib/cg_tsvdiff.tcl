proc tsvdiff_file {file1 file2 rcomments} {
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
		append error "header diff\n"
		append error "<extrafields: [list_lremove $header1 $common]\n"
		append error "---\n"
		append error ">extrafields: [list_lremove $header2 $common]\n"
	}
	set temp1 [tempfile]
	set temp2 [tempfile]
	if {[catch {
		cg select -rc $rcomments -f $common $file1 $temp1
	}]} {
		if {[catch {exec diff $file1 $file2} msg]} {
			puts stderr $msg
			return
		}
	}
	if {[catch {
		cg select -rc $rcomments -f $common $file2 $temp2
	}]} {
		if {[catch {exec diff $file1 $file2} msg]} {
			puts stderr $msg
			return
		}
	}
	if {[catch {exec diff $temp1 $temp2} result]} {
		append error "header\n  [join $common \t]\n"
		append error $result
	}
	if {$error ne ""} {
		puts stderr "diff $file1 $file2"
		puts stderr $error
	}
}

proc tsvdiff_file_brief {file1 file2 rcomments} {
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

proc tsvdiff {file1 file2 rcomments exclude brief} {
	if {![file isdir $file1] && ![file isdir $file1]} {
		if {[file extension [gzroot $file1]] in ".tsv .sft .hsmetrics"} {
			if {$brief} {
				tsvdiff_file_brief $file1 $file2 $rcomments
			} else {
				tsvdiff_file $file1 $file2 $rcomments
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
		tsvdiff [gzfile $file1/$file] [gzfile $file2/$file] $rcomments $exclude $brief
	}
}

proc cg_tsvdiff args {
	set rcomments 1
	set exclude {}
	set brief 0
	cg_options tsvdiff args {
		-c {
			if {[true $value]} {set rcomments 0} else {set rcomments 1}
		}
		-x - --exclude {
			lappend exclude $value
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
	tsvdiff $file1 $file2 $rcomments $exclude $brief
}

if 0 {
set file1 tmp/mastr_120477_MT_private/ceph1347_02_36640/sreg-gatk-crsbwa-ceph1347_02_36640.tsv.lz4
set file2 expected/mastr_120477_MT_private/ceph1347_02_36640/sreg-gatk-crsbwa-ceph1347_02_36640.tsv
cg tsvdiff tmp/mastr_120477_MT_private/ceph1347_02_36640/sreg-gatk-crsbwa-ceph1347_02_36640.tsv.lz4 expected/mastr_120477_MT_private/ceph1347_02_36640/sreg-gatk-crsbwa-ceph1347_02_36640.tsv

set dir1 tmp/mastr_120477_MT_private/ceph1347_02_36640
set dir2 expected/mastr_120477_MT_private/ceph1347_02_36640


}