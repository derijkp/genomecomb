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
	global errors
	set f1 [gzopen $file1]
	set header1 [tsv_open $f1 comment1]
	catch {close $f1}
	set f2 [gzopen $file2]
	set header2 [tsv_open $f2 comment2]
	catch {close $f2}
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
	set tempdir [tempdir]
	set temp1 $tempdir/[file tail $file1]
	set temp2 $tempdir/[file tail $file2]
	if {$temp2 eq $temp1} {append temp2 -b}
	if {$type eq "xl"} {
		set error1 [xlong $file1 $temp1]
		set msg1 "error in xl conversion"
		set error2 [xlong $file2 $temp2]
		set msg2 "error in xl conversion"
	} elseif {$splitlines} {
		set error1 [catch {
			exec cg select -rc $rcomments -f $common $file1 | awk {{print $0"\n----"}} > $temp1
		} msg1]
		set error2 [catch {
			exec cg select -rc $rcomments -f $common $file2 | awk {{print $0"\n----"}} > $temp2
		} msg2]
	} else {
		set error1 [catch {
			cg select -rc $rcomments -f $common $file1 $temp1
		} msg1]
		set error2 [catch {
			cg select -rc $rcomments -f $common $file2 $temp2
		} msg2]
	}
	if {$error1} {
		error "Could not convert file $file1 for comparison: $msg1"
	}
	if {$error2} {
		error "Could not convert file $file2 for comparison: $msg2"
	}
	if {$diffprog ne ""} {
		exec {*}$diffprog $temp1 $temp2
	} else {
		if {[catch {exec diff {*}$diffopts $temp1 $temp2} result]} {
			append error "header\n  [join $common \t]\n"
			append error $result
		}
		if {$error ne ""} {
			incr errors
			puts "diff $file1 $file2"
			puts $error
		}
	}
}

proc tsvdiff_file_brief {file1 file2 rcomments type fields diffopts splitlines diffprog} {
	global errors
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
		incr errors
		puts "Files differ: [string_change $file1 {{ } {\ }}] [string_change $file2 {{ } {\ }}]"
		return
	}
	set temp1 [tempfile]
	set temp2 [tempfile]
	if {[catch {
		cg select -rc $rcomments -f $common $file1 $temp1
	}]} {
		if {[catch {exec diff -q $file1 $file2} msg]} {
			incr errors
			puts "Files differ: [string_change $file1 {{ } {\ }}] [string_change $file2 {{ } {\ }}]"
			return
		}
	}
	if {[catch {
		cg select -rc $rcomments -f $common $file2 $temp2
	}]} {
		if {[catch {exec diff -q $file1 $file2} msg]} {
			incr errors
			puts "Files differ: [string_change $file1 {{ } {\ }}] [string_change $file2 {{ } {\ }}]"
			return
		}
	}
	if {[catch {exec diff -q $temp1 $temp2} result]} {
		incr errors
		puts "Files differ: [string_change $file1 {{ } {\ }}] [string_change $file2 {{ } {\ }}]"
		return
	}
}

proc cg_tsvdiff args {
	global errors
	set rcomments 1
	set fields {}
	set brief 0
	set type diff
	set splitlines 0
	set diffprog {}
	set diffopts {}
	set excludeopts {}
	cg_options tsvdiff args {
		-c {
			if {[true $value]} {set rcomments 0} else {set rcomments 1}
		}
		-f - --fields {
			set fields $value
		}
		-x - --exclude {
			lappend excludeopts -x $value
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
	} {dir1 dir2}
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
	set errors 0
	if {![catch {exec diff -r {*}$diffopts {*}$excludeopts --brief $dir1 $dir2} diff]} {exit 0}
	set len1 [string length $dir1]
	set last1 [expr {$len1 - 1}]
	set len2 [string length $dir2]
	set last2 [expr {$len2 - 1}]
	set diff [split [string trim $diff] \n]
	foreach line $diff {
		if {[regexp {^Only in (.*): (.*)$} $line temp dir file]} {
			set file $dir/$file
			if {[string range $dir 0 $last1] eq $dir1} {
				set post [string range $file $len1 end]
				set file2 [gzfile [gzroot $dir2$post]]
				if {[file exists $file2]} {
					if {!$brief} {
						tsvdiff_file $file $file2 $rcomments $type $fields $diffopts $splitlines $diffprog
					} else {
						tsvdiff_file_brief $file $file2 $rcomments $type $fields $diffopts $splitlines $diffprog
					}
				} else {
					incr errors
					puts "Only in 1: $file"
				}
			} else {
				set post [string range $file $len2 end]
				set file1 [gzfile [gzroot $dir1$post]]
				# skip if file exists (only compare if we come across it in 1)
				if {![file exists $file1]} {
					incr errors
					puts "Only in 2: $file"
				}
			}
		} elseif {[regexp {^Files (.*) and (.*) differ$} $line temp file1 file2]} {
			if {[file extension [gzroot $file1]] in ".tsv .sft .hsmetrics"} {
				if {$brief} {
					tsvdiff_file_brief $file1 $file2 $rcomments $type $fields $diffopts $splitlines $diffprog
				} else {
					tsvdiff_file $file1 $file2 $rcomments $type $fields $diffopts $splitlines $diffprog
				}
			} else {
				if {$brief} {
					incr errors
					puts "Files differ: [string_change $file1 {{ } {\ }}] [string_change $file2 {{ } {\ }}]"
				} else {
					if {[catch {exec diff $file1 $file2} msg]} {
						incr errors
						puts "diff $file1 $file2"
						puts $msg
					}
				}
			}
		} elseif {$line eq "child process exited abnormally"} {
		} else {
			incr errors
			puts $line
		}
	}
	if {$errors} {exit 1}
}

if 0 {
set file1 tmp/mastr_120477_MT_private/ceph1347_02_36640/sreg-gatk-crsbwa-ceph1347_02_36640.tsv.lz4
set file2 expected/mastr_120477_MT_private/ceph1347_02_36640/sreg-gatk-crsbwa-ceph1347_02_36640.tsv
cg tsvdiff tmp/mastr_120477_MT_private/ceph1347_02_36640/sreg-gatk-crsbwa-ceph1347_02_36640.tsv.lz4 expected/mastr_120477_MT_private/ceph1347_02_36640/sreg-gatk-crsbwa-ceph1347_02_36640.tsv

set dir1 tmp/mastr_120477_MT_private/ceph1347_02_36640
set dir2 expected/mastr_120477_MT_private/ceph1347_02_36640


}