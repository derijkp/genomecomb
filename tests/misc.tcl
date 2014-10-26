#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

set keepdir [pwd]

file mkdir tmp
file mkdir tmp/test
file_write tmp/test/test1 1
file_write tmp/test/test2 2
file mkdir tmp/test/subdir
file_write tmp/test/subdir/test3 3
file mkdir tmp/out

proc dirinfo {dir} {
	set result {}
	foreach file [lsort -dict [glob $dir/*]] {
		if {[catch {file link $file} link]} {set link -}
		if {[file isdir $file]} {
			set data [list dir [dirinfo $file]]
		} else {
			set data [file_read $file]
		}
		lappend result [list $file $link $data]
	}
	return $result
}

test cplinked {relative} {
	file delete -force tmp/out
	file mkdir tmp/out
	exec cg cplinked tmp/test tmp/out
	set result [dirinfo tmp/out/test]
	file delete -force tmp/out/test
	set result
} {{tmp/out/test/subdir - {dir {{tmp/out/test/subdir/test3 ../../../test/subdir/test3 3}}}} {tmp/out/test/test1 ../../test/test1 1} {tmp/out/test/test2 ../../test/test2 2}} 

test cplinked {over existing} {
	file delete -force tmp/out
	file mkdir tmp/out
	exec cg cplinked tmp/test tmp/out
	exec cg cplinked tmp/test tmp/out
	set result [dirinfo tmp/out/test]
	file delete -force tmp/out/test
	set result
} {{tmp/out/test/subdir - {dir {{tmp/out/test/subdir/test3 ../../../test/subdir/test3 3}}}} {tmp/out/test/test1 ../../test/test1 1} {tmp/out/test/test2 ../../test/test2 2}} 

test cplinked {relative, new dir} {
	file delete -force tmp/out
	file mkdir tmp/out
	exec cg cplinked tmp/test tmp/out/test2
	set result [dirinfo tmp/out/test2]
	file delete -force tmp/out/test2
	set result
} {{tmp/out/test2/subdir - {dir {{tmp/out/test2/subdir/test3 ../../../test/subdir/test3 3}}}} {tmp/out/test2/test1 ../../test/test1 1} {tmp/out/test2/test2 ../../test/test2 2}} 

test cplinked {relative, new dir 2} {
	file delete -force tmp/out
	exec cg cplinked tmp/test tmp/out/test
	set result [dirinfo tmp/out/test]
	file delete -force tmp/out/test
	set result
} {{tmp/out/test/subdir - {dir {{tmp/out/test/subdir/test3 ../../../test/subdir/test3 3}}}} {tmp/out/test/test1 ../../test/test1 1} {tmp/out/test/test2 ../../test/test2 2}} 

test cplinked {relative extra level} {
	file delete -force tmp/out
	file mkdir tmp/out/test2
	exec cg cplinked tmp/test tmp/out/test2/testl
	set result [dirinfo tmp/out/test2/testl]
	file delete -force tmp/out
	set result
} {{tmp/out/test2/testl/subdir - {dir {{tmp/out/test2/testl/subdir/test3 ../../../../test/subdir/test3 3}}}} {tmp/out/test2/testl/test1 ../../../test/test1 1} {tmp/out/test2/testl/test2 ../../../test/test2 2}} 

test cplinked {absolute} {
	file delete -force /tmp/test2
	exec cg cplinked tmp/test /tmp/test2
	set result [dirinfo /tmp/test2]
	file delete -force /tmp/test2
	set result
} [list [list /tmp/test2/subdir - [list dir [list [list /tmp/test2/subdir/test3 $keepdir/tmp/test/subdir/test3 3]]]] [list /tmp/test2/test1 $keepdir/tmp/test/test1 1] [list /tmp/test2/test2 $keepdir/tmp/test/test2 2]]

test cplinked {backup existing} {
	file delete -force tmp/out
	file mkdir tmp/out/test
	file_write tmp/out/test/test1 pre
	exec cg cplinked tmp/test tmp/out 2> /dev/null
	set result [dirinfo tmp/out/test]
	file delete -force tmp/out/test
	set result
} {{tmp/out/test/subdir - {dir {{tmp/out/test/subdir/test3 ../../../test/subdir/test3 3}}}} {tmp/out/test/test1 ../../test/test1 1} {tmp/out/test/test1.old1 - pre} {tmp/out/test/test2 ../../test/test2 2}} 

test cplinked {backup existing 2} {
	file delete -force tmp/out
	file mkdir tmp/out/test
	file_write tmp/out/test/test1 pre
	file_write tmp/out/test/test1.old1 old
	exec cg cplinked tmp/test tmp/out 2> /dev/null
	set result [dirinfo tmp/out/test]
	file delete -force tmp/out/test
	set result
} {{tmp/out/test/subdir - {dir {{tmp/out/test/subdir/test3 ../../../test/subdir/test3 3}}}} {tmp/out/test/test1 ../../test/test1 1} {tmp/out/test/test1.old1 - old} {tmp/out/test/test1.old2 - pre} {tmp/out/test/test2 ../../test/test2 2}} 

test cplinked {file} {
	file delete -force tmp/out
	file mkdir tmp/out
	exec cg cplinked tmp/test/test1 tmp/out/testl
	set result [dirinfo tmp/out/]
	file delete -force tmp/out/test
	set result
} {{tmp/out/testl ../test/test1 1}} 

test cplinked {file replace link} {
	file delete -force tmp/out
	file mkdir tmp/out
	exec cg cplinked tmp/test/test1 tmp/out/testl
	exec cg cplinked tmp/test/test2 tmp/out/testl
	set result [dirinfo tmp/out/]
	file delete -force tmp/out/test
	set result
} {{tmp/out/testl ../test/test2 2}} 

test cplinked {file replace file} {
	file delete -force tmp/out
	file mkdir tmp/out
	file_write tmp/out/testl pre
	exec cg cplinked tmp/test/test1 tmp/out/testl
	set result [dirinfo tmp/out/]
	file delete -force tmp/out/test
	set result
} {{tmp/out/testl ../test/test1 1} {tmp/out/testl.old1 - pre}} 

test distr2chr {basic} {
	test_cleantmp
	exec distr2chr tmp/distrvars1- 0 < data/vars1.sft
	list [lsort -dict [glob tmp/*]] [file_read tmp/distrvars1-chr1]
} {{tmp/distrvars1-chr1 tmp/distrvars1-chr2 tmp/distrvars1-chromosome} {chr1	4000	4001	snp	G	A	A	G	1	v	A	G	0	v	4
chr1	4001	4002	snp	A	G,C	G	G	1	v	G	G	0	v	1;2,3;4
chr1	4099	5000	snp	C	T	T	T	47	v	T	T	35	v	1,2
chr1	5000	5010	del	AGCGTGGCAA		AGCGTGGCAA		32	v	-	-	41	u	1;2
chr1	5020	5021	snp	G	C	G	C	54	v	G	G	52	r	3
}} 

test distr2chr {header} {
	test_cleantmp
	exec distr2chr tmp/distrvars1- 0 1 < data/vars1.sft
	list [lsort -dict [glob tmp/*]] [file_read tmp/distrvars1-chr1]
} {{tmp/distrvars1-chr1 tmp/distrvars1-chr2} {chr1	4000	4001	snp	G	A	A	G	1	v	A	G	0	v	4
chr1	4001	4002	snp	A	G,C	G	G	1	v	G	G	0	v	1;2,3;4
chr1	4099	5000	snp	C	T	T	T	47	v	T	T	35	v	1,2
chr1	5000	5010	del	AGCGTGGCAA		AGCGTGGCAA		32	v	-	-	41	u	1;2
chr1	5020	5021	snp	G	C	G	C	54	v	G	G	52	r	3
}} 

test distr2chr {many} {
	test_cleantmp
	exec distr2chr tmp/distrvars1- 0 < data/manychr.tsv
	exec cat {*}[glob tmp/distrvars1-*] > tmp/cat.tsv
	exec sort tmp/cat.tsv > tmp/scat.tsv
	exec sort data/manychr.tsv > tmp/check.tsv
	set result {}
	lappend result [exec diff tmp/scat.tsv tmp/check.tsv]
	set files [glob tmp/distrvars1-*]
	foreach file $files {
		set temp [lindex [split $file -] end]
		catch {exec grep -v $temp $file} r
		if {$r ne "child process exited abnormally"} {error "file $file error: $r"}
	}
	lappend result [llength $files]
	lappend result [file_read tmp/distrvars1-chr2]
	set result
} {{} 81 {chr2	2
chr2	2-2
}} 

test file_absolute {relative path} {
	set keep [pwd]
	cd /tmp
	set result [file_absolute abc]
	cd $keep
	set result
} /tmp/abc

test file_absolute {relative path with ..} {
	set keep [pwd]
	cd /usr/bin
	set result [file_absolute ../abc]
	cd $keep
	set result
} /usr/abc

test file_absolute {absolute with .. and .} {
	file_absolute /abc/def/../ghi/./j
} /abc/ghi/j

test file_absolute {empty} {
	file_absolute /abc//def
} /abc/def

test file_absolute {error} {
	file_absolute /abc/def/../../../ghi/./j
} {file_absolute error: cannot .. past root} error

cd $keepdir
file delete -force {*}[glob tmp/*]

testsummarize
