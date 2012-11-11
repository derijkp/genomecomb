#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

set keepdir [pwd]
test_cleantmp

file mkdir tmp
file mkdir tmp/test
file_write tmp/test/test1 1
file_write tmp/test/test2 2
file mkdir tmp/out

proc dirinfo {dir} {
	set result {}
	foreach file [lsort -dict [glob $dir/*]] {
		if {[catch {file link $file} link]} {set link -}
		lappend result [list $file $link [file_read $file]]
	}
	return $result
}

test cplinked {relative} {
	file mkdir tmp/out
	exec cg cplinked tmp/test tmp/out
	set result [dirinfo tmp/out/test]
	file delete -force tmp/out/test
	set result
} {{tmp/out/test/test1 ../../test/test1 1} {tmp/out/test/test2 ../../test/test2 2}} 

test cplinked {relative, new dir} {
	file delete -force tmp/out
	file mkdir tmp/out
	exec cg cplinked tmp/test tmp/out/test2
	set result [dirinfo tmp/out/test2]
	file delete -force tmp/out/test2
	set result
} {{tmp/out/test2/test1 ../../test/test1 1} {tmp/out/test2/test2 ../../test/test2 2}} 

test cplinked {relative, new dir 2} {
	file delete -force tmp/out
	exec cg cplinked tmp/test tmp/out/test
	set result [dirinfo tmp/out/test]
	file delete -force tmp/out/test
	set result
} {{tmp/out/test/test1 ../../test/test1 1} {tmp/out/test/test2 ../../test/test2 2}} 

test cplinked {relative extra level} {
	file delete -force tmp/out
	file mkdir tmp/out/test2
	exec cg cplinked tmp/test tmp/out/test2/testl
	set result [dirinfo tmp/out/test2/testl]
	file delete -force tmp/out
	set result
} {{tmp/out/test2/testl/test1 ../../../test/test1 1} {tmp/out/test2/testl/test2 ../../../test/test2 2}} 

test cplinked {absolute} {
	file delete -force /tmp/test2
	exec cg cplinked tmp/test /tmp/test2
	set result [dirinfo /tmp/test2]
	file delete -force /tmp/test2
	set result
} [list [list /tmp/test2/test1 $keepdir/tmp/test/test1 1] [list /tmp/test2/test2 $keepdir/tmp/test/test2 2]]

test cplinked {backup existing} {
	file delete -force tmp/out
	file mkdir tmp/out/test
	file_write tmp/out/test/test1 pre
	exec cg cplinked tmp/test tmp/out 2> /dev/null
	set result [dirinfo tmp/out/test]
	file delete -force tmp/out/test
	set result
} {{tmp/out/test/test1 ../../test/test1 1} {tmp/out/test/test1.old1 - pre} {tmp/out/test/test2 ../../test/test2 2}} 

test cplinked {backup existing 2} {
	file delete -force tmp/out
	file mkdir tmp/out/test
	file_write tmp/out/test/test1 pre
	file_write tmp/out/test/test1.old1 old
	exec cg cplinked tmp/test tmp/out 2> /dev/null
	set result [dirinfo tmp/out/test]
	file delete -force tmp/out/test
	set result
} {{tmp/out/test/test1 ../../test/test1 1} {tmp/out/test/test1.old1 - old} {tmp/out/test/test1.old2 - pre} {tmp/out/test/test2 ../../test/test2 2}} 

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
	list [glob tmp/*] [file_read tmp/distrvars1-chr1]
} {{tmp/distrvars1-chromosome tmp/distrvars1-chr1 tmp/distrvars1-chr2} {chr1	4000	4001	snp	G	A	A	G	1	v	A	G	0	v	4
chr1	4001	4002	snp	A	G,C	G	G	1	v	G	G	0	v	1;2,3;4
chr1	4099	5000	snp	C	T	T	T	47	v	T	T	35	v	1,2
chr1	5000	5010	del	AGCGTGGCAA		AGCGTGGCAA		32	v	-	-	41	u	1;2
chr1	5020	5021	snp	G	C	G	C	54	v	G	G	52	r	3
}} 

cd $keepdir
file delete -force {*}[glob tmp/*]

testsummarize
