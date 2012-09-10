#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

set keepdir [pwd]
catch {file delete -force {*}[glob tmp/*]}

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

cd $keepdir
file delete -force {*}[glob tmp/*]

testsummarize
