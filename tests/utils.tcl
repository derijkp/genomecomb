#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

# cg select -q 'oneof($locus,1224,1226,1598520,2816907,2816911,2818387,3038185,3042129,5019241,5054443,5211957,9605719,10679396,19364044,21233191,22632054,21542413)' cgNA19240/fannotvar-cgNA19240.tsv > testvars.tsv
# exit

test mklink {basic} {
	test_cleantmp
	file_write tmp/test a
	cg mklink tmp/test tmp/test2
	file_read tmp/test2
} a

test mklink {relative} {
	test_cleantmp
	file_write tmp/test a
	mkdir tmp/dir
	cg mklink tmp/test tmp/dir/test2
	file link tmp/dir/test2
} ../test

test mklink {absolute} {
	test_cleantmp
	file_write tmp/test a
	mkdir tmp/dir
	cg mklink -absolute 1 tmp/test tmp/dir/test2
	file link tmp/dir/test2
} /*/tmp/test match

test mklink {do not overwrite real file} {
	test_cleantmp
	file_write tmp/test a
	file_write tmp/test2 b
	cg mklink tmp/test tmp/test2
	file_read tmp/test2
} {cg mklink error: destination exists and is not a link: *} error match

test mklink {dest is a link to non-existing file} {
	test_cleantmp
	file_write tmp/test a
	exec ln -s unexisting tmp/test2
	cg mklink tmp/test tmp/test2
	file_read tmp/test2
} a

testsummarize
