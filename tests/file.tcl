#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test shadow {cg shadow_mkdir} {
	test_cleantmp
	unset -nocomplain ::env(SHADOWDIR)
	cg shadow_mkdir -shadowdir tmp/shadow tmp/dir
	set shadowdir [glob tmp/shadow/*]
	bsort [list [file link tmp/dir] [file link $shadowdir/shadow_source]]
} {*/tmp/dir */tmp/shadow/shadow.*-dir} match

test shadow {cg shadow_mkdir SHADOWDIR} {
	test_cleantmp
	set ::env(SHADOWDIR) tmp/shadow
	cg shadow_mkdir tmp/dir
	set shadowdir [glob tmp/shadow/*]
	bsort [list [file link tmp/dir] [file link $shadowdir/shadow_source]]
} {*/tmp/dir */tmp/shadow/shadow.*-dir} match

test shadow {cg shadow_mkdir no shadowdir} {
	test_cleantmp
	unset -nocomplain ::env(SHADOWDIR)
	cg shadow_mkdir tmp/dir
	list [file exists tmp/dir] [catch {file link tmp/dir}] [glob -nocomplain tmp/shadow/*]
} {1 1 {}}

test shadow {cg shadow_delete} {
	test_cleantmp
	cg shadow_mkdir -shadowdir tmp/shadow tmp/dir
	set shadowdir [glob tmp/shadow/*]
	set result [bsort [list [file link tmp/dir] [file link $shadowdir/shadow_source]]]
	cg shadow_delete tmp/dir
	lappend result [file exists $shadowdir] [file exists tmp/dir]
} {*/tmp/dir */tmp/shadow/shadow.*-dir 0 0} match

test shadow {cg rm with shadows} {
	test_cleantmp
	cg shadow_mkdir -shadowdir tmp/shadow tmp/dir
	set shadowdir [glob tmp/shadow/*]
	set result [bsort [list [file link tmp/dir] [file link $shadowdir/shadow_source]]]
	cg rm tmp/dir
	lappend result [file exists $shadowdir] [file exists tmp/dir]
} {*/tmp/dir */tmp/shadow/shadow.*-dir 0 0} match

test shadow {cg shadow_clean} {
	test_cleantmp
	set shadowdir $::testdir/tmp/shadow
	cg shadow_mkdir -shadowdir $shadowdir tmp/dir
	cg shadow_mkdir -shadowdir $shadowdir tmp/dir2
	set result {}
	set shadowdir1 [file link tmp/dir]
	set shadowdir2 [file link tmp/dir2]
	lappend result $shadowdir1 $shadowdir2
	file delete tmp/dir
	lappend result [file exists $shadowdir1] [file exists $shadowdir2]
	unset  -nocomplain set ::env(SHADOWDIR)
	cg shadow_clean $shadowdir
	lappend result [file exists $shadowdir1] [file exists $shadowdir2]
} {*/tmp/shadow/shadow.*-dir*/tmp/shadow/shadow.*-dir2 1 1 0 1} match

test shadow {cg shadow_clean env var} {
	test_cleantmp
	set shadowdir $::testdir/tmp/shadow
	cg shadow_mkdir -shadowdir $shadowdir tmp/dir
	cg shadow_mkdir -shadowdir $shadowdir tmp/dir2
	set result {}
	set shadowdir1 [file link tmp/dir]
	set shadowdir2 [file link tmp/dir2]
	lappend result $shadowdir1 $shadowdir2
	file delete tmp/dir
	lappend result [file exists $shadowdir1] [file exists $shadowdir2]
	set ::env(SHADOWDIR) $shadowdir
	cg shadow_clean
	lappend result [file exists $shadowdir1] [file exists $shadowdir2]
} {*/tmp/shadow/shadow.*-dir*/tmp/shadow/shadow.*-dir2 1 1 0 1} match

test gzfiles {gzfiles basic} {
	test_cleantmp
	file_write tmp/test1 test1
	cg zst -stack 1 tmp/test1
	file_write tmp/test1 test1
	cg lz4 tmp/test1
	file_write tmp/test2 test2
	gzfiles tmp/test1
} {tmp/test1.zst}

test gzfiles {gzfiles complex} {
	test_cleantmp
	file_write tmp/s test
	cg bz2 tmp/s
	file_write tmp/s test
	cg zst tmp/s
	file_write tmp/s test
	cg bgz tmp/s
	file_write tmp/test1 test1
	cg rz tmp/test1
	file_write tmp/test1 test1
	cg lz4 tmp/test1
	file_write tmp/test2 test2
	file_write tmp/try_a2 test
	cg zst tmp/try_a2
	file_write tmp/try_a2 test
	cg lz4 tmp/try_a2
	file_write tmp/try_a1 test
	gzfiles tmp/try* tmp/test* tmp/s
} {tmp/try_a1 tmp/try_a2.zst tmp/test1.lz4 tmp/test2 tmp/s.zst}

test jobgzfiles {jobgzfiles complex - not testing "job" part} {
	test_cleantmp
	file_write tmp/s test
	cg bz2 tmp/s
	file_write tmp/s test
	cg zst tmp/s
	file_write tmp/s test
	cg bgz tmp/s
	file_write tmp/test1 test1
	cg rz tmp/test1
	file_write tmp/test1 test1
	cg lz4 tmp/test1
	file_write tmp/test2 test2
	file_write tmp/try_a2 test
	cg zst tmp/try_a2
	file_write tmp/try_a2 test
	cg lz4 tmp/try_a2
	file_write tmp/try_a1 test
	jobgzfiles tmp/try* tmp/test* tmp/s
} {tmp/try_a1 tmp/try_a2.zst tmp/test1.lz4 tmp/test2 tmp/s.zst}

testsummarize
