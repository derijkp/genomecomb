#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

# $Format: "set version $ProjectMajorVersion$.$ProjectMinorVersion$.$ProjectPatchLevel$"$
set version 0.98.3

set keeppath $::env(PATH)
set testlist [list /home/peter/build/genomecomb-$version-Linux-i686 /home/peter/build/genomecomb-$version-Linux-x86_64]
set cgversion [lindex $testlist 0]

foreach cgversion $testlist {
	puts "\n\n==================== Testing $cgversion ===================="
	set ::env(PATH) $cgversion:$::keeppath
	# source annot.tcl
	source all.tcl
}
