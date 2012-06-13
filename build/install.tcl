#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

package require pkgtools
cd [pkgtools::startdir]

# settings
# --------

set libfiles {pkgIndex.tcl init.tcl libext}
set shareddatafiles {}
set headers {}
set libbinaries [::pkgtools::findlib [file dir [pkgtools::startdir]] genomecomb]
puts "libbinaries: $libbinaries pkgtools::startdir:[pkgtools::startdir]"
set binaries {}

# standard
# --------
pkgtools::install $argv
