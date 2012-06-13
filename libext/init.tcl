# Initialisation of the genomecomb package
#
# Copyright (c) 2012 Peter De Rijk (VIB/UA)
#
# See the file "README.txt" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.
#
# =============================================================

namespace eval genomecomb {}

# $Format: "set ::genomecomb::version $ProjectMajorVersion$.$ProjectMinorVersion$"$
set ::genomecomb::version 0.8
# $Format: "set ::genomecomb::patchlevel $ProjectPatchLevel$"$
set ::genomecomb::patchlevel 6

package provide genomecomb $::genomecomb::version

package require pkgtools
pkgtools::init $genomecomb::dir genomecomb tsv_selectc [file join $genomecomb::dir]

#
# The lib dir contains the Tcl code defining the public genomecomb 
# functions. The lib dir is added to the auto_path so that
# these functions will be loaded on demand.
#

lappend auto_path [file join ${genomecomb::dir} libext]

if {![info exists appdir]} {
	set appdir ${genomecomb::dir}
	lappend auto_path $appdir/lib $appdir/lib-exp
	set env(PATH) $appdir/bin:$appdir/extern:[file dir [file dir $appdir]]/bin:$env(PATH)
}

