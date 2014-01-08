#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc multireg {compar_file file} {
	global cache comparposs1 mergeposs1 comparposs2 mergeposs2 dummy1 dummy2 restposs1 restposs2

	set name [file root [file tail [gzfile $file]]]
	catch {close $f1}; catch {close $f2}; catch {close $o}
	set f2 [gzopen $file]
	set poss2 [open_region $f2 h2]
	set num 0
	if {![file exists $compar_file]} {
		close $f2
		set h2base [list_sub $h2 $poss2]
		cg checksort $file
		cg select \
			-f "chromosome=\$[lindex $h2base 0] begin=\$[lindex $h2base 1] end=\$[lindex $h2base 2] $name=1" \
			$file $compar_file
		return
	}
	set f1 [open $compar_file]
	set poss1 [open_region $f1 h1]
	if {[inlist $h1 $name]} {
		error "$name already present in $compar_file"
	}
	set o [open $compar_file.temp w]
	puts $o [join $h1 \t]\t$name
	set dummy1 [list_fill [expr {[llength $h1]-3}] 0]
	close $o
	# putslog "multireg $compar_file $poss1 $dummy1 $file $poss2 >> $compar_file.temp"
	exec multireg $compar_file {*}$poss1 [join $dummy1 \t] $file {*}$poss2 >> $compar_file.temp 2>@stderr
	catch {file rename -force $compar_file $compar_file.old}
	file rename -force $compar_file.temp $compar_file	
}

proc cg_multireg {args} {
	if {([llength $args] < 1)} {
		errorformat multireg
		exit 1
	}
	foreach {compar_file} $args break
	if {[file exists $compar_file]} {
		set f [gzopen $compar_file]
		set header [tsv_open $f]
		close $f
	} else {
		set header {chromosome begin end}
	}
	set files [lrange $args 1 end]
	foreach file $files {
		set name [file root [file tail [gzfile $file]]]
		if {[inlist $header $name]} {
			putslog "*** Skipping $file: $name already in $compar_file ***"
			continue
		}
		putslog "Adding $file to $compar_file"
		multireg $compar_file $file
	}
}

if {[info exists argv0] && [file tail [info script]] eq [file tail $argv0]} {
	package require pkgtools
	set appdir [file dir [pkgtools::startdir]]
	lappend auto_path $appdir/lib
	append env(PATH) :[file dir [file dir $appdir]]/bin:$appdir/bin
	package require Extral
	set ::base [file tail [info script]]
	cg_multireg {*}$argv
}
