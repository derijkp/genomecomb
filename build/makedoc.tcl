#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

package require pkgtools
cd [pkgtools::startdir]

# settings
# --------

set gcdir [file dir [pkgtools::startdir]]

# set gcdir [pwd]

package require Extral

set destdir $gcdir/docs

# make html procs
# ---------------

proc convert.html {file destfile} {
	if {![checkdeps $destfile $file] && ![newtemplate]} return
	puts "convert $file $destfile"
	set c [file_read $file]
	set title {}
	regexp -nocase {<title>([^<>]*)</title>} $c temp title
	regsub -all -nocase {</?html>} $c {} c
	regsub -all -nocase {<header>.*</header>} $c {} c
	writepage $file $destfile $c $title
}

proc convert.wiki {file destfile} {
	global tempfile
	if {![checkdeps $destfile $file] && ![newtemplate]} return
	puts "convert $file $destfile"
	set c [file_read $file]
	set title [file tail [file root $file]]
	puts subst:[regsub -all -nocase {\[\[cg_([^]|]+)\]\]} $c {[[cg_\1|cg \1]]} c]
	file_write $tempfile $c	
	set c [exec nme --autourllink --body --easylink {$.html} --xref < $tempfile]
	regsub -all {<a[ \n]+href="\.\.} $c {<a href="../} c
	writepage $file $destfile $c $title
}

proc checkdeps {destfile args} {
	if {$::force} {return 1}
	if {![file exists $destfile]} {return 1}
	set dtime [file mtime $destfile]
	foreach src $args {
		set error [catch {file mtime $src} stime]
		if {$error || ($dtime < $stime)} {
			return 1
		}
	}
	return 0
}

proc init {src dest} {
	global header footer srcroot destroot newtemplate force tempfile
	if {$force || [checkdeps $dest $src/docs-src/templates/header.html $src/docs-src/templates/footer.html]} {
		set newtemplate 1
	} else {
		set newtemplate 0
	}
	set srcroot $src
	set destroot $dest
	set header [file_read $src/docs-src/templates/header.html]
	set footer [file_read $src/docs-src/templates/footer.html]
	set tempfile $src/build/tempfile
	puts "tempfile=$tempfile"
}

proc newtemplate {} {
	return $::newtemplate
}

proc convertlinks {file data} {
	regsub -all {[^/]+} $file .. pre
	set pre [string range $pre 0 end-2]
	regsub -all -nocase {href *= *"/} $data "href=\"$pre" data
	regsub -all -nocase {src *= *"/} $data "src=\"$pre" data
	return $data
}

proc getheader {file {title {}}} {
	global header
	set nheader $header
	if {$title ne ""} {
		regsub -nocase {<title>([^<>]*)</title>} $header "<title>$title</title>" nheader
	}
	return [convertlinks $file $nheader]
}

proc getfooter {file} {
	global footer
	return [convertlinks $file $footer]
}

proc writepage {file destfile c {title {}}} {
	set o [open $destfile w]
	puts $o [getheader $file $title]
	puts $o $c
	puts $o [getfooter $file]
	close $o
}

proc tohtml {file dest} {
	set ext [file extension $file]
	if {[info commands convert$ext] eq ""} {
		set destfile $dest/[file tail $file]
		puts "copy $file $destfile"
		file copy -force $file $destfile
	}
	if {[file isdir $dest]} {
		set destfile $dest/[file root [file tail $file]].html
	} else {
		set destfile $dest
	}
	convert$ext $file $destfile
}

# extract help for reference section
# ----------------------------------

append env(PATH) :$gcdir/docs-src

file delete -force $destdir.old
if {[file exists $destdir]} {file rename $destdir $destdir.old}
file mkdir $destdir

# clean
set force 1
init $gcdir $destdir

set files [glob $gcdir/lib/*.wiki]
foreach file $files {
	puts $file
	tohtml $file $destdir
}

# make the reference overview for docs without using cg
# by loading and using several parts of it
set appdir $gcdir
file delete -force $destdir/reference.wiki
source $gcdir/lib/cg_help.tcl
source $gcdir/lib/tools.tcl
# standard version won't work, so use this hack
proc version {args} {
	set c [file_read $::gcdir/build/version.tcl]
	regexp {set version ([^ \n]+)} $c temp version
	return $version
}
set help [helptext_overview]

regsub -all {\* +([^:]+):} $help {* [[cg_\1|\1]]:} help
file_write $destdir/reference.wiki $help
tohtml $destdir/reference.wiki $destdir

# extract help for howto section
file mkdir $destdir
set files [glob $gcdir/help/*.wiki]
foreach file $files {
	puts $file
	tohtml $file $destdir
}

file delete -force $destdir/css
file copy -force $gcdir/docs-src/css $destdir
file delete -force $destdir/img
file copy -force $gcdir/docs-src/img $destdir
set files [glob $gcdir/docs-src/*.html]
foreach file $files {
	puts $file
	tohtml $file $destdir
}
