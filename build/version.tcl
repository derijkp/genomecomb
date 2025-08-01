#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

set version 0.114.0
set extversion 0.114.0

# standard
# --------
package require pkgtools
if {$argv ne ""} {
	set version $argv
}
puts "version: setting version to $version"
pkgtools::version $version

if {[info exists ::srcdir]} {
	set srcdir $::srcdir
} else {
	set srcdir [file dir [pkgtools::startdir]]
}
proc file_change {file args} {
	puts "updating version in $file"
	set f [open $file]
	set c [read $f]
	close $f
	# set c [string map $args $c]
	foreach {pattern subst} $args {
		regsub -all $pattern $c $subst c
	}
	set o [open $file.temp w]
	puts -nonewline $o $c
	close $o
	file rename -force $file.temp $file
}
file_change $srcdir/help/illumina_rna_workflow_description.txt {genomecomb [0-9.]+,} "genomecomb $version,"
file_change $srcdir/README.md {genomecomb-[0-9.]+-} "genomecomb-$version-" {download/[0-9.]+/} "download/$version/"

puts "version in configure.in is not updated as the extension does not change often, and can stay with older versions"
puts "If source of the extension has changed, change version manually (or using maketea) and run autoconf"

# If we would want to change the version of the extension
proc change_configure.in {version} {
	if {[info exists ::srcdir]} {
		set srcdir $::srcdir
	} else {
		set srcdir [file dir [pkgtools::startdir]]
	}
	puts "version: rewriting $srcdir/configure.in to $version"
	set f [open $srcdir/configure.in]
	set c [read $f]
	close $f
	if {![regsub \
		{AC_INIT\(\[genomecomb\], \[([^]]+)\]\)} \
		$c \
		"AC_INIT\(\[genomecomb\], \[$version\]\)" \
		c
	]} {
		error "Could not replace version with $version in $srcdir/configure.in"
	}
	
	set o [open $srcdir/configure.in w]
	puts -nonewline $o $c
	close $o
}
