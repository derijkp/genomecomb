#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

set version 0.103.0
set extversion 0.103.0

# standard
# --------
package require pkgtools
if {$argv ne ""} {
	set version $argv
}
puts "version: setting version to $version"
pkgtools::version $version
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
