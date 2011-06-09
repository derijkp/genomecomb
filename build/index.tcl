#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

package require Class
package require pkgtools
set dir [file dir [pkgtools::startdir]]
puts "rebuilding $dir/lib/tclIndex"
Class::auto_mkindex $dir/lib
puts "rebuilding $dir/lib-exp/tclIndex"
Class::auto_mkindex $dir/lib-exp
