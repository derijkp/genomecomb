#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

package require pkgtools
if {![llength $argv]} {
	error "format is: findsymlinks.tcl directory"
}
set startdir [lindex $argv 0]
puts "updating $startdir/update_symlinks.sh"

proc update_symlinks {dir resultvar} {
	upvar $resultvar result
	foreach file [glob -nocomplain $dir/*] {
		if {[file tail $file] eq "_darcs"} continue
		if {![catch {file link $file} link]} {
			lappend result [list ln -sf $link $file]
		} elseif {[file isdir $file]} {
			update_symlinks $file result
		} else {
#			if {[file executable $file]} {
#				lappend result [list chmod ugo+x $file]
#			}
			set permissions [file attributes $file -permissions]
			lappend result [list chmod $permissions $file]
		}
	}
}

cd $startdir
update_symlinks . result
set o [open $startdir/bin/update_symlinks.sh w]
puts $o "#!/bin/sh"
puts $o [join $result \n]
close $o
exec chmod u+x $startdir/bin/update_symlinks.sh
