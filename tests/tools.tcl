package require Extral
catch {tk appname test}

package require pkgtools
namespace import pkgtools::*
package require Extral

# pkgtools::testleak 100

set keeppath $::env(PATH)
set script [info script] ; if {$script eq ""} {set script ./t}
set testdir [file dir [file normalize $script]]
set appdir [file dir [file dir [file normalize $script]]]
append ::env(PATH) :$appdir/bin
# putsvars ::env(PATH)
set env(SCRATCHDIR) [file dir [tempdir]]

proc test_cleantmp {} {
	foreach file [glob -nocomplain tmp/*] {
		catch {file attributes $file -permissions ugo+xw}
		catch {file delete -force $file}
	}
	cg indexclean
}

proc write_tab {file data {comment {}}} {
	set data [split [string trim $data] \n]
	set f [open $file w]
	if {$comment ne ""} {puts $f $comment}
	foreach line $data {
		puts $f [join $line \t]
	}
	close $f
}

proc diff_tab {file data {comment {}}} {
	set tempfile [tempfile]
	set data [split [string trim $data] \n]
	set f [open $tempfile w]
	if {$comment ne ""} {puts $f $comment}
	foreach line $data {
		puts $f [join $line \t]
	}
	close $f
}

lappend auto_path $appdir/lib $appdir/lib-exp $appdir/libext

file mkdir tmp
test_cleantmp
set dbopt {}
set dboptt {}
