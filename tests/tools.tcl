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

lappend auto_path $appdir/lib
