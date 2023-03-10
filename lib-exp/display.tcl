#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

lappend auto_path ~/dev/genomecomb/lib
package require ClassyTk
package require Extral
package require dbi_sqlite3

Widget subclass cgdisplay

cgdisplay method init {args} {
	super init
	private $object start end
	Classy::Canvas $object.c
	pack $object.c -fill both -expand yes
	set start 0
	set end 10000
	dbi_sqlite3 $object.db
	if {"$args" != ""} {$object configure {*}$args}
	Classy::todo $object redraw
	return $object
}

cgdisplay method open {file} {
	$object.db open $file
}

cgdisplay method redraw {} {
return
	private $object start end drawna
	set bins [lsort -integer [region2bins $start $end]]
	foreach bin $bins {
		if {[info exists drawna($bin)]} continue
		set hits [$object.db exec {
			select strand1,start1,end1,strand2,start2,end2,type
			from hits where bin = ?
		} $bin]
		set drawna($bin) 1
	}
}

if 0 {
cd /complgen/sv/
set file sv79-20s.sqlite
set object .d
cgdisplay .d
pack .d -fill both -expand yes
.d open $file
set db $object.db

.d redraw


}
