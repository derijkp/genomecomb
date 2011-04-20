proc cg_help {} {
global appdir
puts {
== Format ==
cg action ....

== Actions == }
set files [dirglob $appdir/lib/ cg_*.help]
unset -nocomplain a
foreach file $files {
	set action [string range [file root [file tail $file]] 3 end]
	set h [helpparts $action]
	set category {}
	catch {dict get $h Category} category
	set descr {}
	set item " * $action"
	if {![catch {dict get $h Summary} summary]} {
		append item ": $summary"
	}
	lappend a($category) $item
}
set categories [array names a]
set pre {Conversion Annotation Query Regions Structural}
set categories [list_concat [list_common $pre $categories] [list_lremove $categories $pre]]
foreach category $categories {
	puts " === $category ==="
	puts [join $a($category) \n]\n
}

puts " === Other ==="
puts { * select, graph, multicompar, regsubtract, regjoin, regcommon, makeregions, makeprimers, ...
}
}

proc help {action} {
	set help [file_read $::appdir/lib/cg_$action.help]
	puts $help
}

proc help_actions {} {
	global appdir
	set files [lsort -dict [dirglob $appdir/lib/ cg_*.help]]
	set list {}
	foreach file $files {
		set action [string range [file root [file tail $file]] 3 end]
		lappend list $action
	}
	return $list
}

proc errorformat {action} {
	set help [helpparts $action]
	puts stderr "\nERROR: Wrong number of arguments, correct format is:"
	puts stderr [dict get $help Format]
	puts stderr "\nFor more help, use:\ncg $action -h\n"
}

proc helpparts {action} {
	set help [file_read $::appdir/lib/cg_$action.help]
	regsub -all {[ \n\t]*== *([^=]+?) *==[ \n\t]*} $help {@@@@\1@@@@} help
	set result [lrange [string_split $help @@@@] 1 end]
	return $result
}
