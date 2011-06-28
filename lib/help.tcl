#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg_help {{item {}}} {
global appdir
if {$item eq "distr"} {
	set files [glob -nocomplain $appdir/lib/cg_*.wiki]
} elseif {$item ne ""} {
	help $item
	exit
} else {
	set files [glob -nocomplain $appdir/lib/cg_*.wiki $appdir/lib-exp/cg_*.wiki]
}
puts {
= Reference =

== Format ==
cg action ....

== Actions == }
unset -nocomplain a
foreach file $files {
	set action [string range [file root [file tail $file]] 3 end]
	set h [helpparts $action]
	set category {}
	catch {dict get $h Category} category
	set category [string trim $category]
	set descr {}
	set item " * $action"
	if {![catch {dict get $h Summary} summary]} {
		append item ": $summary"
	}
	lappend a($category) $item
}
set categories [array names a]
set pre {Conversion Annotation Compare Query Regions Structural}
set categories [list_concat [list_common $pre $categories] [list_lremove $categories $pre]]
foreach category $categories {
	puts " === $category ==="
	puts [join $a($category) \n]\n
}

puts " === Other ==="
puts { * select, graph, multicompar, regsubtract, regjoin, regcommon, makeregions, makeprimers, ...
}
}

proc help_get {action} {
	if {[file exists $::appdir/lib/cg_$action.wiki]} {
		set help [file_read $::appdir/lib/cg_$action.wiki]
	} elseif {[file exists $::appdir/lib-exp/cg_$action.wiki]} {
		set help [file_read $::appdir/lib-exp/cg_$action.wiki]
	} elseif {[file exists $::appdir/docs/$action.wiki]} {
		set help [file_read $::appdir/docs/$action.wiki]
	} else {
		puts stderr "Unknown help subject \"$action\""
		puts stderr "Known subjects are"
		puts stderr "Docs: [help_docs]"
		puts stderr "Commands: [help_actions]"
		puts stderr "Use without arguments for overview"
		exit 1
	}
}

proc help {action} {
	set help [help_get $action]
	set cyan "\033\[1;36m"
	set normal "\033\[0m"
#	regsub -all {(^|\n)\; *([^:\n]+):} $help "\n${cyan}\\2$normal:" help
	regsub -all {\*\*([^\n*]+)\*\*} $help "${cyan}\\1$normal" help
	puts $help
}

proc help_actions {} {
	global appdir
	set files [lsort -dict [list_concat [dirglob $appdir/lib/ cg_*.wiki] [dirglob $appdir/lib-exp/ cg_*.wiki]]]
	set list {}
	foreach file $files {
		set action [string range [file root [file tail $file]] 3 end]
		lappend list $action
	}
	return $list
}

proc help_docs {} {
	global appdir
	set files [lsort -dict [dirglob $appdir/docs/ *.wiki]]
	set list {}
	foreach file $files {
		set action [string range [file root [file tail $file]] 0 end]
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
	set help [help_get $action]
	regsub -all {[ \n\t]*== *([^=]+?) *==[ \n\t]*} $help {@@@@\1@@@@} help
	set result [lrange [string_split $help @@@@] 1 end]
	return $result
}
