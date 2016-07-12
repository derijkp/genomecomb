#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc help_get {action} {
	if {[file exists $::appdir/lib/cg_$action.wiki]} {
		set help [file_read $::appdir/lib/cg_$action.wiki]
	} elseif {[file exists $::appdir/lib-exp/cg_$action.wiki]} {
		set help [file_read $::appdir/lib-exp/cg_$action.wiki]
	} elseif {[file exists $::appdir/plugins/cg_$action.wiki]} {
		set help [file_read $::appdir/plugins/cg_$action.wiki]
	} elseif {[file exists $::appdir/docs/$action.wiki]} {
		set help [file_read $::appdir/docs/$action.wiki]
	} else {
		puts stderr "Unknown help topic \"$action\", known topics are:"
		puts stderr "Docs:\n[help_docs]"
		puts stderr "Commands:\n[help_actions]"
		puts stderr "Use without arguments for overview"
		exit 1
	}
}

proc help_actions {} {
	global appdir
	set files [ssort -natural [list_concat [dirglob $appdir/lib/ cg_*.wiki] [dirglob $appdir/lib-exp/ cg_*.wiki]]]
	set list {}
	foreach file $files {
		set action [string range [file root [file tail $file]] 3 end]
		lappend list $action
	}
	return $list
}

proc help_docs {} {
	global appdir
	set files [ssort -natural [dirglob $appdir/docs/ *.wiki]]
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

proc help_formatw {foutput text width mode indent {format 1}} {
# putsvars text width indent format mode
#foreach var {text width indent format mode} {
#	puts $foutput [list set $var [get $var]]
#}
	if {$text eq ""} {return ""}
	if {$format} {
		set green "\033\[1;32m"
		set yellow "\033\[1;33m"
		set cyan "\033\[1;36m"
		set normal "\033\[0m"
	} else {
		set green ""
		set yellow ""
		set cyan ""
		set normal ""
	}
	set w [expr {$width - $indent}]
	set pre [string_fill " " $indent]
	set collect {}
	foreach line [split [string trim $text] \n] {
		# keep {{{literals}}} to put back in later
		set list [regexp -inline -all {\{\{\{.*?\}\}\}} $line]
		regsub -all {\{\{\{.*?\}\}\}} $line {{{{}}}} line
		# **bold**
		regsub -all {\*\*([^\n*]+)\*\*} $line "${yellow}\\1$normal" line
		# \ escaped characters
		regsub -all {\\([^A-Za-z0-9])} $line {\1} line
		# [[link]]
		regsub -all {\[\[([^|]*)\]\]} $line "${cyan}\\1$normal" line
		# [[link|text]]
		regsub -all {\[\[([^|]*)\|([^|]*)\]\]} $line "${cyan}\\2$normal" line
		# place extracted literals back
		foreach value $list {
			regsub {\{\{\{.*?\}\}\}} $line ${green}[string range $value 3 end-3]$normal line
		}
		if {$collect ne ""} {append collect " "}
		append collect [string trim $line]
		# format (separate lines and indent)
		while {[string length $collect] >= $w} {
			set pos [string last " " $collect $w]
			if {$pos == -1} {set pos [string last "\t" $collect $w]}
			if {$pos == -1} {set pos $w}
			append newtemp "$pre[string range $collect 0 [expr {$pos-1}]]\n"
			incr pos
			set collect [string range $collect $pos end]
			if {$mode eq "l"} {
				set w [expr {$width - $indent - 2}]
				set pre "$pre  "
				set mode {}
			}
		}
	}
	if {$collect ne ""} {append newtemp "$pre$collect\n"}
	puts -nonewline $foutput $newtemp
}

proc help {action {format 1}} {
	set help [help_get $action]
	if {$format} {
		set bold "\033\[1;1m"
		set underline "\033\[1;4m"
		set green "\033\[1;32m"
		set yellow "\033\[1;33m"
		set cyan "\033\[1;36m"
		set normal "\033\[0m"
		set foutput [open "| less -r" w]
	} else {
		set bold ""
		set underline ""
		set green ""
		set yellow ""
		set cyan ""
		set normal ""
		set foutput stdout
	}
	set width 80
	if {![catch {exec stty -a} result] && [regexp {columns (\d+)} $result temp cols]} {
		if {$cols > 10} {set width $cols}
	}

	set indent 0
	set mode {}
	set output {}
	set collect {}
	foreach line [split $help \n] {
# puts $foutput [list set line [get line]]
		set fchr [string index $line 0]
		if {$mode eq "b"} {
			if {$line eq "\}\}\}"} {
				puts -nonewline $foutput $normal
				set mode $model
			} else {
				puts -nonewline $foutput [string_fill " " $indent]${line}\n
			}
		} elseif {$fchr eq ""} {
			help_formatw $foutput $collect $width $mode $indent $format
			puts -nonewline $foutput \n
			set collect {}
			set mode {}
			set indent 0
		} elseif {$fchr eq "="} {
			help_formatw $foutput $collect $width $mode $indent $format
			set collect {}
			set mode {}
			set indent 0
			if {[string index $line 3] eq "="} {
				# == heading4 ==
				puts -nonewline $foutput [string trim $line " ="]\n
			} elseif {[string index $line 2] eq "="} {
				# == heading3 ==
				puts -nonewline $foutput ${cyan}[string trim $line " ="]$normal\n
			} elseif {[string index $line 1] eq "="} {
				# == heading2 ==
				puts -nonewline $foutput ${underline}[string trim $line " ="]$normal\n
			} else {
				# = heading1 =
				puts -nonewline $foutput ${underline}${green}[string trim $line " ="]$normal\n
			}
		} elseif {[regexp {^\; } $line]} {
			# ; terms : definitions
			help_formatw $foutput $collect $width $mode $indent $format
			set collect {}
			set indent 5
			if {[regexp {;? *([^\n]+?) *: +(.*)} $line tmp term def]} {
				set mode d
				regsub -all {\*\*([^\n*]+)\*\*} $term "${yellow}\\1$normal" term
				puts -nonewline $foutput "  ${yellow}$term$normal\n"
				append collect $def\n
			} else {
				append collect "  * $line\n"
			}
		} elseif {[regexp {^\* } $line]} {
			# ; terms : definitions
			help_formatw $foutput $collect $width $mode $indent $format
			set collect {}
			set indent 2
			set mode l
			append collect $line\n
		} elseif {$line eq "\{\{\{"} {
			# {{{literal}}}
			help_formatw $foutput $collect $width $mode $indent $format
			set collect {}
			set model $mode
			set mode b
			puts -nonewline $foutput ${green}
		} elseif {[string range $line end-1 end] eq "  "} {
			append collect $line\n
			help_formatw $foutput $collect $width $mode $indent $format
			set collect {}
		} else {
			append collect $line\n
		}
	}
	help_formatw $foutput $collect $width $mode $indent $format

	if {$format} {
		close $foutput
	}
}

proc cg_help {args} {
	global appdir
	set format 1
	cg_options help args {
		-format {set format $value}
	} 0 1
	set item [lindex $args 0]
	if {$item eq "distr"} {
		set files [glob -nocomplain $appdir/lib/cg_*.wiki]
	} elseif {$item ne ""} {
		help $item $format
		exit
	} else {
		set files [glob -nocomplain $appdir/lib/cg_*.wiki $appdir/lib-exp/cg_*.wiki]
	}
	puts "= Reference ="
	puts ""
	puts "== Format =="
	puts "cg action ...."
	puts ""
	puts "== Actions =="
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
