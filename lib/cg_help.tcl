#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc help_get {action} {
	if {$action eq ""} {
		helptext_overview
	} elseif {[file exists $::appdir/lib/cg_$action.wiki]} {
		set help [file_read $::appdir/lib/cg_$action.wiki]
	} elseif {[file exists $::appdir/lib-exp/cg_$action.wiki]} {
		set help [file_read $::appdir/lib-exp/cg_$action.wiki]
	} elseif {[file exists $::appdir/plugins/cg_$action.wiki]} {
		set help [file_read $::appdir/plugins/cg_$action.wiki]
	} elseif {[file exists $::appdir/docs/$action.wiki]} {
		set help [file_read $::appdir/docs/$action.wiki]
	} elseif {[auto_load helptext_$action]} {
		set help [helptext_$action]
	} elseif {[auto_load ${action}_job]} {
		catch {set help [${action}_job -h]} e
	} elseif {[auto_load cg_$action]} {
		catch {set help [cg_$action -h]} e
	}
	if {![info exists help]} {
		set msg "Unknown help topic \"$action\", known topics are:\n\n"
		append msg "Docs:\n[help_docs]\n\n"
		append msg "Commands:\n[help_actions]\n\n"
		append msg "Use help without arguments for overview"
		error $msg
	} else {
		return $help
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

proc errorformat_calc {action {options {}} {minargs {}} {maxargs {}} {parameters {}}} {
	set out "cg $action"
	if {$options ne ""} {append out " ?options?"}
	set pos 0
	while {$pos < $minargs} {
		set p [lindex $parameters $pos]
		if {$p eq ""} {set p arg}
		append out " $p"
		incr pos
	}
	while {$maxargs eq "..." || $pos < $maxargs} {
		if {$pos >= [llength $parameters]} {
			append out " ..."
			break
		}
		set p [lindex $parameters $pos]
		if {$p eq ""} {set p arg}
		append out " ?$p?"
		incr pos
	}
	if {$options ne ""} {
		append out "\n  with options: $options"
	}
	return $out
}

proc errorformat {action {options {}} {minargs {}} {maxargs {}} {parameters {}}} {
# putsvars action options minargs maxargs parameters
	if {![catch {
		set help [helpparts $action]
	}] && [dict exists $help Format]} {
		set msg "\nERROR: Wrong number of arguments, correct format is:"
		append msg \n[dict get $help Format]
		if {[dict exists $help Options]} {
			set options [dict get $help Options]
			list_unmerge [regexp -all -inline {; *(-[^ ]+)} $options] 1 temp
			if {[llength $temp]} {
				append msg "\n  with options: [join $temp ,]"
			}
		}
		append msg "\n\nFor more help, use:\ncg $action -h\n"
		error $msg
	} else {
		set msg "\nERROR: Wrong number of arguments, correct format is:"
		append msg \n[errorformat_calc $action $options $minargs $maxargs $parameters]
		error $msg
	}
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
	if {$format eq "md"} {
		puts $help
		return
	} elseif {$format} {
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
				regsub {\[\[[^\|]+\|(.*)\]\]} $term {\1} term
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

proc helptext_overview {} {
	global appdir
	set help "= Genomecomb [version genomecomb] Reference =\n\n"
	append help [string trim [deindent {
		== Format ==
		cg subcommand ?options? ....
		
		== Description ==
		This help page gives a (reference style) overview of all genomecomb functions. For an 
		introductory text to genomecomb and its formats, use
		{{{
		cg help intro
		}}}
		All genomecombs functions are called using the cg command and a subcommand. 
		The available subcommands are listed on this page (in categories) with a short description.
		To get further info on how to use the subcommands and their parameters, use
		{{{
		cg help subcommand
		}}}
		or 
		{{{
		cg subcommand -h
		}}}
		
		== Options ==
		The following options are generic and available for all subcommands. They must however
		always preceed the subcommand specific options.
		; -v number (--verbose): Setting this to 1 or 2 (instead of the default 0) makes some subcommands chattier about their progress.
		At the given number is 1, logging messages are shown (warnings, start of subtask, etc.)
		If the number >= 2, progress counters are also shown (for commands that support them)
		; --stack 0/1: When the program returns an error, by default only the error message is shown,
		which is normally ok to show errors in input format, etc.
		if --stack is set to 1, a full stack trace is shown on error (which may be useful to solve 
		errors caused by bugs in the program)
		
		== Available subcommands ==
	}]] \n
	unset -nocomplain a
	set files [glob -nocomplain $appdir/lib/cg_*.wiki $appdir/lib-exp/cg_*.wiki]
	foreach file $files {
		set action [string range [file root [file tail $file]] 3 end]
		set h [helpparts $action]
		set category {}
		if {[catch {dict get $h Category} category]} continue
		set category [string trim $category]
		set descr {}
		set item "; \[\[cg_$action|$action\]\]"
		if {![catch {dict get $h Summary} summary]} {
			append item ": $summary"
		}
		lappend a($category) [list $action $item]
	}
	unset -nocomplain a(Depricated)
	set categories [array names a]
	set pre {Process Query Analysis Regions Annotation Validation tsv Conversion {Format Conversion} Compare Structural Report Info}
	array set preorder {
		Process {process_project}
		Query {select viz multiselect groupby}
		Regions {multireg regselect}
		Conversion {sam_clipamplicons liftover liftregion liftsample liftchain2tsv correctvariants bamreorder9}
	}
	set categories [list_concat [list_common $pre $categories] [list_lremove $categories $pre]]
	foreach category $categories {
		append help "\n=== $category ===\n"
		set list $a($category)
		foreach item [get preorder($category) ""] {
			set pos [lsearch [list_subindex $list 0] $item]
			if {$pos != -1} {
				set line [list_pop list $pos]
				append help [lindex $line 1]\n
			}
		}
		set list [lsort -index 0 $list]
		append help [join [list_subindex $list 1] \n]\n
	}
	return $help
}

proc cg_help {args} {
	global appdir
	set format 1
	set item overview
	cg_options help args {
		-format {set format $value}
	} item 0 1
	if {$item ne ""} {
		help $item $format
	} else {
		help overview
	}
}
