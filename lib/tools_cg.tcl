#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg {cmd args} {
	# puts "cg $cmd $args"
	catch_exec cg $cmd {*}$args
}

# default command line argument parsing for subcommands in genomecomb
# cmd: name of the subcommand it is being used in -> used for error messages and finding the help
# argsVar: variable name with arguments (usually args); This command changes the contrent of args: options that are recognized are removed from args
# def: block with options in the form -option optioncommand .... The value of the option is availble in the variable value
# parameters: list parameters/arguments to subcommand. values will be assigned to the variables in this list. If more args are available, they will stay in args var
# minargs: minimum number of arguments that must be given
# maxargs: maximum number of arguments that can be given
# summary: short help, used if the default help system cannot find the associated wiki format help file (which is search based on cmd)
#
# supports subcommand specific options (given as e.g. -sniffles-n) in generic interfaces like cg process_project by
# if the variable ::specialopt(-$cmd-$option) exists for a supported $option, it will be set to the given value
# 
proc cg_options {cmd argsVar def {parameters {}} {minargs {}} {maxargs ...} {summary {}} {optsVar {}}} {
# putsvars cmd argsVar def parameters minargs maxargs summary
	# we allways need options in the subst command anyway
	set optionlist [list_unmerge $def]
	set options [join [list_remove $optionlist default] ,]
	if {[inlist $optionlist default]} {append options ,...}
	set len [llength $parameters]
	if {$minargs eq ""} {set minargs $len}
	if {$maxargs eq "" || $maxargs eq "-1"} {set maxargs $len}
	upvar $argsVar args
	if {$args eq "-h"} {
		set format [errorformat_calc $cmd $options $minargs $maxargs $parameters]
		set help "= $cmd =\n\n== Format ==\n$format\n\n== Summary ==\n[string trim $summary]\n\n"
		return -code return $help
	}
	if {$maxargs eq "..."} {set test ""} else {set test " || (\$len > $maxargs)"}
	if {[lindex $def end-1] eq "default"} {
		set default [subst {
			if {\[string index \$key 0\] ne "-"} break
			[lindex $def end]
		}]
		set def [lrange $def 0 end-2]
	} else {
		set default [subst {
			if {\[string index \$key 0\] eq "-"} {
				error "error calling cg $cmd: unknown option \\"\$key\\", must be one of: $options"
			}
			break
		}]
	}
	if {$optsVar eq ""} {
		set optsdefault {}
	} else {
		set optsdefault "lappend $optsVar \$key \$value"
	}
	set fullcmd [subst {
		set pos 0
		while 1 {
			set key \[lindex \$$argsVar \$pos\]
			incr pos
			set value \[lindex \$$argsVar \$pos\]
			incr pos
			if {\$key eq "-"} break
			if {\[string range \$key 0 1\] eq "-."} break
			regsub {^--} \$key - key
			switch -glob -- \$key {
				$def
				- {incr pos ; break}
				default \{$default\}
			}
		}
		incr pos -2
		set $argsVar \[lrange \$$argsVar \$pos end\]
		set len \[llength \$$argsVar\]
		if {(\$len < $minargs) $test} {
			[list errorformat $cmd $options $minargs $maxargs $parameters]
		}
		foreach {key value} [list [specialopts -$cmd]] {
			regsub {^--} \$key - key
			switch -glob -- \$key {
				$def
				default \{$optsdefault\}
			}
		}
	}]
	if {[llength $parameters]} {
		append fullcmd "set ::_temp_max \[expr {\[llength \$$argsVar\] - 1}\]\n"
		append fullcmd "if {\$::_temp_max >= 0} {foreach \[lrange [list $parameters] 0 \$::_temp_max\] \$$argsVar break}\n"
		append fullcmd "set $argsVar \[lrange \$args $len end\]"
	}
	uplevel $fullcmd
}

# specialopts can be used to get subcommand specific options (given as e.g. -sniffles-n) in generic interfaces like cg process_project
# key is the first part (e.g. -sniffles). 
# It returns a list of key value ..
proc specialopts {key} {
	global specialopt
	set result {}
	foreach name [array names specialopt *$key-*] {
		if {![regexp -- -?$key\(-.*\) $name temp subkey]} continue
		lappend result $subkey $specialopt($name)
	}
	return $result 
}

proc bgcg_progress {bgexechandleVar args} {
	upvar #0 $bgexechandleVar bgexechandle
	if {![isint $args]} {
		append ::bgerror [lindex $args 0]\n
		return
	}
	if {[catch {
		progress set $args
	}]} {
		puts error
		Extral::bgexec_cancel $bgexechandle
	}
}

proc bgcg {progresscommand channelvar cmd args} {
	# puts "progresscommand cg $args"
	if {[info exists ::stderr_redirect]} {
		set tempfile $::stderr_redirect
	} else {
		set tempfile [tempfile]
	}
	if {[string length $args] < 2000} {
		set ::bgerror {}
		Extral::bgexec -progresscommand [list $progresscommand $channelvar] -no_error_redir -channelvar $channelvar \
				cg $cmd {*}$args 2>@1
		if {$::bgerror ne ""} {error $::bgerror}
	} else {
		set poss [list_concat [list_find -glob $args ">*"] [list_find -glob $args "<*"]]
		if {[llength $poss]} {
			set pos [min $poss]
			set redirect [lrange $args $pos end]
			set code [lrange $args 0 [expr {$pos-1}]]
		} else {
			set code $args
			set redirect {}
		}
		set temprunfile [tempfile]
		file_write $temprunfile [list cg_$cmd {*}$code]\n
		set ::bgerror {}
		Extral::bgexec -progresscommand [list $progresscommand $channelvar] -no_error_redir -channelvar $channelvar \
				cg source $temprunfile {*}$redirect 2>@1
		if {$::bgerror ne ""} {error $::bgerror}
	}
}

proc cmd_getoptions {cmd} {
	if {![catch $cmd temp]} {
		error "Cannot get available options for \"$args\", does not give an error"
	}
	if {![regexp {with options: *([^\n]*)} $temp temp methodoptions]} {
		error "Cannot get available options for \"$args\", error does not contain \"with options: \" followed by available options"
	}
	split $methodoptions ,
}
