#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg {cmd args} {
	# puts "cg $args"
	if {[info exists ::stderr_redirect] && [lsearch -regexp $args 2>.?] == -1} {
		set tempfile $::stderr_redirect
		set redir [list 2> $::stderr_redirect]
	} else {
		set tempfile [tempfile]
		set redir ""
	}
	if {[string length $args] < 2000} {
		set error [catch {exec cg $cmd {*}$args {*}$redir} result]
	} else {
		set temprunfile [tempfile]
		set poss [list_concat [list_find -glob $args ">*"] [list_find -glob $args "<*"]]
		if {[llength $poss]} {
			set pos [min $poss]
			set redirect [lrange $args $pos end]
			set code [lrange $args 0 [expr {$pos-1}]]
		} else {
			set code $args
			set redirect {}
		}
		file_write $temprunfile [list cg_$cmd {*}$code]\n
		set error [catch {exec cg source $temprunfile {*}$redirect {*}$redir} result]
	}
	if {$error} {
		if {[file exists $tempfile]} {
			set errmessage [file_read $tempfile]\n
		} else {
			set errmessage {}
		}
		if {[info exists ::stderr_redirect]} {file delete $tempfile}
		return -code error $errmessage$result
	} else {
		if {[info exists ::stderr_redirect]} {file delete $tempfile}
		return $result
	}
}

proc cg_options {cmd argsVar def {parameters {}} {minargs {}} {maxargs ...} {summary {}}} {
# putsvars cmd argsVar def parameters minargs maxargs summary
	set options [join [list_unmerge $def] ,]
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
				error "unknown option \\"\$key\\", must be one of: $options"
			}
			break
		}]
	}
	set fullcmd [subst {
		set pos 0
		while 1 {
			set key \[lindex \$$argsVar \$pos\]
			incr pos
			set value \[lindex \$$argsVar \$pos\]
			incr pos
			switch -- \$key {
				$def
				-- {incr pos ; break}
				default \{$default\}
			}
		}
		incr pos -2
		set $argsVar \[lrange \$$argsVar \$pos end\]
		set len \[llength \$$argsVar\]
		if {(\$len < $minargs) $test} {
			[list errorformat $cmd $options $minargs $maxargs $parameters]
		}
	}]
	if {[llength $parameters]} {
		append fullcmd "set ::_temp_max \[expr {\[llength \$$argsVar\] - 1}\]\n"
		append fullcmd "if {\$::_temp_max >= 0} {foreach \[lrange [list $parameters] 0 \$::_temp_max\] \$$argsVar break}\n"
		append fullcmd "set $argsVar \[lrange \$args $len end\]"
	}
	uplevel $fullcmd
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

