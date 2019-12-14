#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

package require Extral

proc cg_wish {args} {
	package require Tk
	set tk 1
	package require Tclx
	catch {signal -restart error SIGINT}
	if {[info commands "console"] == "console"} {
		console show
	} else {
		package require ClassyTk
		Classy::cmd
	}
}

proc cg_sh {args} {
	global tcl_platform
	if {[lsearch $args tk] != -1} {
		cg_wish {*}$args
		return
	}
	if {[info commands "console"] == "console"} {
		console show
	} elseif {[lsearch $args el] != -1 && ![catch {package require eltclsh}]} {
		namespace eval el {}
		proc el::echo {msg} {return $msg}
		set ::el::prompt1 {el::echo "cg% "}
		set ::el::prompt2 {el::echo ""}
		# change so tab just gives a tab character, ^O will still do completion
		auto_load el::matches
		rename el::matches el::matches.ori
		proc el::matches {string type} {
			if {$type == 1} {
				# tab returns \t
				set len [string length $string]
				return [list $len $len [list \t {} {}]]
			} else {
				# do expansion on ^O
				el::matches.ori $string $type
			}
		}
		uplevel #0 interactive
	} elseif {[lsearch $args nox] == -1 && ![catch {package require Tclx}]} {
		catch {signal -restart error SIGINT}
		rename cindex tclx_cindex
		uplevel #0 {commandloop -prompt1 {puts -nonewline "% "} -prompt2 {puts -nonewline ""}}
	} else {
		if {![catch {package require Tclx}]} {
			catch {signal -restart error SIGINT}
			rename cindex tclx_cindex
		}
		package require TclReadLine
		uplevel #0 TclReadLine::interact
	}
}

proc cg_source {file args} {
	set file [file_absolute $file]
	set ::argv $args
	uplevel #0 [list source $file]
}

proc cg_exec {commands args} {
	set ::argv $args
	uplevel #0 $commands
}
