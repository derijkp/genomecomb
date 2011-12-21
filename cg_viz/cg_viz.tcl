#!/bin/sh
# the next line restarts using wish \
exec wish "$0" ${1+"$@"}
# ClassyTk Builder v0.2

package require Tk
wm withdraw .
if {"[lindex $argv 0]" == "-builder"} {
	set builder 1
	set argv [lrange $argv 1 end]
}
set script [file normalize [info script]]
if {"$script"==""} {
	set appname classyapp
} else {
	if {"$tcl_platform(platform)"=="unix"} {
		while 1 {
			if {[catch {set script [file normalize [file readlink $script]]}]} break
		}
	}
	set appname [file root [file tail $script]]
}
tk appname $appname
if {[package require ClassyTk] < 1.0} {
	error "version conflict for package \"ClassyTk\": need version 1.0 or later"
}
set pwd [pwd]
lappend Classy::help_path [file join $Classy::appdir help]
lappend auto_path [file join $Classy::appdir lib interface] [file join $Classy::appdir lib code]
set Classy::starterror [catch {eval main $argv} Classy::result]
set Classy::starterrorinfo $errorInfo
if {$Classy::starterror} {
	puts $Classy::result
}
if {[info exists builder]} {
	Classy::Builder .classy__.builder
	raise .classy__.builder
	if {$Classy::starterror} {
		set errorInfo $Classy::starterrorinfo
		bgerror $Classy::result
	}
} elseif {$Classy::starterror} {
	error $Classy::result $Classy::starterrorinfo
}

