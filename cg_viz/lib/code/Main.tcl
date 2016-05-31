proc main args {
	global env appdir cgdir cgbase
	#putsvars ::appdir
	package require genomecomb
	set appdir $::Classy::appdir
	lappend ::auto_path [file join [pwd] $appdir/../lib]
	lappend ::auto_path [file dir $appdir]/lib [file dir $appdir]/lib-exp
	genomecombenv
	mainw .mainw
	focus .mainw
	#Classy::cmd
	if {![llength $args]} {
	} elseif {[file exists [lindex $args 0]]} {
		.mainw opentsv {*}$args
	} else {
		.mainw opendb {*}$args
	}
	loadplugins
}
