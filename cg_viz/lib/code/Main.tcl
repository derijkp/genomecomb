proc main args {
	global env appdir
	#putsvars ::appdir
	set appdir $::Classy::appdir
	lappend ::auto_path [file normalize $appdir/../lib]
	set cgdir [file normalize $appdir/..]
	lappend ::auto_path $cgdir/lib $cgdir/lib-exp
	set env(PATH) $appdir/bin:$appdir/extern:[file dir [file dir $appdir]]/bin:$env(PATH)
	mainw .mainw
	focus .mainw
	#Classy::cmd
	if {![llength $args]} {
	} elseif {[file exists [lindex $args 0]]} {
		.mainw opentsv {*}$args
	} else {
		.mainw opendb {*}$args
	}
}
