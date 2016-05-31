proc plugin_addbutton {type name text command args} {
	set tools [Classy::DynaTool define MainTool]
	append tools "$type \"$name\" \"$text\" [list $command] $args"
	Classy::DynaTool define MainTool $tools
}

proc loadplugins {} {
	foreach file [glob -nocomplain $Classy::appdir/plugins/*.tcl] {
		uplevel #0 source $file
	}
}
