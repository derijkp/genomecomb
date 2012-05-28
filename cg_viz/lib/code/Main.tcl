proc main args {
#putsvars ::appdir
lappend ::auto_path [file normalize $::Classy::appdir/../lib]
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
