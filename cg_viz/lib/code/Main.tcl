proc main args {
#putsvars ::appdir
lappend ::auto_path [file normalize $::Classy::appdir/../lib]
mainw .mainw
focus .mainw
Classy::cmd
.mainw opendb {*}$args
}
