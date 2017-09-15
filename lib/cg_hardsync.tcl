proc cg_hardsync {args} {
	if {[llength $args] < 2} {
		error "wrong # of arguments, correct format is: cg hardsync src ?...? dest"
	}
	set opts {}
	set num 0
	foreach el $args {
		if {[string index $el 0] ne "-"} break
		if {[string index $el 1] eq "-"} break
		lappend opts $el
		incr num
	}
	set args [lrange $args $num end]
	set dest [lindex $args end]
	set args [lrange $args 0 end-1]
	if {![file isdir $dest]} {
		error "destination must be an existing directory"
	}
	foreach src $args {
		set src [file_absolute $src]
		set src [string trimright $src /]
		set project [file tail $src]
		set finaldest $dest/$project
		puts "cp -alf $src $dest"
		hardlink -f $src $dest >@ stdout 2>@ stderr
		puts "rsync -av $opts --delete --link-dest=$src $src/ $finaldest"
		exec rsync -av {*}$opts --delete --link-dest=$src $src/ $finaldest >@ stdout 2>@ stderr]
	}
}
