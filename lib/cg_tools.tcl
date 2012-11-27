proc cg_hardsync {args} {
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
		puts stderr "destination must be an existing directory"
		exit 1
	}
	foreach src $args {
		set src [file normalize $src]
		set src [string trimright $src /]
		set project [file tail $src]
		set finaldest $dest/$project
		puts "cp -alf $src $dest"
		exec cp -alf $src $dest >@ stdout 2>@ stderr
		puts "rsync -av $opts --delete --link-dest=$src $src/ $finaldest"
		eval exec rsync -av $opts [list --delete --link-dest=$src $src/ $finaldest >@ stdout 2>@ stderr]
	}
}

proc cg_rsync {args} {
	set list {}
	foreach el $args {
		lappend list [string trimright $el /]
	}
	puts "rsync -av --delete $list"
	eval exec rsync -av --delete $list >@ stdout 2>@ stderr
}