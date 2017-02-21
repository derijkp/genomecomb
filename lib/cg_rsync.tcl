proc cg_rsync {args} {
	set list {}
	if {[llength $args] < 2} {
		error "wrong # of arguments, correct format is: cg rsync src ?...? dest"
	}
	foreach el $args {
		lappend list [string trimright $el /]
	}
	puts "rsync -av --delete $list"
	eval exec rsync -av --delete $list >@ stdout 2>@ stderr
}
