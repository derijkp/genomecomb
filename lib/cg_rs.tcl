proc cg_rs args {
	set list {}
	foreach el $args {
		lappend list [string trimright $el /]
	}
	puts "rsync -av --delete $list"
	exec rsync -avH --delete {*}$list >@ stdout 2>@ stderr
}
