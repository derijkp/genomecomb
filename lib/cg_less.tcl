proc cg_less {args} {
	set file [list_pop args]
	exec {*}[catprog $file] $file | less {*}$args >@ stdout 2>@ stderr
}
