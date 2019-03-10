proc cg_zindex {file} {
	set type [string range [file extension $file] 1 end]
	if {[auto_load cg_${type}index]} {
		cg_${type}index $file
	}
}
