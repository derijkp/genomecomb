proc liftover_refdb {old new dest dbbuild build} {
	set root [gzroot $old]
	set newroot [gzroot $new]
	cg liftover -split 0 $old $new ${dest}/liftover/${dbbuild}To${build}.over.tsv
	if {[file exists $root.info]} {
		set c [file_read $root.info]
		regsub citation $c "liftover\t${dbbuild}to${build}\ncitation" c
		file_write $newroot.info $c
	}
	if {[file exists $root.opt]} {
		file rename $root.opt $newroot.opt
	}
	file delete $old $old.lz4i
	file delete $root.info
}
