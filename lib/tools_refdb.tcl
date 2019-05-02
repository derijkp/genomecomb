proc liftover_refdb {old new dest dbbuild build {split 1}} {
	set root [gzroot $old]
	set newroot [gzroot $new]
	set liftoverfile ${dest}/liftover/${dbbuild}To${build}.over.tsv
	if {![file exists $liftoverfile]} {
		set liftoverfile ${dest}/liftover/${dbbuild}To[string toupper [string index $build 0]][string range $build 1 end].over.tsv
	}
	cg liftover -split $split $old $new $liftoverfile
	if {[file exists $root.info]} {
		set c [file_read $root.info]
		regsub citation $c "liftover\t${dbbuild}to${build}\ncitation" c
		file_write $newroot.info $c
	}
	if {[file exists $root.opt]} {
		file rename $root.opt $newroot.opt
	}
	file delete $old $old.zsti
	file delete $root.info
}
