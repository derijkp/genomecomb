proc liftoverfile {liftoverfile} {
	if {[file ext $liftoverfile] eq ".chain"} {
		set useliftoverfile [file root $liftoverfile].tsv
		if {![file exists $useliftoverfile]} {
			cg liftchain2tsv $liftoverfile $useliftoverfile
		}
		return $useliftoverfile
	} else {
		return $liftoverfile
	}
}
