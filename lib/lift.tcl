proc openliftoverfile {liftoverfile {headerVar {}} {oldrefVar {}} {newrefVar {}}} {
	if {$headerVar ne ""} {upvar $headerVar header}
	if {$oldrefVar ne ""} {upvar $oldrefVar oldref}
	if {$newrefVar ne ""} {upvar $newrefVar newref}
	if {[file ext $liftoverfile] eq ".chain"} {
		set useliftoverfile [file root $liftoverfile].tsv
		if {![file exists $useliftoverfile]} {
			cg liftchain2tsv $liftoverfile $useliftoverfile
		}
		set liftoverfile $useliftoverfile
	}
	set f [gzopen $liftoverfile]
	set header [tsv_open $f comment]
	if {$header ne {chromosome begin end strand destchromosome destbegin destend deststrand}} {
		exiterror "header of file $liftoverfile should be: chromosome begin end strand destchromosome destbegin destend deststrand"
	}
	set cinfo [comment2dict $comment]
	set oldref [dict get $cinfo ref]
	set newref [dict get $cinfo destref]
	return $f
}

proc closeliftoverfile f {
	gzclose $f
}
