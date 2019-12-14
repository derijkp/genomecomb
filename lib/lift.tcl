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
		error "header of file $liftoverfile should be: chromosome begin end strand destchromosome destbegin destend deststrand"
	}
	set cinfo [comment2dict $comment]
	if {[dict exists $cinfo ref]} {
		set oldref [dict get $cinfo ref]
	} else {
		set liftoverfilebase [lindex [split [file tail $liftoverfile] .] 0]
		if {![regexp {^(.*)(To|2)(.*)} $liftoverfilebase temp oldrefname temp newrefname]} {
			set oldref old
		}
		set oldref [string tolower $oldrefname]
	}
	if {[dict exists $cinfo ref]} {
		set newref [dict get $cinfo destref]
	} else {
		set liftoverfilebase [lindex [split [file tail $liftoverfile] .] 0]
		if {![regexp {^(.*)(To|2)(.*)} $liftoverfilebase temp oldrefname temp newrefname]} {
			set newref new
		}
		set newref [string tolower $newrefname]
	}
	return $f
}

proc closeliftoverfile f {
	gzclose $f
}
