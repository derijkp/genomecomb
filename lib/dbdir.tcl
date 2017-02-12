proc dbdir {args} {
	global dbdir
	if {[llength $args]} {
		set temp [lindex $args 0]
		if {$temp ne ""} {
			set dbdir $temp
		}
	}
	if {![info exists dbdir]} {
		if {[info exists ::env(GENOMECOMB_DBDIR)]} {
			set dbdir $::env(GENOMECOMB_DBDIR)
		} else {
			error "dbdir not specified, use options e.g. (-dbdir) or environment variable GENOMECOMB_DBDIR to set"
		}
	}
	if {![file isdir $dbdir]} {
		error "dbdir $dbdir does not exist, use options e.g. (-dbdir) or environment variable GENOMECOMB_DBDIR to change"
	}
	set dbdir [file_absolute $dbdir]
	return $dbdir
}
