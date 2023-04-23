proc shadow_clean {{shadowbase {}}} {
	global env
	if {$shadowbase ne ""} {
		# shadowbase given directly, do nothing extra
	} elseif {[info exists env(SHADOWBASE)]} {
		set shadowbase $env(SHADOWBASE)
	} else {
		error "shadowbase not given, and could not find env var SHADOWBASE"
	}
	set shadowbase [file_absolute $shadowbase]
	foreach shadowdir [glob $shadowbase/*] {
		set link [file link $shadowdir/shadow_source]
		if {
			[file exists $link] &&
			![catch {file link $link} linklink] &&
			$linklink eq $shadowdir
		} continue
		file delete -force $shadowdir
	}
}

proc cg_shadow_clean {args} {
	set shadowbase {}
	cg_options shadow_mkdir args {
		-shadowbase {
			set shadowbase $value
		}
	} {shadowbase} 0 1 {
		make a shadow dir
	}
	shadow_clean $shadowbase
}

proc shadow_delete {link} {
	if {![file isdir $link]} {
		file delete $link
		return
	}
	if {[catch {file link $link} shadowdir]} {
		file delete -force $link
		return
	}
	if {![file exists $shadowdir/shadow_source]} {
		error "$link links to $shadowdir is not a shadowdir (file $shadowdir/shadow_source does not exist)"
	}
	file delete -force $shadowdir
	file delete $link
}

proc cg_shadow_delete {args} {
	set shadowbase {}
	cg_options shadow_delete args {
	} {link} 1 1 {
		delete a shadow dir
	}
	shadow_delete $link
}

# shadow_mkdir creates a temp or workdir for distributed analysis.
# - it can be accessed from different nodes in cluster
# - it has to be explicitely deleted: it does not get auto-deleted when a job or program finishes
proc shadow_mkdir {link {shadowbase {}}} {
	global env
	if {[file exists $link]} {
		return $link
	}
	if {$shadowbase ne ""} {
		# shadowbase given directly, do nothing extra
	} elseif {[info exists env(SHADOWBASE)]} {
		set shadowbase $env(SHADOWBASE)
	} else {
		# putslog "shadowbase not given, and could not find evn var SHADOWBASE, using $link directly"
		mkdir $link
		return $link
	}
	set shadowbase [file_absolute $shadowbase]
	# find unused dir for shadow and create
	for {set i 0} {$i < 20} {incr i} {
		set shadowdir [file join $shadowbase shadow.[pid]-[Extral::randstring 20]]-[file tail $link]
		if {[file exists $shadowdir]} continue
		if {[catch {
			file mkdir $shadowdir
			if {$::tcl_platform(platform) eq "unix"} {
				file attributes $shadowdir -permissions 0700
			}
			set files [glob -nocomplain $shadowdir/*]
			if {[llength $files]} {
				error "Very fishy: there are files in the temporary directory I just created"
			}
		}]} continue
		break
	}
	# create backlink (for cleaning)
	mklink $link $shadowdir/shadow_source 1
	# create link
	mklink $shadowdir $link 1
	return $link
}

proc cg_shadow_mkdir {args} {
	set shadowbase {}
	cg_options shadow_mkdir args {
		-shadowbase {
			set shadowbase $value
		}
	} {link shadowbase} 1 2 {
		make a shadow dir
	}
	shadow_mkdir $link $shadowbase
}
