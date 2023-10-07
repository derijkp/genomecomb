proc shadow_clean {{shadowdir {}}} {
	global env
	if {$shadowdir ne ""} {
		# shadowdir given directly, do nothing extra
	} elseif {[info exists env(SHADOWDIR)]} {
		set shadowdir $env(SHADOWDIR)
	} else {
		error "shadowdir not given, and could not find env var SHADOWDIR"
	}
	set shadowdir [file_absolute $shadowdir]
	foreach shadow [glob $shadowdir/*] {
		set link [file link $shadow/shadow_source]
		if {
			[file exists $link] &&
			![catch {file link $link} linklink] &&
			$linklink eq $shadow
		} continue
		file delete -force $shadow
	}
}

proc cg_shadow_clean {args} {
	set shadowdir {}
	cg_options shadow_mkdir args {
		-shadowdir {
			set shadowdir $value
		}
	} {shadowdir} 0 1 {
		make a shadow dir
	}
	shadow_clean $shadowdir
}

proc shadow_delete {link} {
	if {![file isdir $link]} {
		file delete $link
		return
	}
	if {[catch {file link $link} shadow]} {
		file delete -force $link
		return
	}
	if {![file exists $shadow/shadow_source]} {
		error "$link links to $shadow is not a shadow (file $shadow/shadow_source does not exist)"
	}
	file delete -force $shadow
	file delete $link
}

proc cg_shadow_delete {args} {
	cg_options shadow_delete args {
	} {link} 1 1 {
		delete a shadow dir
	}
	shadow_delete $link
}

# shadow_mkdir creates a temp or workdir for distributed analysis.
# - it can be accessed from different nodes in cluster
# - it has to be explicitely deleted: it does not get auto-deleted when a job or program finishes
proc shadow_mkdir {link {shadowdir {}}} {
	global env
	if {[file exists $link]} {
		return $link
	}
	if {$shadowdir ne ""} {
		# shadowdir given directly, do nothing extra
	} elseif {[info exists env(SHADOWDIR)]} {
		set shadowdir $env(SHADOWDIR)
	} else {
		# putslog "shadowdir not given, and could not find evn var SHADOWDIR, using $link directly"
		mkdir $link
		return $link
	}
	set shadowdir [file_absolute $shadowdir]
	# find unused dir for shadow and create
	for {set i 0} {$i < 20} {incr i} {
		set shadow [file join $shadowdir shadow.[pid]-[Extral::randstring 20]]-[file tail $link]
		if {[file exists $shadow]} continue
		if {[catch {
			file mkdir $shadow
			if {$::tcl_platform(platform) eq "unix"} {
				file attributes $shadow -permissions 0700
			}
			set files [glob -nocomplain $shadow/*]
			if {[llength $files]} {
				error "Very fishy: there are files in the temporary directory I just created"
			}
		}]} continue
		break
	}
	# create backlink (for cleaning)
	mklink $link $shadow/shadow_source 1
	# create link
	mklink $shadow $link 1
	return $link
}

proc cg_shadow_mkdir {args} {
	set shadowdir {}
	cg_options shadow_mkdir args {
		-shadowdir {
			set shadowdir $value
		}
	} {link shadowdir} 1 2 {
		make a shadow dir
	}
	shadow_mkdir $link $shadowdir
}
